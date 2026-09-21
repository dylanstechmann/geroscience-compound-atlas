"""PubChem PUG REST API client with local disk caching and NCBI rate-limiting compliance."""

import json
import logging
import time
import urllib.parse
from pathlib import Path
from typing import Any

import httpx

logger = logging.getLogger(__name__)

PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


class PubChemClient:
    """Client for resolving chemical names to PubChem CIDs, SMILES, and InChIKeys."""

    def __init__(
        self,
        cache_path: str | Path = "data/interim/pubchem_cache.json",
        rate_limit_sec: float = 0.25,
        timeout_sec: float = 10.0,
    ) -> None:
        self.cache_path = Path(cache_path)
        self.rate_limit_sec = rate_limit_sec
        self.timeout_sec = timeout_sec
        self._last_request_time: float = 0.0
        self._cache: dict[str, Any] = self._load_cache()
        self._client = httpx.Client(timeout=self.timeout_sec)

    def close(self) -> None:
        """Close the underlying HTTP client."""
        self._client.close()

    def _load_cache(self) -> dict[str, Any]:
        """Load cache from disk if it exists."""
        if self.cache_path.exists():
            try:
                with open(self.cache_path, "r", encoding="utf-8") as f:
                    return json.load(f)
            except (json.JSONDecodeError, OSError) as exc:
                logger.warning("Failed to load PubChem cache from %s: %s", self.cache_path, exc)
        return {}

    def _save_cache(self) -> None:
        """Persist cache to disk."""
        self.cache_path.parent.mkdir(parents=True, exist_ok=True)
        try:
            with open(self.cache_path, "w", encoding="utf-8") as f:
                json.dump(self._cache, f, indent=2)
        except (TypeError, OSError) as exc:
            logger.warning("Failed to save PubChem cache to %s: %s", self.cache_path, exc)

    def _wait_for_rate_limit(self) -> None:
        """Enforce NCBI rate limits (<= 5 requests per second)."""
        elapsed = time.time() - self._last_request_time
        if elapsed < self.rate_limit_sec:
            time.sleep(self.rate_limit_sec - elapsed)
        self._last_request_time = time.time()

    def resolve_name(self, name: str) -> dict[str, Any] | None:
        """Resolve a compound name to structured PubChem properties.

        Args:
            name: Common chemical or trade name.

        Returns:
            Dict containing cid, canonical_smiles, inchikey, title, etc., or None if not found.
        """
        clean_name = name.strip()
        cache_key = clean_name.lower()

        if cache_key in self._cache:
            return self._cache[cache_key]

        encoded_name = urllib.parse.quote(clean_name)
        url = (
            f"{PUBCHEM_BASE_URL}/compound/name/{encoded_name}/property/"
            "Title,IUPACName,InChIKey,SMILES,ConnectivitySMILES,MolecularWeight/JSON"
        )

        headers = {"User-Agent": "GeroscienceAtlas/0.1.0 (academic research)"}
        retries = 3

        for attempt in range(retries):
            self._wait_for_rate_limit()
            try:
                response = self._client.get(url, headers=headers)

                if response.status_code == 200:
                    data = response.json()
                    props = data.get("PropertyTable", {}).get("Properties", [])
                    if props:
                        first_prop = props[0]
                        smiles = first_prop.get("SMILES") or first_prop.get("ConnectivitySMILES")
                        result = {
                            "raw_name": clean_name,
                            "cid": first_prop.get("CID"),
                            "title": first_prop.get("Title", clean_name),
                            "iupac_name": first_prop.get("IUPACName"),
                            "inchikey": first_prop.get("InChIKey"),
                            "canonical_smiles": smiles,
                            "molecular_weight": first_prop.get("MolecularWeight"),
                            "resolved": True,
                        }
                        self._cache[cache_key] = result
                        self._save_cache()
                        return result

                elif response.status_code == 404:
                    # Name not found in PubChem
                    result = {
                        "raw_name": clean_name,
                        "resolved": False,
                        "error": "PubChem compound not found (404)",
                    }
                    self._cache[cache_key] = result
                    self._save_cache()
                    return result

                elif response.status_code in (429, 503):
                    # Rate limit or temporary service issue - exponential backoff
                    wait_time = (attempt + 1) * 1.5
                    logger.warning(
                        "PubChem returned HTTP %d for '%s', retrying in %.1fs...",
                        response.status_code,
                        clean_name,
                        wait_time,
                    )
                    time.sleep(wait_time)
                    continue

                else:
                    logger.warning(
                        "PubChem unexpected status %d for '%s': %s",
                        response.status_code,
                        clean_name,
                        response.text[:200],
                    )
                    break

            except (httpx.HTTPError, OSError) as exc:
                logger.warning(
                    "Network error resolving '%s' (attempt %d/%d): %s",
                    clean_name,
                    attempt + 1,
                    retries,
                    exc,
                )
                time.sleep((attempt + 1) * 1.0)

        # If resolution completely failed after retries, do not hard-cache failure so it can be retried later
        return None

    def batch_resolve(self, names: list[str]) -> dict[str, dict[str, Any] | None]:
        """Resolve a list of chemical names."""
        results = {}
        for name in names:
            results[name] = self.resolve_name(name)
        return results

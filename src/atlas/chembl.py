"""ChEMBL Web Services REST API client with local disk caching and rate-limiting."""

import json
import logging
import time
import urllib.parse
from pathlib import Path
from typing import Any

import httpx

logger = logging.getLogger(__name__)

CHEMBL_BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"


class ChEMBLClient:
    """Client for retrieving ChEMBL molecule IDs, bioactivities, and target records."""

    def __init__(
        self,
        cache_path: str | Path = "data/interim/chembl_cache.json",
        rate_limit_sec: float = 0.3,
        timeout_sec: float = 15.0,
    ) -> None:
        self.cache_path = Path(cache_path)
        self.rate_limit_sec = rate_limit_sec
        self.timeout_sec = timeout_sec
        self._last_request_time: float = 0.0
        self._cache: dict[str, Any] = self._load_cache()

    def _load_cache(self) -> dict[str, Any]:
        """Load cached queries from disk."""
        if self.cache_path.exists():
            try:
                with open(self.cache_path, "r", encoding="utf-8") as f:
                    return json.load(f)
            except (json.JSONDecodeError, OSError) as exc:
                logger.warning("Failed to load ChEMBL cache from %s: %s", self.cache_path, exc)
        return {"molecules": {}, "activities": {}, "targets": {}}

    def _save_cache(self) -> None:
        """Persist cache to disk."""
        self.cache_path.parent.mkdir(parents=True, exist_ok=True)
        try:
            with open(self.cache_path, "w", encoding="utf-8") as f:
                json.dump(self._cache, f, indent=2)
        except (TypeError, OSError) as exc:
            logger.warning("Failed to save ChEMBL cache to %s: %s", self.cache_path, exc)

    def _wait_for_rate_limit(self) -> None:
        """Enforce request rate limit."""
        elapsed = time.time() - self._last_request_time
        if elapsed < self.rate_limit_sec:
            time.sleep(self.rate_limit_sec - elapsed)
        self._last_request_time = time.time()

    def get_chembl_id_by_inchikey(self, inchikey: str) -> str | None:
        """Retrieve ChEMBL molecule ID for a given InChIKey."""
        clean_key = inchikey.strip().upper()
        mol_cache = self._cache.setdefault("molecules", {})

        if clean_key in mol_cache:
            return mol_cache[clean_key]

        encoded_key = urllib.parse.quote(clean_key)
        url = f"{CHEMBL_BASE_URL}/molecule?molecule_structures__standard_inchi_key={encoded_key}&format=json"
        headers = {"User-Agent": "GeroscienceAtlas/0.1.0 (academic research)"}

        for attempt in range(3):
            self._wait_for_rate_limit()
            try:
                with httpx.Client(timeout=self.timeout_sec) as client:
                    response = client.get(url, headers=headers)

                if response.status_code == 200:
                    data = response.json()
                    molecules = data.get("molecules", [])
                    if molecules:
                        chembl_id = molecules[0].get("molecule_chembl_id")
                        mol_cache[clean_key] = chembl_id
                        self._save_cache()
                        return chembl_id
                    else:
                        mol_cache[clean_key] = None
                        self._save_cache()
                        return None
                elif response.status_code == 404:
                    mol_cache[clean_key] = None
                    self._save_cache()
                    return None
                elif response.status_code in (429, 503):
                    time.sleep((attempt + 1) * 2.0)
                    continue
                else:
                    logger.warning("ChEMBL HTTP %d for key %s", response.status_code, clean_key)
                    break
            except (httpx.HTTPError, OSError) as exc:
                logger.warning("Network error querying ChEMBL for %s: %s", clean_key, exc)
                time.sleep((attempt + 1) * 1.5)

        return None

    def get_compound_activities(self, chembl_id: str, limit: int = 50) -> list[dict[str, Any]]:
        """Fetch bioactivity rows for a ChEMBL molecule ID."""
        clean_id = chembl_id.strip().upper()
        act_cache = self._cache.setdefault("activities", {})

        if clean_id in act_cache:
            return act_cache[clean_id]

        url = f"{CHEMBL_BASE_URL}/activity?molecule_chembl_id={clean_id}&limit={limit}&format=json"
        headers = {"User-Agent": "GeroscienceAtlas/0.1.0 (academic research)"}

        for attempt in range(3):
            self._wait_for_rate_limit()
            try:
                with httpx.Client(timeout=self.timeout_sec) as client:
                    response = client.get(url, headers=headers)

                if response.status_code == 200:
                    data = response.json()
                    raw_acts = data.get("activities", [])
                    simplified = []
                    for act in raw_acts:
                        simplified.append(
                            {
                                "activity_id": act.get("activity_id"),
                                "assay_chembl_id": act.get("assay_chembl_id"),
                                "assay_type": act.get("assay_type"),
                                "assay_description": act.get("assay_description"),
                                "target_chembl_id": act.get("target_chembl_id"),
                                "target_pref_name": act.get("target_pref_name"),
                                "standard_type": act.get("standard_type"),
                                "standard_relation": act.get("standard_relation"),
                                "standard_value": act.get("standard_value"),
                                "standard_units": act.get("standard_units"),
                                "pchembl_value": act.get("pchembl_value"),
                                "document_chembl_id": act.get("document_chembl_id"),
                            }
                        )
                    act_cache[clean_id] = simplified
                    self._save_cache()
                    return simplified
                elif response.status_code in (429, 503):
                    time.sleep((attempt + 1) * 2.0)
                    continue
                else:
                    break
            except (httpx.HTTPError, OSError) as exc:
                logger.warning("Network error querying activities for %s: %s", clean_id, exc)
                time.sleep((attempt + 1) * 1.5)

        return []

    def get_target_details(self, target_chembl_id: str) -> dict[str, Any] | None:
        """Fetch target metadata (UniProt accession, gene symbol, target type)."""
        clean_target = target_chembl_id.strip().upper()
        target_cache = self._cache.setdefault("targets", {})

        if clean_target in target_cache:
            return target_cache[clean_target]

        url = f"{CHEMBL_BASE_URL}/target/{clean_target}.json"
        headers = {"User-Agent": "GeroscienceAtlas/0.1.0 (academic research)"}

        for attempt in range(3):
            self._wait_for_rate_limit()
            try:
                with httpx.Client(timeout=self.timeout_sec) as client:
                    response = client.get(url, headers=headers)

                if response.status_code == 200:
                    data = response.json()
                    components = data.get("target_components", [])
                    accession = components[0].get("accession") if components else None
                    record = {
                        "target_chembl_id": clean_target,
                        "pref_name": data.get("pref_name"),
                        "target_type": data.get("target_type"),
                        "organism": data.get("organism"),
                        "uniprot_accession": accession,
                    }
                    target_cache[clean_target] = record
                    self._save_cache()
                    return record
                elif response.status_code == 404:
                    target_cache[clean_target] = None
                    self._save_cache()
                    return None
                elif response.status_code in (429, 503):
                    time.sleep((attempt + 1) * 2.0)
                    continue
                else:
                    break
            except (httpx.HTTPError, OSError) as exc:
                logger.warning("Network error querying target %s: %s", clean_target, exc)
                time.sleep((attempt + 1) * 1.5)

        return None

"""Reconstruct published aqueous equilibria; no human PK or sleep prediction.

All concentrations are mol/L. Species coefficients are cumulative formation
constants: [R_n G_m] = beta_nm * [R]**n * [G]**m. Coefficients describe the
source's NMR medium, not plasma. Chosen total concentrations are illustrations.
"""

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parent


@dataclass(frozen=True)
class Species:
    name: str
    host_count: int
    guest_count: int
    log10_beta: float

    def validate(self):
        counts = (self.host_count, self.guest_count)
        if any(type(n) is not int or n < 0 for n in counts) or sum(counts) < 2:
            raise ValueError("Species must contain at least two integer units")
        if not self.name or not math.isfinite(self.log10_beta):
            raise ValueError("Species need a name and a finite formation constant")

    def concentration(self, host, guest):
        # The study's range does not approach float overflow or underflow.
        if (self.host_count and host == 0) or (self.guest_count and guest == 0):
            return 0.0
        log_c = self.log10_beta
        if self.host_count:
            log_c += self.host_count * math.log10(host)
        if self.guest_count:
            log_c += self.guest_count * math.log10(guest)
        try:
            value = 10**log_c
        except OverflowError as error:
            raise ValueError("Species concentration exceeds numeric range") from error
        if not math.isfinite(value):
            raise ValueError("Species concentration exceeds numeric range")
        return value


def solve(host_total, guest_total, species):
    """Solve two mass balances by nested monotone bisection.

At a fixed free-host concentration, the guest balance increases with free
guest. After solving it, the host balance increases with free host (the
positive susceptibility Schur complement). Brackets are [0, total]. This
works for this ideal, two-component mass-action model, without kinetic claims.
"""
    if any(not math.isfinite(x) or not 0 <= x <= 0.01 for x in (host_total, guest_total)):
        raise ValueError("Illustrative totals must be finite and within 0 to 0.01 M")
    species = tuple(species)
    if len({s.name for s in species}) != len(species):
        raise ValueError("Species names must be unique")
    for s in species:
        s.validate()

    def inventory(host, guest):
        values = {s.name: s.concentration(host, guest) for s in species}
        host_mass = host + sum(s.host_count * values[s.name] for s in species)
        guest_mass = guest + sum(s.guest_count * values[s.name] for s in species)
        return host_mass, guest_mass, values

    def guest_at(host):
        lo, hi = 0.0, guest_total
        for _ in range(90):
            mid = (lo + hi) / 2
            if inventory(host, mid)[1] > guest_total:
                hi = mid
            else:
                lo = mid
        return (lo + hi) / 2

    lo, hi = 0.0, host_total
    for _ in range(90):
        mid = (lo + hi) / 2
        guest = guest_at(mid)
        if inventory(mid, guest)[0] > host_total:
            hi = mid
        else:
            lo = mid
    host = (lo + hi) / 2
    guest = guest_at(host)
    host_mass, guest_mass, values = inventory(host, guest)
    captured = sum(s.guest_count * values[s.name] for s in species if s.host_count)
    # Guest dimers are not assumed inactive; only host capture is counted.
    return {
        "host_total_M": host_total,
        "guest_total_M": guest_total,
        "free_host_M": host,
        "free_guest_monomer_M": guest,
        "species_M": values,
        "reconstructed_host_M": host_mass,
        "reconstructed_guest_M": guest_mass,
        "host_captured_guest_M": captured,
        "guest_not_host_captured_fraction": 1 - captured / guest_total if guest_total else None,
    }


def load_models():
    spec = json.loads((ROOT / "binding_constants.json").read_text())
    return {
        model["id"]: tuple(Species(**s) for s in model["species"])
        for model in spec["models"]
    }


def main():
    models = load_models()
    rows = []
    # Arbitrary chemistry scenarios, unrelated to intake, blood PK or dosing.
    for guest_total in (1e-6, 1e-5, 1e-4):
        for ratio in (0, 0.25, 0.5, 1, 2, 4, 8, 16):
            for model_id, species in models.items():
                result = solve(ratio * guest_total, guest_total, species)
                rows.append({
                    "model_id": model_id,
                    "illustrative_guest_total_M": guest_total,
                    "host_to_guest_total_ratio": ratio,
                    "free_host_M": result["free_host_M"],
                    "free_guest_monomer_M": result["free_guest_monomer_M"],
                    "guest_not_host_captured_fraction": result["guest_not_host_captured_fraction"],
                    "host_mass_residual_M": result["reconstructed_host_M"] - ratio * guest_total,
                    "guest_mass_residual_M": result["reconstructed_guest_M"] - guest_total,
                })
    with (ROOT / "illustrative_equilibria.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({"rows": len(rows), "interpretation": "aqueous equilibrium only"}))


if __name__ == "__main__":
    main()

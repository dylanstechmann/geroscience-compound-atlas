"""Synthetic competition at one class of independent binding sites.

No measured affinities, physiological concentrations, oral PK, brain kinetics
or human sleep predictions are encoded. All values use one arbitrary shared
concentration scale; the output is an equilibrium illustration only.
"""

import json
import math
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Guest:
    name: str
    total: float
    kd: float


def _bound_fraction(free_host, kd):
    if free_host == 0:
        return 0.0
    return 1 / (1 + kd / free_host)


def equilibrate(host_total, guests):
    """Solve H_total = H_free + sum(C_i,total * H_free/(Kd_i + H_free)).

    Every guest competes for the same host sites, with one guest per site.
    All Kd values must be positive. No cooperative or irreversible binding.
    """
    guests = tuple(guests)
    if not math.isfinite(host_total) or host_total < 0:
        raise ValueError("Host total must be finite and nonnegative")
    if len({guest.name for guest in guests}) != len(guests):
        raise ValueError("Guest names must be unique")
    for guest in guests:
        if not guest.name:
            raise ValueError("Guests must have names")
        if not math.isfinite(guest.total) or guest.total < 0:
            raise ValueError("Guest totals must be finite and nonnegative")
        if not math.isfinite(guest.kd) or guest.kd <= 0:
            raise ValueError("Kd must be finite and strictly positive")
    try:
        math.fsum([host_total, *(guest.total for guest in guests)])
    except OverflowError as error:
        raise ValueError("Combined concentration exceeds numeric range") from error

    lower, upper = 0.0, host_total
    for _ in range(2400):
        midpoint = lower + (upper - lower) / 2
        if midpoint == lower or midpoint == upper:
            break
        # At nearly full occupancy, summing total - free separately preserves
        # small free concentrations instead of rounding bound/total to one.
        balance_terms = [midpoint, -host_total]
        for guest in guests:
            if midpoint > guest.kd:
                balance_terms.extend((guest.total,
                                      -guest.total / (1 + midpoint / guest.kd)))
            else:
                balance_terms.append(guest.total * _bound_fraction(midpoint, guest.kd))
        if math.fsum(balance_terms) > 0:
            upper = midpoint
        else:
            lower = midpoint
    free_host = lower + (upper - lower) / 2
    records = []
    for guest in guests:
        free_fraction = 1 / (1 + free_host / guest.kd)
        records.append({"name": guest.name, "total": guest.total, "assumed_kd": guest.kd,
                        "free": guest.total * free_fraction,
                        "bound": guest.total * _bound_fraction(free_host, guest.kd),
                        "free_fraction": free_fraction if guest.total > 0 else None})
    return {"host_total": host_total, "free_host": free_host, "guests": records}


def scenarios():
    target = Guest("target", 1, 0.01)
    return [
        ("target_only", 1, (target,)),
        ("equal_affinity_active_metabolite", 1, (target, Guest("active_metabolite", 1, 0.01))),
        ("poor_metabolite_capture", 1, (target, Guest("active_metabolite", 1, 10))),
        ("abundant_tighter_competitor", 1, (target, Guest("competitor", 20, 0.001))),
        ("abundant_weaker_competitor", 1, (target, Guest("competitor", 20, 1))),
        ("more_capacity_with_competitor", 2, (target, Guest("competitor", 20, 0.001))),
    ]


def main():
    records = [{"scenario": name, **equilibrate(host, guests)}
               for name, host, guests in scenarios()]
    output = {
        "status": "synthetic_equilibrium_examples_no_biological_validation",
        "concentration_scale": "arbitrary_shared_dimensionless_scale",
        "limits": [
            "No measured caffeine, metabolite or physiological competitor affinities.",
            "Target and metabolite labels do not assign real pharmacological activity.",
            "No kinetics, protein binding, compartments, oral absorption or sleep model.",
            "Capacity examples are not dose estimates or candidate-selection scores.",
        ],
        "records": records,
    }
    path = Path(__file__).with_name("synthetic_competitive_binding.json")
    path.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"scenarios": len(records), "biological_validation": "none"}))


if __name__ == "__main__":
    main()

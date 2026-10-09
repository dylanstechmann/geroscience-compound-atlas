"""Dimensionless sensitivity analysis and peptide chemistry bookkeeping.

These are mathematical examples, not fitted pharmacokinetics, candidate scores,
oral bioavailability predictions, exposure targets or human-use instructions.
Run from this directory with Python 3.11+; only the standard library is needed.
"""

import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ATOMIC_WEIGHTS = {"C": 12.011, "H": 1.008, "N": 14.007, "O": 15.999}


def mass(formula):
    return round(sum(ATOMIC_WEIGHTS[element] * number for element, number in formula.items()), 3)


def oral_exposure_proxy(systemic_bioavailability, half_life_relative, free_fraction_relative,
                        cellular_access_relative, activity_relative):
    """Linear, equal-volume approximation proportional to F*t_half*fu*access*activity."""
    values = (systemic_bioavailability, half_life_relative, free_fraction_relative,
              cellular_access_relative, activity_relative)
    if any(not math.isfinite(value) or value < 0 for value in values):
        raise ValueError("Parameters must be finite and nonnegative")
    if systemic_bioavailability > 1 or half_life_relative <= 0:
        raise ValueError("Systemic bioavailability must be <=1 and half-life must be positive")
    return math.prod(values)


def free_guest_fraction(host_relative, kd_relative):
    """One-to-one equilibrium, guest total=1; no kinetics or brain compartment.

    Concentrations and Kd are normalized to the same hypothetical guest total.
    A stable quadratic root enforces binding-site capacity and conservation.
    """
    if any(not math.isfinite(value) or value < 0 for value in (host_relative, kd_relative)):
        raise ValueError("Host and Kd must be finite and nonnegative")
    if host_relative == 0:
        return 1.0
    if kd_relative == 0:
        return max(0.0, 1.0 - host_relative)
    total = 1.0 + host_relative + kd_relative
    # Equivalent discriminant avoids cancellation near H=guest=1 and Kd=0.
    discriminant = (1.0 - host_relative) ** 2 + 2.0 * kd_relative * (1.0 + host_relative) + kd_relative ** 2
    bound = 2.0 * host_relative / (total + math.sqrt(discriminant))
    return max(0.0, 1.0 - bound)


def write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    # Formulae represent neutral covalent peptides, excluding counterions.
    specifications = [
        ("SS31_reference", "D-Arg-Dmt-Lys-Phe-NH2", {"C": 32, "H": 49, "N": 9, "O": 5}, 3,
         "Published reference; Dmt is 2,6-dimethyl-L-tyrosine"),
        ("SS31_N_acetyl_hypothesis", "Ac-D-Arg-Dmt-Lys-Phe-NH2", {"C": 34, "H": 51, "N": 9, "O": 6}, 2,
         "Unvalidated modification hypothesis; caps one positive terminal amine"),
        ("SS31_terminal_D_Phe_hypothesis", "D-Arg-Dmt-Lys-D-Phe-NH2", {"C": 32, "H": 49, "N": 9, "O": 5}, 3,
         "Unvalidated C-terminal stereochemistry comparator; degradation, cardiolipin binding and renal clearance unknown"),
        ("SS_Tyr_hypothesis", "D-Arg-Tyr-Lys-Phe-NH2", {"C": 30, "H": 45, "N": 9, "O": 5}, 3,
         "Published SPN4 membrane/cell comparator; endpoint-specific activity reported, human oral PK and cost equivalence unknown"),
        ("SS20_reference", "Phe-D-Arg-Phe-Lys-NH2", {"C": 30, "H": 45, "N": 9, "O": 4}, 3,
         "Published related scaffold; lacks the phenolic side chain; not interchangeable with SS31"),
        ("Epitalon_reference", "Ala-Glu-Asp-Gly-OH", {"C": 14, "H": 22, "N": 4, "O": 9}, -2,
         "Parent reference; no human rejuvenation inference"),
        ("Epitalon_N_acetyl_hypothesis", "Ac-Ala-Glu-Asp-Gly-OH", {"C": 16, "H": 24, "N": 4, "O": 10}, -3,
         "Modification hypothesis; more negative nominal physiological charge"),
        ("Epitalon_N_acetyl_amidate_hypothesis", "Ac-Ala-Glu-Asp-Gly-NH2", {"C": 16, "H": 25, "N": 5, "O": 9}, -2,
         "Distinct covalent structure from acetylated free acid; unvalidated"),
    ]
    candidates = [{"id": key, "sequence_specification": sequence,
                   "neutral_formula": formula, "neutral_average_mass_da": mass(formula),
                   "nominal_charge_near_physiological_ph": charge, "status": status,
                   "human_oral_bioavailability": "unknown", "improved_human_pk": "not_established"}
                  for key, sequence, formula, charge, status in specifications]
    for candidate in candidates:
        if candidate["id"] == "SS_Tyr_hypothesis":
            candidate["published_alias"] = "SPN4"
            candidate["source_url"] = "https://elifesciences.org/articles/75531"
            candidate["legacy_id_note"] = "Identifier retained; published comparator identified in follow-up source review"
    metadata = {
        "status": "chemistry_bookkeeping_and_unvalidated_hypotheses",
        "limits": ["Mass and nominal charge are not measured permeability, potency or half-life.",
                   "Ionization uses ordinary residue states, not measured compound-specific pKa.",
                   "Modification names are specifications, not certificates of identity or synthesis instructions."],
        "records": candidates,
    }
    (ROOT / "peptide_hypotheses.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")

    scenarios = [
        ("synthetic_reference", 0.01, 1, 1, 1, 1),
        ("absorption_only_improved", 0.1, 1, 1, 1, 1),
        ("half_life_only_improved", 0.01, 4, 1, 1, 1),
        ("longer_lived_but_bound_and_poor_cell_access", 0.01, 4, 0.05, 0.5, 1),
        ("absorbed_but_activity_lost", 0.1, 1, 1, 1, 0.05),
    ]
    exposure = []
    baseline = oral_exposure_proxy(*scenarios[0][1:])
    for name, absorbed, half_life, free, access, activity in scenarios:
        proxy = oral_exposure_proxy(absorbed, half_life, free, access, activity)
        exposure.append({"synthetic_scenario": name, "assumed_systemic_oral_bioavailability": absorbed,
                         "half_life_relative": half_life, "free_fraction_relative": free,
                         "cellular_access_relative": access, "activity_relative": activity,
                         "exposure_proxy_relative_to_reference": proxy / baseline})
    write_csv(ROOT / "synthetic_exposure_sensitivity.csv", exposure)
    binding = [{"assumed_host_to_guest_total_ratio": host, "assumed_kd_to_guest_total_ratio": kd,
                "equilibrium_free_guest_fraction": free_guest_fraction(host, kd)}
               for host in (0, 0.25, 0.5, 1, 1.2, 2)
               for kd in (1, 0.1, 0.01, 0.001, 0.000001)]
    write_csv(ROOT / "synthetic_binding_capacity.csv", binding)
    print(json.dumps({"peptide_records": len(candidates), "synthetic_exposure_scenarios": len(exposure),
                      "synthetic_binding_scenarios": len(binding), "biological_validation": "none"}))


if __name__ == "__main__":
    main()

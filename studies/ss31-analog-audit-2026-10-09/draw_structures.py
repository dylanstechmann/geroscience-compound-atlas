"""2D chemical depictions; no conformational or binding prediction."""

from build_structures import ROOT, SPECIFICATIONS, build_molecule
from rdkit.Chem import Draw, rdDepictor


def draw_group(specifications, name):
    molecules = [build_molecule(specification) for specification in specifications]
    for molecule in molecules:
        rdDepictor.Compute2DCoords(molecule)
    legends = [specification.identifier + "\n"
               + ("Published comparator" if specification.status.startswith("published") else
                  "Unvalidated proposal") for specification in specifications]
    options = Draw.MolDrawOptions()
    options.legendFontSize = 18
    options.baseFontSize = 0.65
    options.padding = 0.08
    image = Draw.MolsToGridImage(molecules, molsPerRow=2, subImgSize=(620, 350),
                                 legends=legends, drawOptions=options)
    image.save(ROOT / name)


def main():
    references = [specification for specification in SPECIFICATIONS
                  if specification.status.startswith("published")]
    proposals = [specification for specification in SPECIFICATIONS
                 if not specification.status.startswith("published")]
    draw_group(references, "published_comparators.png")
    draw_group(proposals, "unvalidated_proposals.png")


if __name__ == "__main__":
    main()

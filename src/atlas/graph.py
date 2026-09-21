"""Knowledge graph structures and in-memory containers for geroscience compounds and edges."""

from .models import Compound, EvidenceEdge, EvidenceGrade, HallmarkSlug

# Canonical ordering of evidence grades for filtering by minimum grade
_GRADE_ORDER = [
    EvidenceGrade.E0,
    EvidenceGrade.E1,
    EvidenceGrade.E2,
    EvidenceGrade.E3,
    EvidenceGrade.E4,
]


class EvidenceGraph:
    """In-memory evidence graph linking compounds, targets, and hallmarks."""

    def __init__(self) -> None:
        self.compounds: dict[str, Compound] = {}
        self.edges: list[EvidenceEdge] = []

    def add_compound(self, compound: Compound) -> None:
        """Add or update a compound in the graph indexed by its InChIKey."""
        self.compounds[compound.inchikey] = compound

    def add_edge(self, edge: EvidenceEdge) -> None:
        """Add an evidence edge to the graph.

        Raises:
            ValueError: If edge.grade is missing or invalid.
        """
        if not edge.grade:
            raise ValueError("Every EvidenceEdge must have an assigned EvidenceGrade.")
        self.edges.append(edge)

    def get_edges_for_compound(
        self, inchikey: str, min_grade: EvidenceGrade | None = None
    ) -> list[EvidenceEdge]:
        """Retrieve edges for a compound, optionally filtered by minimum evidence grade."""
        min_idx = _GRADE_ORDER.index(min_grade) if min_grade else 0

        matching = []
        for edge in self.edges:
            if edge.compound_inchikey == inchikey and _GRADE_ORDER.index(edge.grade) >= min_idx:
                matching.append(edge)
        return matching

    def get_compounds_by_hallmark(
        self, hallmark: HallmarkSlug, min_grade: EvidenceGrade | None = None
    ) -> list[str]:
        """Return distinct compound InChIKeys associated with a hallmark."""
        min_idx = _GRADE_ORDER.index(min_grade) if min_grade else 0

        results = set()
        for edge in self.edges:
            if edge.hallmark == hallmark and _GRADE_ORDER.index(edge.grade) >= min_idx:
                results.add(edge.compound_inchikey)
        return sorted(results)

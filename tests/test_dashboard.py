"""Unit tests for the interactive visualization dashboard and SVG rendering."""

from pathlib import Path

from viz.dashboard import build_dashboard_html, generate_molecule_svg


def test_generate_molecule_svg_valid():
    # Ethanol
    svg = generate_molecule_svg("CCO")
    assert "<svg" in svg
    assert "</svg>" in svg


def test_generate_molecule_svg_none_and_empty():
    svg_none = generate_molecule_svg(None)
    assert "No Structure Available" in svg_none

    svg_empty = generate_molecule_svg("   ")
    assert "No Structure Available" in svg_empty


def test_generate_molecule_svg_invalid():
    svg_invalid = generate_molecule_svg("INVALID_NOT_A_SMILES_999")
    assert "Unparseable Structure" in svg_invalid


def test_build_dashboard_html(tmp_path: Path):
    out_file = tmp_path / "test_dashboard.html"

    # Generate using current repository artifacts
    result_path = build_dashboard_html(output_html=out_file)

    assert result_path.exists()
    assert (
        result_path.stat().st_size > 100_000
    )  # Standalone HTML with SVGs and cards should be substantial

    content = result_path.read_text(encoding="utf-8")

    # Assert core structural sections and tabs exist
    assert "<!DOCTYPE html>" in content
    assert "Geroscience Compound Atlas + Predictive Bench" in content
    assert 'id="tab-cards"' in content
    assert 'id="tab-holes"' in content
    assert 'id="tab-bench"' in content
    assert 'id="tab-errors"' in content

    # Assert specific seed compounds are present
    assert "Rapamycin" in content or "rapamycin" in content
    assert "Metformin" in content or "metformin" in content
    assert "BPC-157" in content
    assert "Epitalon" in content
    assert "FakeGeroCompoundX99" in content

    # Assert benchmark metrics are embedded
    assert "0.9735" in content  # Scaffold baseline AUROC
    assert "0.9434" in content  # Scaffold contender AUROC

    # Assert scientific humility callouts
    assert "QED heuristic" in content
    assert "Coverage Holes" in content

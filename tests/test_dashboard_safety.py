"""Regression checks for safe evidence-card rendering in the standalone dashboard."""

from pathlib import Path

from viz.dashboard import build_dashboard_html


def test_dashboard_escapes_evidence_fields_and_restricts_link_protocols(tmp_path: Path):
    output = build_dashboard_html(output_html=tmp_path / "safe-dashboard.html")
    html = output.read_text(encoding="utf-8")

    assert "function escapeHtml(value)" in html
    assert "function safeHttpUrl(value)" in html
    assert "['http:', 'https:'].includes(url.protocol)" in html
    assert "escapeHtml(e.hallmark" in html
    assert "escapeHtml(e.notes || '')" in html
    assert "Context: ${context.map(escapeHtml).join(' · ')}" in html
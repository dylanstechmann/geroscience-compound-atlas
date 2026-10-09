"""Publication boundary checks for stale snapshots and embedded data."""

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from build_site import ROOT, load_snapshot, publish, render, site_outputs


class PublicationChecks(unittest.TestCase):
    def test_embedded_closing_tag_cannot_escape_json(self):
        catalog = {"text": "</script><script>alert('source')</script>&", "unicode": "β"}
        html = render(catalog, "<script type='application/json'>@@CATALOG@@</script>@@SOURCE_BASE@@")
        payload = html.split("application/json'>", 1)[1].split("</script>", 1)[0]
        self.assertEqual(json.loads(payload), catalog)
        self.assertEqual(html.count("</script>"), 1)
        self.assertNotIn("<script>alert", html)

    def test_stale_catalog_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in ("input_manifest.json", "programs.json"):
                (root / name).write_text((ROOT / name).read_text())
            (root / "research_catalog.json").write_text('{"stale": true}')
            with self.assertRaisesRegex(ValueError, "Catalog is stale"):
                load_snapshot(root)

    def test_saved_site_matches_catalog_and_template(self):
        publish(check=True)
        outputs = site_outputs()
        self.assertEqual(json.loads(outputs["research-catalog.json"]), load_snapshot())
        self.assertNotIn("@@CATALOG@@", outputs["research.html"])

    def test_failed_validation_leaves_existing_site_untouched(self):
        with tempfile.TemporaryDirectory() as directory:
            site = Path(directory)
            sentinel = site / "research.html"
            sentinel.write_text("existing")
            with patch("build_site.site_outputs", side_effect=ValueError("stale")), self.assertRaisesRegex(ValueError, "stale"):
                publish(site=site)
            self.assertEqual(sentinel.read_text(), "existing")
            self.assertFalse((site / "research-catalog.json").exists())


if __name__ == "__main__":
    unittest.main()

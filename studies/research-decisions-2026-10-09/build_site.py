"""Publish a deterministic, qualified research snapshot, rejecting stale catalogs."""

import argparse
import json

from build_catalog import REPO, ROOT, build

SOURCE_COMMIT = "47d1b94fbaf76c43d2faa77d16bc168debc2c9f0"
SOURCE_BASE = "https://github.com/dylanstechmann/geroscience-compound-atlas/blob/" + SOURCE_COMMIT + "/"


def load_snapshot(root=ROOT):
    expected = build(json.loads((root / "input_manifest.json").read_text()),
                     json.loads((root / "programs.json").read_text()))
    catalog = json.loads((root / "research_catalog.json").read_text())
    if catalog != expected:
        raise ValueError("Catalog is stale: review inputs and rebuild before publishing")
    return catalog


def render(catalog, template):
    # Prevent even a source string containing a closing script tag from ending
    # the application/json block. Browser rendering uses textContent throughout.
    payload = json.dumps(catalog, ensure_ascii=True, allow_nan=False)
    for character, escape in (("&", "\\u0026"), ("<", "\\u003c"), (">", "\\u003e")):
        payload = payload.replace(character, escape)
    if template.count("@@CATALOG@@") != 1 or template.count("@@SOURCE_BASE@@") != 1:
        raise ValueError("Template must contain exactly one of each data placeholder")
    return template.replace("@@SOURCE_BASE@@", SOURCE_BASE).replace("@@CATALOG@@", payload)


def site_outputs(root=ROOT):
    catalog = load_snapshot(root)
    return {
        "research.html": render(catalog, (root / "research_page.html").read_text()),
        "research-catalog.json": json.dumps(catalog, indent=2, allow_nan=False) + "\n",
    }


def publish(check=False, site=REPO / "site", root=ROOT):
    outputs = site_outputs(root)
    if check:
        for name, expected in outputs.items():
            target = site / name
            if not target.is_file() or target.read_text() != expected:
                raise ValueError(f"Published output missing or stale: {name}")
    else:
        site.mkdir(parents=True, exist_ok=True)
        for name, content in outputs.items():
            (site / name).write_text(content, encoding="utf-8")
    return list(outputs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Check saved outputs without rewriting them")
    args = parser.parse_args()
    print(json.dumps({"checked" if args.check else "written": publish(check=args.check)}))

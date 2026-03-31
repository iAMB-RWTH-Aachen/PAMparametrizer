#!/usr/bin/env python3
"""
Generate Markdown API stubs (mkdocstrings directives) for all modules under a package.

Writes docs/api/<module>.md with:
  # Module XYZ
  ::: full.import.path.to.module
    options: ...
"""

from __future__ import annotations
import sys
import pkgutil
import importlib
import argparse
import pathlib
import yaml  # pip install pyyaml
from typing import List, Iterable

# --------- USER CONFIG ----------
SRC_ROOT = pathlib.Path("src")
DOCS_API_DIR = pathlib.Path("docs") / "api"
MKDOCS_YML = pathlib.Path("mkdocs.yml")
PACKAGE_PREFIX = "PAMparametrizer"
EXCLUDE_PREFIXES = ("test", "tests", "_")
MKDOCSTRINGS_OPTIONS = {
    "members": True,
    "show_source": True,
    "selection": {"members": True}
}
# --------------------------------

def find_package_modules(package_root: str) -> List[str]:
    """
    Walk the package prefix and return dotted module names.

    Requires the package to be importable (pip install -e . or PYTHONPATH includes src/).
    """
    try:
        pkg = importlib.import_module(package_root)
    except Exception as e:
        raise RuntimeError(f"Could not import package '{package_root}': {e}")

    package_path = pathlib.Path(pkg.__file__).resolve().parent
    modules = []

    for finder, name, ispkg in pkgutil.walk_packages([str(package_path)], prefix=package_root + "."):
        # name is the dotted module name
        # filter private/test modules
        short = name.split(".")[-1]
        if short.startswith(EXCLUDE_PREFIXES):
            continue
        modules.append(name)
    # Put package __init__ first (if present)
    modules.sort()
    if package_root + ".__init__" in modules:
        modules.remove(package_root + ".__init__")
        modules.insert(0, package_root)
    return modules

def ensure_api_dir(path: pathlib.Path):
    path.mkdir(parents=True, exist_ok=True)

def mkdocstrings_block(module_name: str, options: dict | None = None) -> str:
    """
    Create a mkdocstrings directive block for a module.
    """
    if options is None:
        options = MKDOCSTRINGS_OPTIONS
    # render YAML-ish options under the directive (mkdocstrings supports nested options)
    # We want something like:
    # ::: my.module
    #     options:
    #       members: true
    #       show_source: true
    opt_lines = []
    for k, v in options.items():
        if isinstance(v, dict):
            opt_lines.append(f"      {k}:")
            for kk, vv in v.items():
                opt_lines.append(f"        {kk}: {yaml.safe_dump(vv, default_flow_style=True).strip()}")
        else:
            opt_lines.append(f"      {k}: {str(v).lower() if isinstance(v, bool) else v}")
    opt_block = "\n".join(opt_lines)
    block = f"""::: {module_name}
    options:
{opt_block}
"""
    return block

def write_module_md(module_name: str, docs_api_dir: pathlib.Path):
    """Create docs/api/<module_name>.md (module dots replaced by /)"""
    rel_path = module_name.replace(".", "/") + ".md"  # nested folders for packages
    out_path = docs_api_dir / rel_path
    out_path.parent.mkdir(parents=True, exist_ok=True)
    title = module_name
    content = f"# {title}\n\n" + mkdocstrings_block(module_name)
    out_path.write_text(content, encoding="utf8")
    print(f"WROTE {out_path}")

def update_mkdocs_nav(modules: Iterable[str], mkdocs_yml_path: pathlib.Path, api_dir_name: str = "api"):
    """
    Insert or update an "API" section in mkdocs.yml that points to docs/api/* files.

    This attempts to be minimally invasive: if `nav` exists it appends / replaces an "API" entry.
    BACKUP mkdocs.yml before running this function.
    """
    if not mkdocs_yml_path.exists():
        print(f"No {mkdocs_yml_path} found; skipping mkdocs nav update.")
        return

    cfg = yaml.safe_load(mkdocs_yml_path.read_text(encoding="utf8"))
    nav = cfg.get("nav", [])

    # Build API nav entries: top-level grouping "API" -> list of pages
    api_entries = []
    for mod in modules:
        # create friendly path: module -> api/<module...>.md
        label = mod.replace(f"{PACKAGE_PREFIX}.", "")
        md_path = f"{api_dir_name}/{mod.replace('.', '/')}.md"
        api_entries.append({label: md_path})

    # Replace or append API section
    new_nav = []
    replaced = False
    for entry in nav:
        if isinstance(entry, dict) and "API" in entry:
            new_nav.append({"API": api_entries})
            replaced = True
        else:
            new_nav.append(entry)
    if not replaced:
        new_nav.append({"API": api_entries})
    cfg["nav"] = new_nav

    # backup then write
    backup = mkdocs_yml_path.with_suffix(".yml.bak")
    mkdocs_yml_path.rename(backup)
    mkdocs_yml_path.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf8")
    print(f"mkdocs.yml updated (backup written to {backup})")

def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--package", default=PACKAGE_PREFIX, help="Top-level package to document (importable)")
    p.add_argument("--src-root", default=str(SRC_ROOT), help="Source root folder (for resolving packages)")
    p.add_argument("--out", default=str(DOCS_API_DIR), help="Output docs/api folder")
    p.add_argument("--update-mkdocs", action="store_true", help="Update mkdocs.yml navigation")
    args = p.parse_args(argv)

    # Make sure src is on path for import
    src_root = pathlib.Path(args.src_root)
    if str(src_root) not in sys.path:
        sys.path.insert(0, str(src_root))

    out_dir = pathlib.Path(args.out)
    ensure_api_dir(out_dir)

    modules = find_package_modules(args.package)
    print(f"Found {len(modules)} modules in package {args.package}")

    for m in modules:
        write_module_md(m, out_dir)

    if args.update_mkdocs:
        update_mkdocs_nav(modules, MKDOCS_YML)

if __name__ == "__main__":
    main()
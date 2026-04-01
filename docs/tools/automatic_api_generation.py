#!/usr/bin/env python3
"""
Generate Sphinx API rst files (one per top‑level sub‑package) from a Python
package and ensure that index.rst contains a reference to the generated
api.rst.

The script:
  1. walks the package,
  2. creates a markdown stub (mkdocstrings) for each module (optional),
  3. groups those markdown files by the first sub‑package after the top‑level
     package,
  4. writes one <folder>.rst file that contains a toctree of the markdown files,
  5. writes an api.rst that lists all <folder>.rst files,
  6. updates index.rst so that the api.rst entry appears in the main navigation.
"""

from __future__ import annotations
import sys
import pkgutil
import importlib
import argparse
import pathlib
import yaml  # pip install pyyaml
from typing import List, Dict

# ----------------------------------------------------------------------
# USER SETTINGS – adapt to your repository layout
# ----------------------------------------------------------------------
SRC_ROOT = pathlib.Path("Modules")  # source code root
DOCS_SRC = pathlib.Path("docs") / "readthedocs"  # Sphinx source dir
DOCS_API_DIR = DOCS_SRC / "api"  # where we write the generated rst files
PACKAGE_PREFIX = "PAMparametrizer"  # top‑level importable package
EXCLUDE_PREFIXES = ("test", "tests", "_")  # modules we ignore

# mkdocstrings options (used only if you request the markdown stubs)
MKDOCSTRINGS_OPTIONS = {
    "members": True,
    "show_source": True,
    "selection": {"members": True}
}


# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------
def find_package_modules(package_root: str) -> List[str]:
    """Return a list of dotted module names inside *package_root*."""
    try:
        pkg = importlib.import_module(package_root)
    except Exception as exc:
        raise RuntimeError(f"Could not import package '{package_root}': {exc}")

    pkg_path = pathlib.Path(pkg.__file__).resolve().parent
    modules: List[str] = []

    for _, name, _ in pkgutil.walk_packages([str(pkg_path)], prefix=package_root + "."):
        short = name.split(".")[-1]
        if short.startswith(EXCLUDE_PREFIXES):
            continue
        modules.append(name)

    modules.sort()
    init_name = package_root + ".__init__"
    if init_name in modules:
        modules.remove(init_name)
        modules.insert(0, package_root)
    return modules


def ensure_dir(path: pathlib.Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def mkdocstrings_block(module_name: str, options: dict | None = None) -> str:
    """Render a mkdocstrings directive block for *module_name*."""
    if options is None:
        options = MKDOCSTRINGS_OPTIONS

    opt_lines: List[str] = []
    for k, v in options.items():
        if isinstance(v, dict):
            opt_lines.append(f"      {k}:")
            for kk, vv in v.items():
                yaml_val = yaml.safe_dump(vv, default_flow_style=True).strip()
                opt_lines.append(f"        {kk}: {yaml_val}")
        else:
            val = str(v).lower() if isinstance(v, bool) else v
            opt_lines.append(f"      {k}: {val}")
    options_block = "\n".join(opt_lines)

    # Now the f‑string contains only plain variables – no back‑slashes.
    return f"""::: {module_name}
        options:
    {options_block}
    """

def write_api_rst(
    folder_names: Iterable[str],
    api_dir: pathlib.Path,
    docs_src: pathlib.Path,
) -> pathlib.Path:
    """
    Create the top‑level ``api.rst`` file.

    *folder_names* – names of the generated ``<folder>.rst`` files (e.g.
    ``PAM_parametrizer``).
    The toctree entries are written as paths **relative to the Sphinx source
    directory** (``docs_src``), i.e. ``api/PAM_parametrizer.rst``.
    """
    # ------------------------------------------------------------------
    # Build the list of paths that Sphinx will follow.
    # ------------------------------------------------------------------
    rel_paths = [f"api/{name}.rst" for name in sorted(folder_names)]

    # ------------------------------------------------------------------
    # Compose the toctree block.
    # ------------------------------------------------------------------
    toctree = [
        ".. toctree::",
        "   :maxdepth: 2",
        "   :caption: API reference",
        "",
    ]
    toctree.extend([f"   {p}" for p in rel_paths])
    toctree.append("")   # required trailing blank line

    # ------------------------------------------------------------------
    # Full file content (title + toctree)
    # ------------------------------------------------------------------
    title = "API reference"
    content = f"{title}\n{'=' * len(title)}\n\n" + "\n".join(toctree)

    out_path = docs_src / "api.rst"
    out_path.write_text(content, encoding="utf8")
    print(f"WROTE {out_path}")
    return out_path


def write_module_md(module_name: str, docs_api_dir: pathlib.Path) -> pathlib.Path:
    """Create docs/api/<module>.md containing a mkdocstrings block."""
    out_path = docs_api_dir / (module_name.replace(".", "/") + ".md")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    content = f"# {module_name}\n\n" + mkdocstrings_block(module_name)
    out_path.write_text(content, encoding="utf8")
    print(f"WROTE {out_path}")
    return out_path


def write_folder_rst(
    folder_name: str,
    md_files: List[pathlib.Path],
    api_dir: pathlib.Path,          # <-- the folder that contains the .rst files
) -> pathlib.Path:
    """
    Create ``api/<folder_name>.rst`` that contains a toctree of the supplied
    markdown files.  The paths inside the toctree are **relative to the folder
    rst file itself** (i.e. relative to ``api_dir``), because the file lives
    inside that directory.
    """
    # ---------- sanity checks ----------
    for p in md_files:
        if not p.is_file():
            raise FileNotFoundError(f"Expected markdown file not found: {p}")

    # ---------- compute paths **relative to api_dir** ----------
    # Example:
    #   p = docs/readthedocs/api/PAM_parametrizer/KcatConstraintConfig.md
    #   api_dir = docs/readthedocs/api
    #   result -> "PAM_parametrizer/KcatConstraintConfig.md"
    rel_paths = [p.relative_to(api_dir).as_posix() for p in sorted(md_files)]

    # ---------- build the toctree ----------
    toctree = [
        ".. toctree::",
        "   :maxdepth: 2",
        f"   :caption: {folder_name}",
        "",
    ]
    toctree.extend([f"   {p}" for p in rel_paths])
    toctree.append("")   # final blank line required by reST

    # ---------- full file content ----------
    content = f"{folder_name}\n{'=' * len(folder_name)}\n\n" + "\n".join(toctree)

    out_path = api_dir / f"{folder_name}.rst"
    out_path.write_text(content, encoding="utf8")
    print(f"WROTE {out_path}")
    return out_path

def update_index_rst(index_path: pathlib.Path, api_rst_name: str = "api.rst") -> None:
    """
    Ensure that *api.rst* appears inside the first ``.. toctree::`` block of
    ``index.rst``.  The function is idempotent – running it repeatedly will not
    duplicate the line.
    """
    lines = index_path.read_text(encoding="utf8").splitlines()
    new_lines: List[str] = []
    in_toctree = False
    added = False

    for line in lines:
        stripped = line.strip()
        if stripped.startswith(".. toctree::"):
            in_toctree = True
            new_lines.append(line)
            continue

        if in_toctree:
            if stripped.startswith(":"):
                new_lines.append(line)
                continue
            # if not stripped:  # blank line → end of block
            #     if not added:
            #         new_lines.append(f"   {api_rst_name}")
            #         added = True
            #     new_lines.append(line)
            #     in_toctree = False
            #     continue

            if stripped == api_rst_name:
                added = True
            new_lines.append(line)
            continue

        new_lines.append(line)

    if not added:
        new_lines.append("")
        new_lines.append(".. toctree::")
        new_lines.append("   :maxdepth: 1")
        new_lines.append(f"   {api_rst_name}")
        new_lines.append("")

    index_path.write_text("\n".join(new_lines) + "\n", encoding="utf8")
    print(f"Updated {index_path} – ensured entry for {api_rst_name}")


def write_all_md_stubs(
        package_root: str,
        src_root: pathlib.Path,
        output_dir: pathlib.Path,
) -> List[pathlib.Path]:
    """
    Walk ``src_root`` (the directory that contains the Python package) and,
    for every importable module that lives under ``package_root``, write a
    Markdown file with a mkdocstrings directive to ``output_dir``.

    The output path mirrors the dotted module name:
        PAMparametrizer.foo.bar  →  output_dir / "PAMparametrizer/foo/bar.md"

    Returns a list with the absolute Paths of the files that have been written.
    """
    # Ensure the output folder exists
    ensure_dir(output_dir)

    written_files: List[pathlib.Path] = []

    # ``pkgutil.iter_modules`` gives us a clean list of *files* that are importable.
    # We use ``walk_packages`` because it also traverses sub‑packages.
    try:
        pkg = importlib.import_module(package_root)
    except Exception as exc:
        raise RuntimeError(f"Cannot import package '{package_root}': {exc}")

    package_path = pathlib.Path(pkg.__file__).resolve().parent

    for _, mod_name, is_pkg in pkgutil.walk_packages([str(package_path)], prefix=package_root + "."):
        # Skip test / private modules – same filter as before
        short_name = mod_name.split(".")[-1]
        if short_name.startswith(EXCLUDE_PREFIXES):
            continue

        # Build the destination markdown file (dots → slashes)
        rel_md_path = pathlib.Path(*mod_name.split(".")).with_suffix(".md")
        out_path = output_dir / rel_md_path
        out_path.parent.mkdir(parents=True, exist_ok=True)

        # Content – title + mkdocstrings block
        title = f"# {mod_name}\n\n"
        block = mkdocstrings_block(mod_name)  # re‑uses the helper we already have
        out_path.write_text(title + block, encoding="utf8")
        written_files.append(out_path)

    print(f"Generated {len(written_files)} mkdocstrings markdown stubs under {output_dir}")
    return written_files

def write_api_index(folders: Iterable[str], docs_src: pathlib.Path) -> pathlib.Path:
    """
    Create a top‑level ``api.rst`` that contains a toctree pointing to each
    ``api/<folder>.rst`` file generated earlier.
    """
    lines = [
        "API reference",
        "==============",
        "",
        ".. toctree::",
        "   :maxdepth: 2",
        "   :caption: API reference",
        "",
    ]
    for name in sorted(folders):
        lines.append(f"   api/{name}.rst")
    lines.append("")

    out_path = docs_src / "api.rst"
    out_path.write_text("\n".join(lines), encoding="utf8")
    print(f"WROTE {out_path}")
    return out_path
# ----------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------
def main(argv: List[str] | None = None) -> None:
    # --------------------------------------------------------------
    # Argument parsing (unchanged)
    # --------------------------------------------------------------
    parser = argparse.ArgumentParser(
        description="Generate Sphinx API .rst files from a Python package."
    )
    parser.add_argument("--package", default=PACKAGE_PREFIX,
                        help="Top‑level package to document (must be importable).")
    parser.add_argument("--src-root", default=str(SRC_ROOT),
                        help="Root directory that contains the package source code.")
    parser.add_argument("--docs-src", default=str(DOCS_SRC),
                        help="Sphinx source directory (contains index.rst).")
    parser.add_argument("--out", default=str(DOCS_API_DIR),
                        help="Output directory for generated API rst files (inside docs‑src).")
    parser.add_argument("--write-modules-md", action="store_true",
                        help="Also emit a Markdown stub for every Python module (mkdocstrings).")
    args = parser.parse_args(argv)

    # --------------------------------------------------------------
    # Resolve everything to absolute paths
    # --------------------------------------------------------------
    src_root = pathlib.Path(args.src_root).resolve()
    if str(src_root) not in sys.path:
        sys.path.insert(0, str(src_root))

    docs_src = pathlib.Path(args.docs_src).resolve()
    api_dir   = pathlib.Path(args.out).resolve()
    ensure_dir(api_dir)

    # --------------------------------------------------------------
    # 0️⃣  (optional) generate **all** markdown stubs for every module
    # --------------------------------------------------------------
    # The generated files will be placed under docs_src / "api"
    # (which is exactly ``api_dir`` – we already created it above).
    if args.write_modules_md:
        # This will also return the list of paths; we do not need it later,
        # but you could keep it if you want to do extra checks.
        _ = write_all_md_stubs(args.package, src_root, api_dir)

    # --------------------------------------------------------------
    # 1️⃣  Discover modules (used for building the folder → md map)
    # --------------------------------------------------------------
    modules = find_package_modules(args.package)
    print(f"Found {len(modules)} modules in package {args.package}")

    # --------------------------------------------------------------
    # 2️⃣  Build folder → markdown map (same as before)
    # --------------------------------------------------------------
    md_by_folder: Dict[str, List[pathlib.Path]] = {}
    for mod in modules:
        # If we already generated the stub above, the file already exists.
        # Otherwise we create an empty placeholder – this is harmless.
        md_path = api_dir / (mod.replace(".", "/") + ".md")
        md_by_folder.setdefault(mod.split(".")[1] if "." in mod else mod, []).append(md_path)

    # --------------------------------------------------------------
    # 3️⃣  Write one .rst per top‑level folder (unchanged)
    # --------------------------------------------------------------
    generated_folders: List[str] = []
    for folder, md_files in md_by_folder.items():
        write_folder_rst(folder, md_files, api_dir)   # note: now only 3 args
        generated_folders.append(folder)

    # --------------------------------------------------------------
    # 4️⃣  Write the central api.rst that links the folder pages
    # --------------------------------------------------------------
    write_api_index(generated_folders, docs_src)

    # --------------------------------------------------------------
    # 5️⃣  Ensure index.rst contains a reference to api.rst
    # --------------------------------------------------------------
    index_path = docs_src / "index.rst"
    if not index_path.is_file():
        raise FileNotFoundError(f"index.rst not found at {index_path}")
    update_index_rst(index_path, api_rst_name="api.rst")

    print("\n✅  Generation finished – run `make html` again.\n")

if __name__ == "__main__":
    main()
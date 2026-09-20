from pathlib import Path
import subprocess
import sys
try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "python"))

project = "disprs"
release = tomllib.loads((ROOT / "Cargo.toml").read_text())["package"]["version"]
subprocess.run(
    ["cargo", "run", "--locked", "--quiet", "--manifest-path", str(ROOT / "cli" / "Cargo.toml"),
     "--", "help", str(ROOT / "docs" / "ref" / "cli" / "_generated"), "--format", "rst"],
    check=True,
)
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.githubpages",
    "myst_parser",
    "sphinxcontrib_rust",
    "sphinx_fortran_domain",
    "breathe",
]
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
myst_enable_extensions = ["colon_fence"]
myst_ref_domains = ["rust"]
myst_heading_anchors = 2
exclude_patterns = ["_build"]
html_theme = "sphinx_book_theme"
html_title = f"disprs {release}"

rust_crates = {"disprs": str(ROOT)}
rust_doc_dir = str(ROOT / "docs" / "ref" / "api" / "_generated" / "rust")
rust_rustdoc_fmt = "md"
rust_visibility = "pub"

fortran_sources = [str(ROOT / "fortran" / "*_ext.f90")]
fortran_lexer = "ford"

breathe_projects_source = {
    "disprs": (str(ROOT / "include"), ["disprs.h", "dftd3.h", "dftd4.h"])
}
breathe_default_project = "disprs"
breathe_domain_by_extension = {"h": "c"}
breathe_doxygen_config_options = {
    "EXTRACT_ALL": "YES",
    "ENABLE_PREPROCESSING": "YES",
    "MACRO_EXPANSION": "YES",
    "EXPAND_ONLY_PREDEF": "YES",
    "EXPAND_AS_DEFINED": "DFTD3_DELETE_ADAPTER DFTD4_DELETE_ADAPTER",
    "EXCLUDE_SYMBOLS": "DFTD3_DELETE_ADAPTER DFTD4_DELETE_ADAPTER",
    "OPTIMIZE_OUTPUT_FOR_C": "YES",
}
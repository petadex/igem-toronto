"""Compile notebook(s) to markdown with every image embedded inline.

`jupyter nbconvert --to markdown --ExtractOutputPreprocessor.enabled=False`
inlines code-cell *output* images (e.g. a rendered plot) as base64 data
URIs, but two other image sources still land as separate files:

- markdown-cell attachments (`![](attachment:x.png)`, e.g. issue 51) are
  extracted into a `<name>_files/` folder regardless of that flag.
- images referenced on disk via plain `<img src="../resources/...">` /
  `![]()` tags (e.g. issue 74, 75, 78) are left untouched by nbconvert, so
  the compiled markdown ends up with broken links once moved into
  `compiled/<date>/`.

This script runs nbconvert and then rewrites every such reference
(nbconvert-extracted attachments and externally-referenced files alike) as
an inline base64 data URI, so the compiled `.mdx` is fully self-contained
and no `_files/` folder is needed.

Usage:
    python compile_notebook.py 260324_issue51_clusterquantification
    python compile_notebook.py notebooks/260324_issue51_clusterquantification.ipynb
    python compile_notebook.py --all
"""

import argparse
import base64
import mimetypes
import re
import shutil
import subprocess
import sys
import urllib.parse
from pathlib import Path

NOTEBOOKS_DIR = Path(__file__).resolve().parent
COMPILED_DIR = NOTEBOOKS_DIR / "compiled"

MARKDOWN_IMAGE_RE = re.compile(r"!\[([^\]]*)\]\(([^)]+)\)")
HTML_IMAGE_RE = re.compile(r'(<img\b[^>]*\bsrc=")([^"]+)(")')


def date_prefix(notebook_name: str) -> str:
    return notebook_name.split("_")[0]

def issue_number(notebook_name: str) -> str:
    name = notebook_name.split("_")[1]
    number = re.search(r"\d+", name)
    return number.group(0) if number else name


def run_nbconvert(notebook_path: Path, output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "jupyter",
            "nbconvert",
            "--to",
            "markdown",
            str(notebook_path),
            f"--output-dir={output_dir}",
            "--ExtractOutputPreprocessor.enabled=False",
        ],
        check=True,
    )
    # change the suffix to .mdx to match the expected output for the site
    md_path = output_dir / f"{notebook_path.stem}.md"
    md_path.rename(output_dir / f"{notebook_path.stem}.mdx")
    return output_dir / f"{notebook_path.stem}.mdx"


def is_external_reference(src: str) -> bool:
    if src.startswith(("http://", "https://", "data:", "#")):
        return False
    return True


def inline_images(md_path: Path, notebook_dir: Path) -> None:
    """Rewrite every file-based image reference as a base64 data URI.

    Covers both nbconvert-extracted attachments (already copied into
    `<name>_files/`) and externally-referenced images on disk (still
    pointing at their original, notebook-relative path). Cell-output
    images are already inlined by nbconvert itself and are skipped via
    `is_external_reference`.
    """
    files_dir = md_path.with_name(f"{md_path.stem}_files")
    text = md_path.read_text(encoding="utf-8")
    inlined: dict[str, str] = {}

    def resolve_and_inline(raw_src: str) -> str | None:
        decoded = urllib.parse.unquote(raw_src)
        if raw_src.startswith(f"{files_dir.name}/"):
            source_path = (md_path.parent / decoded).resolve()
        else:
            source_path = (notebook_dir / decoded).resolve()
        if not source_path.is_file():
            return None

        cache_key = str(source_path)
        if cache_key in inlined:
            return inlined[cache_key]

        mime, _ = mimetypes.guess_type(source_path.name)
        mime = mime or "application/octet-stream"
        encoded = base64.b64encode(source_path.read_bytes()).decode("ascii")
        data_uri = f"data:{mime};base64,{encoded}"
        inlined[cache_key] = data_uri
        return data_uri

    def replace_markdown(match: re.Match) -> str:
        alt, src = match.group(1), match.group(2)
        if not is_external_reference(src):
            return match.group(0)
        new_ref = resolve_and_inline(src)
        return f"![{alt}]({new_ref})" if new_ref else match.group(0)

    def replace_html(match: re.Match) -> str:
        prefix, src, suffix = match.group(1), match.group(2), match.group(3)
        if not is_external_reference(src):
            return match.group(0)
        new_ref = resolve_and_inline(src)
        return f"{prefix}{new_ref}{suffix}" if new_ref else match.group(0)

    text = MARKDOWN_IMAGE_RE.sub(replace_markdown, text)
    text = HTML_IMAGE_RE.sub(replace_html, text)
    md_path.write_text(text, encoding="utf-8")

    if files_dir.is_dir():
        shutil.rmtree(files_dir)

    if inlined:
        print(f"  inlined {len(inlined)} image(s) as base64 data URIs")

def add_header(md_path: Path, notebook_path: Path) -> None:
    """Add a header to the compiled markdown file with the notebook name and date."""
    with notebook_path.open("r", encoding="utf-8") as f:
        # get the first line of the notebook, which should be a markdown cell with the title
        first_line = f.readline().strip()
        if first_line.startswith("#"):
            title = first_line.lstrip("#").strip()
        else:
            title = notebook_path.stem
    header = f"""---
    title: {title}
    section: "Dry Lab"
    path: "/dry-lab/issue/{issue_number(notebook_path.stem)}"
    navTitle: "Issue {issue_number(notebook_path.stem)}"
    order: 20
    description: ""
    owners: ["Dry Lab"]
    updated: "2026-05-16"
    status: "draft"
---
    """
    with md_path.open("r+", encoding="utf-8") as f:
        content = f.read()
        f.seek(0)
        f.write(header + "\n" + content)


def compile_notebook(notebook_path: Path) -> None:
    if not notebook_path.is_file():
        raise FileNotFoundError(f"Notebook not found: {notebook_path}")

    output_dir = COMPILED_DIR / date_prefix(notebook_path.stem)
    print(f"Converting {notebook_path.name} -> {output_dir}")
    md_path = run_nbconvert(notebook_path, output_dir)
    inline_images(md_path, notebook_path.parent)
    add_header(md_path, notebook_path)


def resolve_notebook_arg(arg: str) -> Path:
    path = Path(arg)
    if path.suffix != ".ipynb":
        path = path.with_suffix(".ipynb")
    if not path.is_absolute() and not path.is_file():
        path = NOTEBOOKS_DIR / path.name
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "notebooks",
        nargs="*",
        help="Notebook name(s) or path(s), with or without .ipynb extension",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Compile every .ipynb file directly under notebooks/",
    )
    args = parser.parse_args()

    if args.all:
        targets = sorted(NOTEBOOKS_DIR.glob("*.ipynb"))
    elif args.notebooks:
        targets = [resolve_notebook_arg(arg) for arg in args.notebooks]
    else:
        parser.error("provide notebook name(s) or --all")
        return

    for notebook_path in targets:
        try:
            compile_notebook(notebook_path)
        except (FileNotFoundError, subprocess.CalledProcessError) as exc:
            print(f"  FAILED: {exc}", file=sys.stderr)


if __name__ == "__main__":
    main()

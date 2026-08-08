"""Compile notebook(s) to markdown and bundle every referenced image locally.

`jupyter nbconvert --to markdown` only extracts images that are embedded in
the notebook itself (code-cell outputs, markdown-cell attachments) into the
`<name>_files/` folder. Some notebooks instead reference images on disk via
plain `<img src="../resources/...">` / `![]()` tags (e.g. issue 74, 75, 78);
nbconvert leaves those paths untouched, so the compiled markdown ends up
with broken links once moved into `compiled/<date>/`.

This script runs nbconvert and then copies any such externally-referenced
image into the notebook's `_files/` folder, rewriting the reference to match
nbconvert's own convention (`<name>_files/<basename>`).

Usage:
    python compile_notebook.py 260324_issue51_clusterquantification
    python compile_notebook.py notebooks/260324_issue51_clusterquantification.ipynb
    python compile_notebook.py --all
"""

import argparse
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
        ],
        check=True,
    )
    return output_dir / f"{notebook_path.stem}.md"


def is_external_reference(src: str) -> bool:
    if src.startswith(("http://", "https://", "data:", "#")):
        return False
    return True


def bundle_external_images(md_path: Path, notebook_dir: Path) -> None:
    files_dir = md_path.with_name(f"{md_path.stem}_files")
    text = md_path.read_text(encoding="utf-8")
    copied: dict[str, str] = {}

    def resolve_and_copy(raw_src: str) -> str | None:
        # already pointing at the bundled files dir, leave it alone
        if raw_src.startswith(f"{files_dir.name}/"):
            return None

        decoded = urllib.parse.unquote(raw_src)
        source_path = (notebook_dir / decoded).resolve()
        if not source_path.is_file():
            return None

        cache_key = str(source_path)
        if cache_key in copied:
            return copied[cache_key]

        files_dir.mkdir(exist_ok=True)
        dest_name = source_path.name
        dest_path = files_dir / dest_name
        suffix = 1
        while dest_path.exists() and dest_path.resolve() != source_path:
            suffix += 1
            dest_path = files_dir / f"{source_path.stem}_{suffix}{source_path.suffix}"
            dest_name = dest_path.name
        if not dest_path.exists():
            shutil.copyfile(source_path, dest_path)

        new_ref = f"{files_dir.name}/{dest_name}"
        copied[cache_key] = new_ref
        return new_ref

    def replace_markdown(match: re.Match) -> str:
        alt, src = match.group(1), match.group(2)
        if not is_external_reference(src):
            return match.group(0)
        new_ref = resolve_and_copy(src)
        return f"![{alt}]({new_ref})" if new_ref else match.group(0)

    def replace_html(match: re.Match) -> str:
        prefix, src, suffix = match.group(1), match.group(2), match.group(3)
        if not is_external_reference(src):
            return match.group(0)
        new_ref = resolve_and_copy(src)
        return f"{prefix}{new_ref}{suffix}" if new_ref else match.group(0)

    text = MARKDOWN_IMAGE_RE.sub(replace_markdown, text)
    text = HTML_IMAGE_RE.sub(replace_html, text)
    md_path.write_text(text, encoding="utf-8")

    if copied:
        print(f"  bundled {len(copied)} external image(s) into {files_dir.name}/")


def compile_notebook(notebook_path: Path) -> None:
    if not notebook_path.is_file():
        raise FileNotFoundError(f"Notebook not found: {notebook_path}")

    output_dir = COMPILED_DIR / date_prefix(notebook_path.stem)
    print(f"Converting {notebook_path.name} -> {output_dir}")
    md_path = run_nbconvert(notebook_path, output_dir)
    bundle_external_images(md_path, notebook_path.parent)


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

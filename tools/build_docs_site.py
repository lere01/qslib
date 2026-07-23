"""Build the local Markdown book and Rust API reference into one site."""

from pathlib import Path
import os
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
OWNERSHIP_MARKER = ".qslib-generated-site"
DESIGN_TAG = "v0.1.0"


def ensure_design_system() -> Path:
    """Fetch the pinned house design system (github.com/lere01/design)."""
    design = ROOT / "design"
    if not design.is_dir():
        subprocess.run(
            "curl -sSfL https://raw.githubusercontent.com/lere01/design/"
            f"{DESIGN_TAG}/ci/fetch-design.sh | sh -s -- {DESIGN_TAG} design",
            cwd=ROOT,
            shell=True,
            check=True,
        )
    return design


def main() -> int:
    destination = Path(sys.argv[1]) if len(sys.argv) == 2 else ROOT / "target" / "qslib-site"
    destination = destination.resolve()
    if destination == ROOT or destination.parent == destination:
        raise SystemExit("refusing to replace the project root or filesystem root")
    if destination.name != "qslib-site" and not destination.name.startswith("qslib-site-"):
        raise SystemExit(
            "documentation output must be a dedicated qslib-site or qslib-site-* directory"
        )
    marker = destination / OWNERSHIP_MARKER
    if destination.exists():
        if not destination.is_dir():
            raise SystemExit(f"documentation output is not a directory: {destination}")
        if not marker.is_file() or marker.read_text(encoding="utf-8") != "qslib-docs-v1\n":
            raise SystemExit(
                f"refusing to delete an unowned documentation directory: {destination}"
            )
        shutil.rmtree(destination)
    destination.mkdir(parents=True)
    design = ensure_design_system()
    subprocess.run(["mdbook", "build", "-d", str(destination)], cwd=ROOT, check=True)
    # mdBook copies the additional-css files but not the font binaries.
    fonts = destination / "design" / "fonts"
    fonts.mkdir(parents=True, exist_ok=True)
    for woff in (design / "fonts").glob("*.woff2"):
        shutil.copy2(woff, fonts / woff.name)
    (destination / OWNERSHIP_MARKER).write_text("qslib-docs-v1\n", encoding="utf-8")
    adapter = design / "adapters" / "rustdoc"
    env = os.environ.copy()
    env["RUSTDOCFLAGS"] = " ".join(
        part
        for part in (
            env.get("RUSTDOCFLAGS", ""),
            f"--html-in-header {adapter / 'header.html'}",
            f"--html-before-content {adapter / 'before-content.html'}",
            f"--html-after-content {adapter / 'after-content.html'}",
        )
        if part
    )
    subprocess.run(
        ["cargo", "doc", "--workspace", "--no-deps", "--all-features"],
        cwd=ROOT,
        check=True,
        env=env,
    )
    api = destination / "api"
    shutil.copytree(ROOT / "target" / "doc", api)
    if not (api / "qslib" / "index.html").is_file():
        raise SystemExit("Rust API reference was not copied into the combined site")
    print(f"documentation site: {destination}")
    print(f"Rust API reference: {api / 'qslib' / 'index.html'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

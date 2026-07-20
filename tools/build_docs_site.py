"""Build the local Markdown book and Rust API reference into one site."""

from pathlib import Path
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
OWNERSHIP_MARKER = ".qslib-generated-site"


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
    subprocess.run(["mdbook", "build", "-d", str(destination)], cwd=ROOT, check=True)
    (destination / OWNERSHIP_MARKER).write_text("qslib-docs-v1\n", encoding="utf-8")
    subprocess.run(
        ["cargo", "doc", "--workspace", "--no-deps", "--all-features"],
        cwd=ROOT,
        check=True,
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

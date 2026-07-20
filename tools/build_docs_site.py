"""Build the local Markdown book and Rust API reference into one site."""

from pathlib import Path
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    destination = Path(sys.argv[1]) if len(sys.argv) == 2 else ROOT / "target" / "qslib-site"
    destination = destination.resolve()
    if destination.exists():
        shutil.rmtree(destination)
    subprocess.run(["mdbook", "build", "-d", str(destination)], cwd=ROOT, check=True)
    subprocess.run(
        ["cargo", "doc", "--workspace", "--no-deps", "--all-features"],
        cwd=ROOT,
        check=True,
    )
    api = destination / "api"
    shutil.copytree(ROOT / "target" / "doc", api)
    print(f"documentation site: {destination}")
    print(f"Rust API reference: {api / 'qslib' / 'index.html'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

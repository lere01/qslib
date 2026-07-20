"""Check relative Markdown links without requiring a documentation toolchain."""

from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parents[1]
LINK = re.compile(r"\[[^\]]+\]\(([^)]+)\)")


def main() -> int:
    errors = []
    for source in [ROOT / "README.md", *sorted((ROOT / "docs").rglob("*.md"))]:
        for target in LINK.findall(source.read_text(encoding="utf-8")):
            if target.startswith(("http://", "https://", "mailto:", "#")):
                continue
            path = target.split("#", 1)[0].split("?", 1)[0]
            if not path:
                continue
            resolved = (source.parent / path).resolve()
            if not resolved.exists():
                errors.append(f"{source.relative_to(ROOT)}: missing {target}")
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print("markdown links ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

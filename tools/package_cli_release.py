"""Create a self-contained qslib CLI archive for one host platform.

The release workflow builds the executable on the target runner and calls
this module with the repository documents. Keeping archive creation in Python
gives Unix and Windows runners the same path-validation and file-set rules.
"""

from __future__ import annotations

import argparse
import re
import shutil
import tarfile
import tempfile
import zipfile
from pathlib import Path


_VERSION = re.compile(r"^\d+\.\d+\.\d+(?:[-+][0-9A-Za-z.-]+)?$")
_PLATFORM = re.compile(r"^[0-9A-Za-z][0-9A-Za-z._-]*$")


def _validate_component(value: str, pattern: re.Pattern[str], label: str) -> None:
    if not pattern.fullmatch(value):
        raise ValueError(f"invalid {label}: {value!r}")


def build_archive(binary: Path, version: str, platform: str, output: Path) -> Path:
    """Package ``binary`` and release documents, returning the archive path."""

    _validate_component(version, _VERSION, "version")
    _validate_component(platform, _PLATFORM, "platform")
    binary = binary.resolve()
    output = output.resolve()
    if not binary.is_file():
        raise FileNotFoundError(binary)
    output.mkdir(parents=True, exist_ok=True)

    root_name = f"qslib-{version}-{platform}"
    archive = output / (
        f"{root_name}.zip" if platform.lower().startswith("windows") else f"{root_name}.tar.gz"
    )
    with tempfile.TemporaryDirectory(prefix="qslib-cli-", dir=output) as temporary:
        staging = Path(temporary) / root_name
        staging.mkdir()
        executable_name = "qslib.exe" if platform.lower().startswith("windows") else "qslib"
        shutil.copy2(binary, staging / executable_name)
        for document in ("README.md", "LICENSE", "RELEASE_NOTES.md"):
            source = binary.parent / document
            if not source.is_file():
                source = Path.cwd() / document
            if source.is_file():
                shutil.copy2(source, staging / document)

        if archive.suffix == ".zip":
            with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as handle:
                for path in sorted(staging.iterdir()):
                    handle.write(path, f"{root_name}/{path.name}")
        else:
            with tarfile.open(archive, "w:gz") as handle:
                for path in sorted(staging.iterdir()):
                    handle.add(path, arcname=f"{root_name}/{path.name}")
    return archive


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--version", required=True)
    parser.add_argument("--platform", required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    archive = build_archive(
        arguments.binary,
        arguments.version,
        arguments.platform,
        arguments.output,
    )
    print(archive)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Tests for the portable CLI release-archive builder."""

from __future__ import annotations

import io
import tarfile
import tempfile
import unittest
import zipfile
from pathlib import Path

from package_cli_release import build_archive


class CliReleaseArchiveTests(unittest.TestCase):
    def _fixture(self, root: Path, *, windows: bool) -> Path:
        binary = root / ("qslib.exe" if windows else "qslib")
        binary.write_bytes(b"test executable")
        (root / "README.md").write_text("readme\n", encoding="utf-8")
        (root / "LICENSE").write_text("license\n", encoding="utf-8")
        (root / "RELEASE_NOTES.md").write_text("notes\n", encoding="utf-8")
        return binary

    def test_unix_archive_contains_binary_and_release_documents(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            binary = self._fixture(root, windows=False)
            archive = build_archive(binary, "1.0.0", "linux-x86_64", root)

            self.assertEqual(archive.name, "qslib-1.0.0-linux-x86_64.tar.gz")
            with tarfile.open(archive, "r:gz") as handle:
                self.assertEqual(
                    sorted(handle.getnames()),
                    [
                        "qslib-1.0.0-linux-x86_64/LICENSE",
                        "qslib-1.0.0-linux-x86_64/README.md",
                        "qslib-1.0.0-linux-x86_64/RELEASE_NOTES.md",
                        "qslib-1.0.0-linux-x86_64/qslib",
                    ],
                )

    def test_windows_archive_uses_zip_and_exe_name(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            binary = self._fixture(root, windows=True)
            archive = build_archive(binary, "1.0.0", "windows-x86_64", root)

            self.assertEqual(archive.name, "qslib-1.0.0-windows-x86_64.zip")
            with zipfile.ZipFile(archive) as handle:
                self.assertEqual(
                    sorted(handle.namelist()),
                    [
                        "qslib-1.0.0-windows-x86_64/LICENSE",
                        "qslib-1.0.0-windows-x86_64/README.md",
                        "qslib-1.0.0-windows-x86_64/RELEASE_NOTES.md",
                        "qslib-1.0.0-windows-x86_64/qslib.exe",
                    ],
                )

    def test_invalid_version_is_rejected_before_writing(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            binary = self._fixture(root, windows=False)
            with self.assertRaises(ValueError):
                build_archive(binary, "../unsafe", "linux-x86_64", root)


if __name__ == "__main__":
    unittest.main()

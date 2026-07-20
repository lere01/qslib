"""Contract tests for third-party distribution metadata.

These tests intentionally inspect the files that package managers and release
automation consume. They do not publish anything or contact a registry.
"""

from __future__ import annotations

import json
import subprocess
import tomllib
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class DistributionContractTests(unittest.TestCase):
    def test_cargo_distribution_set_has_owner_gated_publish_metadata(self) -> None:
        result = subprocess.run(
            ["cargo", "metadata", "--locked", "--no-deps", "--format-version", "1"],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        packages = {package["name"]: package for package in json.loads(result.stdout)["packages"]}
        expected = {
            "qslib-quantum",
            "qslib-quantum-core",
            "qslib-quantum-exact",
            "qslib-quantum-io",
            "qslib-quantum-sse",
            "qslib-quantum-variational",
            "qslib-quantum-cli",
        }
        self.assertTrue(expected.issubset(packages))
        for name in expected:
            package = packages[name]
            self.assertEqual(package["license"], "Apache-2.0")
            self.assertTrue(package["description"])
            self.assertIn(
                package["repository"],
                {"https://github.com/lere01/qslib", "https://github.com/lere01/qslib.git"},
            )
            self.assertEqual(package["publish"], [])

    def test_python_distribution_has_installable_project_metadata(self) -> None:
        metadata = tomllib.loads(
            (ROOT / "crates" / "qslib-python" / "pyproject.toml").read_text()
        )["project"]

        self.assertEqual(metadata["name"], "qslib-quantum")
        self.assertEqual(metadata["requires-python"], ">=3.12")
        self.assertEqual(metadata["readme"], "PYTHON_README.md")
        self.assertEqual(metadata["license"]["file"], "LICENSE-APACHE")
        self.assertEqual(metadata["urls"]["Repository"], "https://github.com/lere01/qslib")
        self.assertEqual(metadata["urls"]["Documentation"], "https://lere01.github.io/qslib/")
        self.assertIn(
            "License :: OSI Approved :: Apache Software License",
            metadata["classifiers"],
        )
        self.assertIn("Programming Language :: Python :: 3 :: Only", metadata["classifiers"])

    def test_python_binding_name_matches_distribution_contract(self) -> None:
        build = tomllib.loads(
            (ROOT / "crates" / "qslib-python" / "pyproject.toml").read_text()
        )["tool"]["maturin"]
        self.assertEqual(build["module-name"], "qslib_quantum")

    def test_release_workflow_keeps_registry_publication_guarded(self) -> None:
        workflow = (ROOT / ".github" / "workflows" / "release.yml").read_text()
        self.assertIn("publish_pypi", workflow)
        self.assertIn("inputs.publish", workflow)
        self.assertIn("gh release create", workflow)
        self.assertNotIn("maturin sdist --locked", workflow)


if __name__ == "__main__":
    unittest.main()

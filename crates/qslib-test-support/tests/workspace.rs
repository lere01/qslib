use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::process::Command;

use serde_json::Value;

const EXPECTED_PACKAGES: [&str; 9] = [
    "qslib-quantum",
    "qslib-quantum-cli",
    "qslib-quantum-core",
    "qslib-quantum-exact",
    "qslib-quantum-io",
    "qslib-quantum-python",
    "qslib-quantum-sse",
    "qslib-quantum-variational",
    "qslib-test-support",
];

#[test]
fn workspace_metadata_matches_the_accepted_package_and_target_map() {
    let metadata = cargo_metadata();
    let packages = package_map(&metadata);
    let actual: BTreeSet<_> = packages.keys().copied().collect();
    let expected: BTreeSet<_> = EXPECTED_PACKAGES.into_iter().collect();
    assert_eq!(actual, expected);

    // Registry publication is allowed for the facade and capability crates.
    // The Python binding ships through PyPI and test support stays private.
    let registry_excluded: BTreeSet<&str> = ["qslib-quantum-python", "qslib-test-support"]
        .into_iter()
        .collect();
    for (name, package) in &packages {
        assert_eq!(package["version"], "0.2.0");
        assert_eq!(package["edition"], "2024");
        assert_eq!(package["rust_version"], "1.85");
        assert_eq!(package["license"], "Apache-2.0");
        assert_eq!(package["repository"], "https://github.com/lere01/qslib.git");
        if registry_excluded.contains(name) {
            assert_eq!(
                package["publish"],
                serde_json::json!([]),
                "{name} must stay excluded from registry publication"
            );
        } else {
            assert_eq!(
                package["publish"],
                Value::Null,
                "{name} must be publishable to crates.io"
            );
        }
        // Every crate page must carry its own role-scoped README. The Python
        // binding ships its PyPI readme through pyproject.toml instead.
        let manifest_dir = Path::new(package["manifest_path"].as_str().expect("manifest path"))
            .parent()
            .expect("manifest directory")
            .to_path_buf();
        if *name == "qslib-quantum-python" {
            assert!(
                manifest_dir.join("PYTHON_README.md").is_file(),
                "{name} must keep its PyPI readme"
            );
            // The PyPI distribution and the Python contract test carry their
            // own version strings; both must track the workspace version so a
            // release bump cannot fail late in the wheel smoke test.
            let version = package["version"].as_str().expect("package version");
            let pyproject = std::fs::read_to_string(manifest_dir.join("pyproject.toml"))
                .expect("read pyproject.toml");
            assert!(
                pyproject.contains(&format!("version = \"{version}\"")),
                "{name} pyproject.toml must declare version {version}"
            );
            let contract = std::fs::read_to_string(manifest_dir.join("tests/python_contract.py"))
                .expect("read python contract test");
            assert!(
                contract.contains(&format!("qslib.__version__ == \"{version}\"")),
                "{name} contract test must pin __version__ {version}"
            );
        } else {
            assert_eq!(
                package["readme"], "README.md",
                "{name} must declare a crate README"
            );
            assert!(
                manifest_dir.join("README.md").is_file(),
                "{name} must have a README.md next to its manifest"
            );
        }
    }

    assert_target(packages["qslib-quantum"], "qslib", "lib");
    let root_targets = packages["qslib-quantum"]["targets"]
        .as_array()
        .expect("root targets");
    assert_eq!(
        root_targets
            .iter()
            .filter(|target| target["kind"] == serde_json::json!(["lib"]))
            .count(),
        1
    );
    assert!(!root_targets.iter().any(|target| {
        target["kind"]
            .as_array()
            .is_some_and(|kinds| kinds.iter().any(|kind| kind == "bin"))
    }));
    assert_target(packages["qslib-quantum-cli"], "qslib", "bin");
    assert_target(packages["qslib-quantum-cli"], "qslib_cli", "lib");
    assert_target(packages["qslib-quantum-python"], "qslib_quantum", "rlib");
    assert_target(packages["qslib-test-support"], "qslib_test_support", "lib");

    let manifest =
        std::fs::read_to_string(workspace_root().join("Cargo.toml")).expect("read root manifest");
    assert!(manifest.contains("resolver = \"3\""));
}

#[test]
fn dependency_edges_follow_the_accepted_one_way_architecture() {
    let metadata = cargo_metadata();
    let packages = package_map(&metadata);
    assert_eq!(
        normal_local_dependencies(packages["qslib-quantum-core"]),
        set(&[])
    );
    for algorithm in [
        "qslib-quantum-exact",
        "qslib-quantum-io",
        "qslib-quantum-sse",
        "qslib-quantum-variational",
    ] {
        assert_eq!(
            normal_local_dependencies(packages[algorithm]),
            set(&["qslib-quantum-core"]),
            "unexpected edge from {algorithm}"
        );
    }
    assert_eq!(
        normal_local_dependencies(packages["qslib-quantum"]),
        set(&[
            "qslib-quantum-core",
            "qslib-quantum-exact",
            "qslib-quantum-io",
            "qslib-quantum-sse",
            "qslib-quantum-variational",
        ])
    );
    for interface in ["qslib-quantum-cli", "qslib-quantum-python"] {
        assert_eq!(
            normal_local_dependencies(packages[interface]),
            set(&["qslib-quantum"]),
            "interface must consume only the facade"
        );
    }
    for (name, package) in &packages {
        if *name != "qslib-test-support" {
            assert!(
                !normal_local_dependencies(package).contains("qslib-test-support"),
                "production package {name} depends on test support"
            );
        }
    }
}

#[test]
fn facade_features_are_additive_and_core_only_is_lightweight() {
    let metadata = cargo_metadata();
    let packages = package_map(&metadata);
    let features = packages["qslib-quantum"]["features"]
        .as_object()
        .expect("facade features");
    assert_eq!(features["default"], serde_json::json!([]));
    assert_eq!(features["exact"], serde_json::json!(["dep:qslib-exact"]));
    assert_eq!(
        features["variational"],
        serde_json::json!(["dep:qslib-variational"])
    );
    assert_eq!(features["sse"], serde_json::json!(["dep:qslib-sse"]));
    assert_eq!(features["io"], serde_json::json!(["dep:qslib-io"]));
    assert_eq!(
        features["full"],
        serde_json::json!(["exact", "variational", "sse", "io"])
    );

    let output = Command::new(env!("CARGO"))
        .args([
            "tree",
            "--locked",
            "-p",
            "qslib-quantum",
            "--no-default-features",
            "--prefix",
            "none",
        ])
        .current_dir(workspace_root())
        .output()
        .expect("run cargo tree");
    assert!(output.status.success(), "cargo tree failed");
    let tree = String::from_utf8(output.stdout).expect("cargo tree UTF-8");
    assert!(tree.contains("qslib-quantum v0.2.0"));
    assert!(tree.contains("qslib-quantum-core v0.2.0"));
    for forbidden in [
        "qslib-quantum-exact",
        "qslib-quantum-io",
        "qslib-quantum-sse",
        "qslib-quantum-variational",
        "rayon",
        "pyo3",
        "parquet",
        "arrow",
    ] {
        assert!(
            !tree.contains(forbidden),
            "core-only tree includes {forbidden}:\n{tree}"
        );
    }
}

fn cargo_metadata() -> Value {
    let output = Command::new(env!("CARGO"))
        .args(["metadata", "--locked", "--no-deps", "--format-version", "1"])
        .current_dir(workspace_root())
        .output()
        .expect("run cargo metadata");
    assert!(output.status.success(), "cargo metadata failed");
    serde_json::from_slice(&output.stdout).expect("parse cargo metadata")
}

fn package_map(metadata: &Value) -> BTreeMap<&str, &Value> {
    metadata["packages"]
        .as_array()
        .expect("metadata packages")
        .iter()
        .map(|package| (package["name"].as_str().expect("package name"), package))
        .collect()
}

fn normal_local_dependencies(package: &Value) -> BTreeSet<&str> {
    package["dependencies"]
        .as_array()
        .expect("package dependencies")
        .iter()
        .filter(|dependency| dependency["kind"].is_null() && dependency["path"].is_string())
        .map(|dependency| dependency["name"].as_str().expect("dependency name"))
        .collect()
}

fn assert_target(package: &Value, name: &str, kind: &str) {
    let found = package["targets"]
        .as_array()
        .expect("package targets")
        .iter()
        .any(|target| {
            target["name"] == name
                && target["kind"]
                    .as_array()
                    .is_some_and(|kinds| kinds.iter().any(|candidate| candidate == kind))
        });
    assert!(found, "missing {kind} target {name}");
}

fn set<'a>(items: &[&'a str]) -> BTreeSet<&'a str> {
    items.iter().copied().collect()
}

fn workspace_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("canonical workspace root")
}

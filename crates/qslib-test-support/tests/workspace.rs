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

    for package in packages.values() {
        assert_eq!(package["version"], "0.1.0");
        assert_eq!(package["edition"], "2024");
        assert_eq!(package["rust_version"], "1.85");
        assert_eq!(package["license"], "Apache-2.0");
        assert_eq!(package["repository"], "https://github.com/lere01/qslib.git");
        assert_eq!(package["publish"], serde_json::json!([]));
    }

    assert_target(packages["qslib-quantum"], "qslib", "lib");
    assert_eq!(
        packages["qslib-quantum"]["targets"]
            .as_array()
            .unwrap()
            .len(),
        1
    );
    assert_target(packages["qslib-quantum-cli"], "qslib", "bin");
    assert_target(packages["qslib-quantum-python"], "qslib_quantum", "lib");
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
    let package_lines: Vec<_> = tree
        .lines()
        .filter(|line| !line.trim().is_empty())
        .collect();
    assert_eq!(
        package_lines.len(),
        2,
        "core-only tree was not lightweight:\n{tree}"
    );
    assert!(tree.contains("qslib-quantum v0.1.0"));
    assert!(tree.contains("qslib-quantum-core v0.1.0"));
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

use qslib_cli::run;
use serde_json::Value;
use std::path::PathBuf;

fn example(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../examples")
        .join(name)
        .display()
        .to_string()
}

fn run_json(command: &[&str]) -> Value {
    let args = command
        .iter()
        .map(|value| (*value).to_owned())
        .collect::<Vec<_>>();
    serde_json::from_str(&run(&args).expect("example command should succeed")).expect("JSON output")
}

#[test]
fn documented_exact_examples_execute_with_stable_semantics() {
    let ground = run_json(&["exact", "ground-state", &example("tfim_4.yaml"), "--json"]);
    assert_eq!(ground["model"], "tfim");
    assert_eq!(ground["site_count"], 4);
    assert!(ground["residual"].as_f64().unwrap() < 1.0e-10);

    let heisenberg = run_json(&[
        "model",
        "validate",
        &example("heisenberg_disordered_4.yaml"),
        "--json",
    ]);
    assert_eq!(heisenberg["model"], "heisenberg");
    assert_eq!(heisenberg["interaction_terms"], 6);
    assert_eq!(heisenberg["provenance"]["site_order"], "row_major");
    assert!(
        heisenberg["resolved_interactions"]
            .as_array()
            .unwrap()
            .iter()
            .any(|term| term["coefficient"] == -0.5)
    );

    let evolved = run_json(&[
        "exact",
        "evolve",
        &example("tfim_4.yaml"),
        "--t-max",
        "0.1",
        "--json",
    ]);
    assert_eq!(evolved["imaginary_time"], false);
    assert!((evolved["norm"].as_f64().unwrap() - 1.0).abs() < 1.0e-10);
}

#[test]
fn documented_sse_smoke_example_executes() {
    let result = run_json(&["sse", "run", &example("tfim_thermal_4.yaml"), "--json"]);
    assert_eq!(result["model"], "tfim");
    assert_eq!(result["chains"], 2);
    assert!(result["energy_per_site"].as_f64().unwrap().is_finite());
}

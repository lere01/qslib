use std::collections::BTreeSet;

use qslib_test_support::{
    CONVENTION_SCHEMA, FIXTURE_SCHEMA, Fixture, FixtureKind, load_conformance_fixtures,
    required_fixture_kinds, validate_fixture_set,
};

#[test]
fn every_required_neutral_fixture_has_a_unique_valid_envelope() {
    let fixtures = load_conformance_fixtures().expect("the committed fixtures must load");
    let required: BTreeSet<_> = required_fixture_kinds().into_iter().collect();
    let actual: BTreeSet<_> = fixtures.iter().map(|fixture| fixture.kind).collect();

    assert_eq!(
        actual, required,
        "one fixture is required for every M1 kind"
    );
    assert_eq!(
        fixtures.len(),
        required.len(),
        "fixture kinds must be unique"
    );

    for fixture in fixtures {
        assert_eq!(fixture.fixture_schema, FIXTURE_SCHEMA);
        assert_eq!(fixture.convention_schema, CONVENTION_SCHEMA);
        assert!(!fixture.id.trim().is_empty());
        assert!(!fixture.claim.trim().is_empty());
        assert!(!fixture.oracle.method.trim().is_empty());
        assert!(fixture.oracle.source.starts_with("docs/conventions.md#24"));
        assert!(!fixture.oracle.source_revision.trim().is_empty());
        assert!(!fixture.oracle.derivation.trim().is_empty());
        assert!(fixture.oracle.independent_of_production);
        assert!(!fixture.authorship.prepared_by.trim().is_empty());
        assert!(!fixture.authorship.reviewed_by.trim().is_empty());
        assert!(fixture.conventions.is_object());
        assert!(fixture.data.is_object());
        fixture
            .validate()
            .expect("a committed fixture must validate");
    }
}

#[test]
fn malformed_or_ambiguous_fixture_envelopes_are_rejected() {
    let bit_fixture = include_str!("../../../fixtures/conformance/v1/bit_packing.json");
    let rotation_fixture =
        include_str!("../../../fixtures/conformance/v1/basis_rotation_spectrum.json");
    let tfim_fixture = include_str!("../../../fixtures/conformance/v1/tfim_one_bond.json");

    assert_error_contains(
        &bit_fixture.replace(FIXTURE_SCHEMA, "future-fixture-v2"),
        "unsupported fixture schema",
    );
    assert_error_contains(
        &bit_fixture.replace(CONVENTION_SCHEMA, "future-conventions-v2"),
        "unsupported convention schema",
    );
    assert_error_contains(
        &bit_fixture.replace(
            "\"independent_of_production\": true",
            "\"independent_of_production\": false",
        ),
        "independent of production",
    );
    assert_error_contains(
        &rotation_fixture.replace(
            "\"absolute_tolerance\": 1e-12",
            "\"absolute_tolerance\": 0.0",
        ),
        "absolute tolerance",
    );
    assert_error_contains(
        &tfim_fixture.replace("\"shape\": [4, 4]", "\"shape\": [5, 4]"),
        "matrix shape",
    );
    assert_error_contains(
        &bit_fixture.replace("\"data\": {", "\"unexpected\": true, \"data\": {"),
        "missing or unknown top-level fields",
    );
    assert_error_contains(
        &bit_fixture.replace("\"packed_mask\": 13", "\"packed_mask\": NaN"),
        "invalid fixture JSON",
    );

    let mut missing_provenance: serde_json::Value =
        serde_json::from_str(bit_fixture).expect("fixture JSON");
    missing_provenance["oracle"]
        .as_object_mut()
        .expect("oracle object")
        .remove("source_revision");
    assert_error_contains(
        &serde_json::to_string(&missing_provenance).expect("serialize mutation"),
        "source_revision",
    );
}

#[test]
fn manifest_rejects_wrong_checksums_and_duplicate_fixture_ids() {
    let manifest = include_str!("../../../fixtures/conformance/v1/manifest.json");
    let files = fixture_files();
    let wrong_checksum = manifest.replacen(
        "36d146c7c4b58ec61bb07d64e54726fe6721e0f14803eb483953e3ee0eb335dc",
        "0000000000000000000000000000000000000000000000000000000000000000",
        1,
    );
    let error = validate_fixture_set(&wrong_checksum, &files).expect_err("checksum must fail");
    assert!(error.to_string().contains("BLAKE3 mismatch"));

    let duplicate_id = "tfim_two_site_hadamard_spectrum";
    let mut owned_files: Vec<(String, Vec<u8>)> = files
        .iter()
        .map(|(path, bytes)| ((*path).to_owned(), bytes.to_vec()))
        .collect();
    let bit_bytes = &mut owned_files
        .iter_mut()
        .find(|(path, _)| path == "bit_packing.json")
        .expect("bit fixture")
        .1;
    let bit_text = std::str::from_utf8(bit_bytes).expect("UTF-8 fixture");
    *bit_bytes = bit_text
        .replace("bit_packing_1011", duplicate_id)
        .into_bytes();
    let digest = blake3::hash(bit_bytes).to_hex().to_string();
    let mut manifest_value: serde_json::Value = serde_json::from_str(manifest).expect("manifest");
    let bit_entry = manifest_value["entries"]
        .as_array_mut()
        .expect("entries")
        .iter_mut()
        .find(|entry| entry["path"] == "bit_packing.json")
        .expect("bit manifest entry");
    bit_entry["id"] = duplicate_id.into();
    bit_entry["blake3"] = digest.into();
    let owned_views: Vec<(&str, &[u8])> = owned_files
        .iter()
        .map(|(path, bytes)| (path.as_str(), bytes.as_slice()))
        .collect();
    let error = validate_fixture_set(
        &serde_json::to_string(&manifest_value).expect("serialize manifest"),
        &owned_views,
    )
    .expect_err("duplicate fixture id must fail");
    assert!(error.to_string().contains("duplicate fixture id"));
}

#[test]
fn documented_json_schema_is_well_formed() {
    let schema = include_str!("../../../fixtures/conformance/v1/_schema.json");
    let parsed: serde_json::Value = serde_json::from_str(schema).expect("valid JSON schema");
    assert_eq!(
        parsed["$schema"],
        "https://json-schema.org/draft/2020-12/schema"
    );
    assert_eq!(parsed["additionalProperties"], false);
}

#[test]
fn fixture_payloads_preserve_the_release_critical_reference_values() {
    let fixtures = load_conformance_fixtures().expect("the committed fixtures must load");
    let by_kind = |kind: FixtureKind| {
        fixtures
            .iter()
            .find(|fixture| fixture.kind == kind)
            .expect("required fixture kind")
    };

    assert_eq!(by_kind(FixtureKind::BitPacking).data["packed_mask"], 13);
    assert_eq!(
        by_kind(FixtureKind::HeisenbergHeterogeneous).data["diagonal_energy"],
        0.25
    );
    assert_eq!(
        by_kind(FixtureKind::RydbergTwoSite).data["matrix"]["entries"][15]["re"],
        -1.0
    );
    assert_eq!(
        by_kind(FixtureKind::ObservableNormalization).data["connected_zz_for_every_pair"],
        0.0
    );
}

fn assert_error_contains(json: &str, expected: &str) {
    let error = Fixture::from_json(json).expect_err("mutated fixture must fail");
    assert!(
        error.to_string().contains(expected),
        "expected error containing {expected:?}, got {error}"
    );
}

fn fixture_files() -> Vec<(&'static str, &'static [u8])> {
    vec![
        (
            "basis_rotation_spectrum.json",
            include_bytes!("../../../fixtures/conformance/v1/basis_rotation_spectrum.json"),
        ),
        (
            "bit_packing.json",
            include_bytes!("../../../fixtures/conformance/v1/bit_packing.json"),
        ),
        (
            "heisenberg_heterogeneous.json",
            include_bytes!("../../../fixtures/conformance/v1/heisenberg_heterogeneous.json"),
        ),
        (
            "heisenberg_one_bond.json",
            include_bytes!("../../../fixtures/conformance/v1/heisenberg_one_bond.json"),
        ),
        (
            "observable_normalization.json",
            include_bytes!("../../../fixtures/conformance/v1/observable_normalization.json"),
        ),
        (
            "rectangular_indexing.json",
            include_bytes!("../../../fixtures/conformance/v1/rectangular_indexing.json"),
        ),
        (
            "rydberg_two_site.json",
            include_bytes!("../../../fixtures/conformance/v1/rydberg_two_site.json"),
        ),
        (
            "tfim_one_bond.json",
            include_bytes!("../../../fixtures/conformance/v1/tfim_one_bond.json"),
        ),
    ]
}

fn main() {
    let fixture_root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join("fixtures/conformance/v1");
    for name in [
        "basis_rotation_spectrum.json",
        "bit_packing.json",
        "heisenberg_heterogeneous.json",
        "heisenberg_one_bond.json",
        "observable_normalization.json",
        "rectangular_indexing.json",
        "rydberg_two_site.json",
        "tfim_one_bond.json",
    ] {
        let bytes = std::fs::read(fixture_root.join(name)).expect("read fixture");
        println!("{}  {name}", blake3::hash(&bytes).to_hex());
    }
}

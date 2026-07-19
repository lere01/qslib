use qslib::{BasisBit, BasisState, SiteCount, SiteId};

#[test]
fn stable_core_types_are_available_from_the_facade() {
    let sites = SiteCount::new(2).expect("sites");
    sites.validate(SiteId::new(1)).expect("site");
    let state = BasisState::from_bits(&[BasisBit::Zero, BasisBit::One]).expect("state");
    assert_eq!(state.hamming_weight(), 1);
}

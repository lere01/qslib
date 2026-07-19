use qslib_core::Complex64;
use qslib_variational::{
    ComplexWeightedMoments, RealizationEstimate, StatisticsError, WeightedMoments, autocorrelation,
    disorder_average, disorder_average_with_uncertainty, r_hat,
};

#[test]
fn weighted_online_moments_match_batch_and_merge_without_losing_weights() {
    let values = [(1.0, 1.0), (3.0, 2.0), (5.0, 1.0)];
    let mut online = WeightedMoments::new().unwrap();
    for (value, weight) in values {
        online.update(value, weight).unwrap();
    }
    assert!((online.mean().unwrap() - 3.0).abs() < 1.0e-12);
    assert!((online.variance().unwrap() - 2.0).abs() < 1.0e-12);
    let mut left = WeightedMoments::new().unwrap();
    left.update(1.0, 1.0).unwrap();
    let mut right = WeightedMoments::new().unwrap();
    right.update(3.0, 2.0).unwrap();
    right.update(5.0, 1.0).unwrap();
    left.merge(&right).unwrap();
    assert!((left.mean().unwrap() - online.mean().unwrap()).abs() < 1.0e-12);
    assert!((left.variance().unwrap() - online.variance().unwrap()).abs() < 1.0e-12);
    assert!(WeightedMoments::new().unwrap().mean().is_err());
    assert!(WeightedMoments::new().unwrap().update(1.0, -1.0).is_err());
}

#[test]
fn moment_overflow_is_error_atomic_for_update_and_merge() {
    let mut accumulator = WeightedMoments::new().unwrap();
    accumulator.update(1.0, f64::MAX).unwrap();
    assert!(matches!(
        accumulator.update(2.0, f64::MAX),
        Err(StatisticsError::NonFinite("accumulated weight"))
    ));
    assert_eq!(accumulator.mean().unwrap(), 1.0);
    assert_eq!(accumulator.weight(), f64::MAX);

    let mut left = WeightedMoments::new().unwrap();
    left.update(1.0, f64::MAX * 0.75).unwrap();
    let mut right = WeightedMoments::new().unwrap();
    right.update(2.0, f64::MAX * 0.75).unwrap();
    assert!(left.merge(&right).is_err());
    assert_eq!(left.mean().unwrap(), 1.0);
}

#[test]
fn complex_moments_preserve_real_and_imaginary_statistics() {
    let mut moments = ComplexWeightedMoments::new().unwrap();
    moments.update(Complex64::new(1.0, 2.0), 1.0).unwrap();
    moments.update(Complex64::new(3.0, 4.0), 1.0).unwrap();
    assert_eq!(moments.mean().unwrap(), Complex64::new(2.0, 3.0));
    assert!(moments.imaginary_variance().unwrap() > 0.0);
}

#[test]
fn complex_moment_overflow_is_component_atomic() {
    let mut moments = ComplexWeightedMoments::new().unwrap();
    moments.update(Complex64::new(0.0, 0.0), 1.0).unwrap();
    assert!(moments.update(Complex64::new(1.0, f64::MAX), 1.0).is_err());
    assert_eq!(moments.mean().unwrap(), Complex64::new(0.0, 0.0));
}

#[test]
fn autocorrelation_reports_correlation_time_and_reduced_ess() {
    let samples = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
    let estimate = autocorrelation(&samples, 4).unwrap();
    assert!((estimate.integrated_time() - 11.0 / 6.0).abs() < 1.0e-12);
    assert!(estimate.effective_sample_size() < samples.len() as f64);
    assert_eq!(
        estimate.algorithm(),
        "geyer_initial_positive_sequence_common_n_tau_floor_1"
    );
}

#[test]
fn autocorrelation_keeps_positive_time_for_alternating_and_short_boundaries() {
    let alternating = autocorrelation(&[1.0, -1.0, 1.0, -1.0], 3).unwrap();
    assert_eq!(alternating.integrated_time(), 1.0);
    assert_eq!(alternating.effective_sample_size(), 4.0);
    let short = autocorrelation(&[0.0, 1.0], 1).unwrap();
    assert!(short.integrated_time().is_finite());
    assert!(short.effective_sample_size() > 0.0);
}

#[test]
fn r_hat_detects_between_chain_failure_and_accepts_identical_chains() {
    let good = vec![vec![1.0, 2.0, 3.0, 4.0], vec![1.0, 2.0, 3.0, 4.0]];
    let good_diagnostic = r_hat(&good).unwrap();
    assert!((good_diagnostic.value() - 1.0).abs() < 1.0e-12);
    assert_eq!(
        good_diagnostic.algorithm(),
        "gelman_rubin_classic_equal_length"
    );
    let bad = vec![vec![1.0, 1.0, 1.0, 1.0], vec![5.0, 5.0, 5.0, 5.0]];
    assert!(r_hat(&bad).unwrap().value() > 1.1);
}

#[test]
fn disorder_average_separates_realization_mean_and_between_realization_error() {
    let summary = disorder_average(&[("a", 1.0, 1.0), ("b", 3.0, 1.0)]).unwrap();
    assert_eq!(summary.mean(), 2.0);
    assert!(summary.between_realization_variance().unwrap() > 0.0);
    assert_eq!(summary.realizations(), &["a".to_string(), "b".to_string()]);
    assert_eq!(summary.sampling_variance(), None);
    let one =
        disorder_average_with_uncertainty(&[
            RealizationEstimate::new("a", 1.0, 1.0, Some(0.25)).unwrap()
        ])
        .unwrap();
    assert_eq!(one.between_realization_variance(), None);
    assert_eq!(one.sampling_variance(), Some(0.25));
    let weighted = disorder_average_with_uncertainty(&[
        RealizationEstimate::new("a", 1.0, 1.0, Some(1.0)).unwrap(),
        RealizationEstimate::new("b", 3.0, 1.0, Some(1.0)).unwrap(),
    ])
    .unwrap();
    assert_eq!(weighted.sampling_variance(), Some(0.5));
    assert!(matches!(
        disorder_average(&[("a", 1.0, 1.0), ("a", 2.0, 1.0)]),
        Err(StatisticsError::DuplicateRealization)
    ));
}

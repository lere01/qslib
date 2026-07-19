use qslib_core::Complex64;
use qslib_variational::{
    ComplexWeightedMoments, WeightedMoments, autocorrelation, disorder_average, r_hat,
};

#[test]
fn weighted_online_moments_match_batch_and_merge_without_losing_weights() {
    let values = [(1.0, 1.0), (3.0, 2.0), (5.0, 1.0)];
    let mut online = WeightedMoments::new().unwrap();
    for (value, weight) in values {
        online.update(value, weight).unwrap();
    }
    assert!((online.mean() - 3.0).abs() < 1.0e-12);
    assert!((online.variance() - 2.0).abs() < 1.0e-12);
    let mut left = WeightedMoments::new().unwrap();
    left.update(1.0, 1.0).unwrap();
    let mut right = WeightedMoments::new().unwrap();
    right.update(3.0, 2.0).unwrap();
    right.update(5.0, 1.0).unwrap();
    left.merge(&right).unwrap();
    assert!((left.mean() - online.mean()).abs() < 1.0e-12);
    assert!((left.variance() - online.variance()).abs() < 1.0e-12);
    assert!(WeightedMoments::new().unwrap().update(1.0, -1.0).is_err());
}

#[test]
fn complex_moments_preserve_real_and_imaginary_statistics() {
    let mut moments = ComplexWeightedMoments::new().unwrap();
    moments.update(Complex64::new(1.0, 2.0), 1.0).unwrap();
    moments.update(Complex64::new(3.0, 4.0), 1.0).unwrap();
    assert_eq!(moments.mean(), Complex64::new(2.0, 3.0));
    assert!(moments.imaginary_variance() > 0.0);
}

#[test]
fn autocorrelation_reports_correlation_time_and_reduced_ess() {
    let samples = vec![0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0];
    let estimate = autocorrelation(&samples, 4).unwrap();
    assert!(estimate.integrated_time() > 1.0);
    assert!(estimate.effective_sample_size() < samples.len() as f64);
}

#[test]
fn r_hat_detects_between_chain_failure_and_accepts_identical_chains() {
    let good = vec![vec![1.0, 2.0, 3.0, 4.0], vec![1.0, 2.0, 3.0, 4.0]];
    assert!((r_hat(&good).unwrap() - 1.0).abs() < 1.0e-12);
    let bad = vec![vec![1.0, 1.0, 1.0, 1.0], vec![5.0, 5.0, 5.0, 5.0]];
    assert!(r_hat(&bad).unwrap() > 1.1);
}

#[test]
fn disorder_average_separates_realization_mean_and_between_realization_error() {
    let summary = disorder_average(&[("a", 1.0, 1.0), ("b", 3.0, 1.0)]).unwrap();
    assert_eq!(summary.mean(), 2.0);
    assert!(summary.between_realization_variance() > 0.0);
    assert_eq!(summary.realizations(), &["a".to_string(), "b".to_string()]);
}

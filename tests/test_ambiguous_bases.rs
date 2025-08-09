use meta_forge::{AmbiguousBaseConfig, AmbiguousBaseStrategy, CanonicalKmer, MinimizerExtractor};

#[test]
fn test_skip_strategy() {
    let config = AmbiguousBaseConfig {
        strategy: AmbiguousBaseStrategy::Skip,
        max_n_count: 0,
        replacement_base: 'A',
        random_probabilities: None,
    };

    // Should work with no N's
    let result = CanonicalKmer::new_with_config("ATCG", &config);
    assert!(result.is_ok());

    // Should fail with N's
    let result = CanonicalKmer::new_with_config("ATCN", &config);
    assert!(result.is_err());
    assert!(result.unwrap_err().to_string().contains("ambiguous bases"));
}

#[test]
fn test_allow_strategy() {
    let config = AmbiguousBaseConfig {
        strategy: AmbiguousBaseStrategy::Allow,
        max_n_count: 2,
        replacement_base: 'A',
        random_probabilities: None,
    };

    // Should work with no N's
    let result = CanonicalKmer::new_with_config("ATCG", &config);
    assert!(result.is_ok());

    // Should work with few N's
    let result = CanonicalKmer::new_with_config("ATCN", &config);
    assert!(result.is_ok());

    // Should work with exactly max N's
    let result = CanonicalKmer::new_with_config("ATNN", &config);
    assert!(result.is_ok());

    // Should fail with too many N's
    let result = CanonicalKmer::new_with_config("ANNN", &config);
    assert!(result.is_err());
    assert!(result
        .unwrap_err()
        .to_string()
        .contains("too many ambiguous bases"));
}

#[test]
fn test_replace_strategy() {
    let config = AmbiguousBaseConfig {
        strategy: AmbiguousBaseStrategy::Replace,
        max_n_count: 0, // Not used in Replace strategy
        replacement_base: 'A',
        random_probabilities: None,
    };

    // Should replace N with A
    let result = CanonicalKmer::new_with_config("ATCN", &config);
    assert!(result.is_ok());
    let kmer = result.unwrap();
    assert!(kmer.sequence.contains('A'));
    assert!(!kmer.sequence.contains('N'));

    // Should replace multiple N's with A
    let result = CanonicalKmer::new_with_config("NNCN", &config);
    assert!(result.is_ok());
    let kmer = result.unwrap();
    assert!(!kmer.sequence.contains('N'));
    assert_eq!(kmer.sequence.chars().filter(|&c| c == 'A').count(), 3); // 2 N's replaced + original A = 3 A's (potentially canonical)
}

#[test]
fn test_context_replace_strategy() {
    let config = AmbiguousBaseConfig {
        strategy: AmbiguousBaseStrategy::ContextReplace,
        max_n_count: 0,
        replacement_base: 'A', // Not used in ContextReplace
        random_probabilities: None,
    };

    // Should replace N with most common base (A appears 3 times)
    let result = CanonicalKmer::new_with_config("AAAN", &config);
    assert!(result.is_ok());
    let kmer = result.unwrap();
    assert!(!kmer.sequence.contains('N'));

    // Should replace N with most common base (G appears twice)
    let result = CanonicalKmer::new_with_config("ATGN", &config);
    assert!(result.is_ok());
    let kmer = result.unwrap();
    assert!(!kmer.sequence.contains('N'));
}

#[test]
fn test_default_behavior() {
    // Default should be Allow with max_n_count = 2
    let result = CanonicalKmer::new("ATCN");
    assert!(result.is_ok());

    let result = CanonicalKmer::new("ATNN");
    assert!(result.is_ok());

    let result = CanonicalKmer::new("ANNN");
    assert!(result.is_err());
}

#[test]
fn test_minimizer_extraction_with_ambiguous_bases() {
    let config = AmbiguousBaseConfig {
        strategy: AmbiguousBaseStrategy::Allow,
        max_n_count: 1,
        replacement_base: 'A',
        random_probabilities: None,
    };

    let extractor = MinimizerExtractor::new_with_config(4, 6, config);
    let sequence = "ATCGNATCG"; // Contains one N

    // Should extract minimizers, skipping k-mers with too many N's
    let result = extractor.extract_minimizers(sequence);
    assert!(result.is_ok());

    let minimizers = result.unwrap();
    assert!(!minimizers.is_empty());
}

#[test]
fn test_original_error_sequence() {
    // Test the exact sequence that was causing the original error
    let problem_sequence = "AATACAAGCATCAAATNNNNNNGCTGTCTG";

    // With Skip strategy (original behavior) - should fail
    let config_skip = AmbiguousBaseConfig {
        strategy: AmbiguousBaseStrategy::Skip,
        max_n_count: 0,
        replacement_base: 'A',
        random_probabilities: None,
    };

    let extractor = MinimizerExtractor::new_with_config(21, 25, config_skip);
    let result = extractor.extract_minimizers(problem_sequence);
    // Should succeed but with fewer minimizers (skipping problematic k-mers)
    assert!(result.is_ok());

    // With Allow strategy - should work better
    let config_allow = AmbiguousBaseConfig {
        strategy: AmbiguousBaseStrategy::Allow,
        max_n_count: 3,
        replacement_base: 'A',
        random_probabilities: None,
    };

    let extractor = MinimizerExtractor::new_with_config(21, 25, config_allow);
    let result = extractor.extract_minimizers(problem_sequence);
    assert!(result.is_ok());

    // With Replace strategy - should work and replace N's
    let config_replace = AmbiguousBaseConfig {
        strategy: AmbiguousBaseStrategy::Replace,
        max_n_count: 0,
        replacement_base: 'A',
        random_probabilities: None,
    };

    let extractor = MinimizerExtractor::new_with_config(21, 25, config_replace);
    let result = extractor.extract_minimizers(problem_sequence);
    assert!(result.is_ok());
    let minimizers = result.unwrap();
    assert!(!minimizers.is_empty());
}

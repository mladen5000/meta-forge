use meta_forge::{AmbiguousBaseConfig, AmbiguousBaseStrategy, CanonicalKmer, MinimizerExtractor};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🧬 Testing K-mer Ambiguous Base Handling");

    // Test the original problematic sequence
    let problem_sequence = "AATACAAGCATCAAATNNNNNNGCTGTCTG";
    println!("Testing sequence: {}", problem_sequence);

    // Test 1: Skip strategy (original behavior)
    println!("\n1. Testing Skip Strategy:");
    let config_skip = AmbiguousBaseConfig {
        strategy: AmbiguousBaseStrategy::Skip,
        max_n_count: 0,
        replacement_base: 'A',
        random_probabilities: None,
    };

    let extractor = MinimizerExtractor::new_with_config(21, 25, config_skip);
    match extractor.extract_minimizers(problem_sequence) {
        Ok(minimizers) => {
            println!(
                "✅ Skip strategy succeeded! Found {} minimizers (skipped problematic k-mers)",
                minimizers.len()
            );
        }
        Err(e) => {
            println!("❌ Skip strategy failed: {}", e);
        }
    }

    // Test 2: Allow strategy
    println!("\n2. Testing Allow Strategy (max 3 N's):");
    let config_allow = AmbiguousBaseConfig {
        strategy: AmbiguousBaseStrategy::Allow,
        max_n_count: 3,
        replacement_base: 'A',
        random_probabilities: None,
    };

    let extractor = MinimizerExtractor::new_with_config(21, 25, config_allow);
    match extractor.extract_minimizers(problem_sequence) {
        Ok(minimizers) => {
            println!(
                "✅ Allow strategy succeeded! Found {} minimizers",
                minimizers.len()
            );
            for (i, min) in minimizers.iter().take(3).enumerate() {
                println!(
                    "   Minimizer {}: {} (pos: {})",
                    i + 1,
                    min.kmer.sequence,
                    min.position
                );
            }
        }
        Err(e) => {
            println!("❌ Allow strategy failed: {}", e);
        }
    }

    // Test 3: Replace strategy
    println!("\n3. Testing Replace Strategy (replace with A):");
    let config_replace = AmbiguousBaseConfig {
        strategy: AmbiguousBaseStrategy::Replace,
        max_n_count: 0,
        replacement_base: 'A',
        random_probabilities: None,
    };

    let extractor = MinimizerExtractor::new_with_config(21, 25, config_replace);
    match extractor.extract_minimizers(problem_sequence) {
        Ok(minimizers) => {
            println!(
                "✅ Replace strategy succeeded! Found {} minimizers",
                minimizers.len()
            );
            for (i, min) in minimizers.iter().take(3).enumerate() {
                println!(
                    "   Minimizer {}: {} (pos: {})",
                    i + 1,
                    min.kmer.sequence,
                    min.position
                );
            }
        }
        Err(e) => {
            println!("❌ Replace strategy failed: {}", e);
        }
    }

    // Test 4: Context replace strategy
    println!("\n4. Testing Context Replace Strategy:");
    let config_context = AmbiguousBaseConfig {
        strategy: AmbiguousBaseStrategy::ContextReplace,
        max_n_count: 0,
        replacement_base: 'A',
        random_probabilities: None,
    };

    let extractor = MinimizerExtractor::new_with_config(21, 25, config_context);
    match extractor.extract_minimizers(problem_sequence) {
        Ok(minimizers) => {
            println!(
                "✅ Context replace strategy succeeded! Found {} minimizers",
                minimizers.len()
            );
            for (i, min) in minimizers.iter().take(3).enumerate() {
                println!(
                    "   Minimizer {}: {} (pos: {})",
                    i + 1,
                    min.kmer.sequence,
                    min.position
                );
            }
        }
        Err(e) => {
            println!("❌ Context replace strategy failed: {}", e);
        }
    }

    // Test 5: Default behavior (should use Allow with max_n_count = 2)
    println!("\n5. Testing Default Behavior:");
    let extractor_default = MinimizerExtractor::new(21, 25);
    match extractor_default.extract_minimizers(problem_sequence) {
        Ok(minimizers) => {
            println!(
                "✅ Default behavior succeeded! Found {} minimizers",
                minimizers.len()
            );
        }
        Err(e) => {
            println!("❌ Default behavior failed: {}", e);
        }
    }

    // Test individual k-mers
    println!("\n6. Testing Individual K-mers:");
    println!("Testing k-mer with few N's: ATCGNATCG");
    match CanonicalKmer::new("ATCGNATCG") {
        Ok(kmer) => println!("✅ K-mer with 1 N succeeded: {}", kmer.sequence),
        Err(e) => println!("❌ K-mer with 1 N failed: {}", e),
    }

    println!("Testing k-mer with many N's: NNNNNNATCG");
    match CanonicalKmer::new("NNNNNNATCG") {
        Ok(kmer) => println!("✅ K-mer with 6 N's succeeded: {}", kmer.sequence),
        Err(e) => println!("❌ K-mer with 6 N's failed: {}", e),
    }

    println!("\n🎉 K-mer ambiguous base handling test completed!");

    Ok(())
}

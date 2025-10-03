//! QC Configuration Presets
//!
//! Provides idiomatic preset configurations for common use cases

use super::{AdapterConfig, QCPipelineConfig, QualityFilterConfig};

/// QC preset levels
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QCPreset {
    /// Strict quality control (Q30, 75bp min, 5% adapter error)
    Strict,
    /// Standard quality control (Q20, 50bp min, 10% adapter error) - default
    Standard,
    /// Lenient quality control (Q15, 30bp min, 15% adapter error)
    Lenient,
    /// No quality control (bypass all filtering)
    None,
}

impl From<QCPreset> for QualityFilterConfig {
    fn from(preset: QCPreset) -> Self {
        match preset {
            QCPreset::Strict => Self {
                min_quality: 30,
                window_size: 4,
                min_window_quality: 30.0,
                min_length: 75,
                min_avg_quality: 30.0,
                quality_offset: 33,
            },
            QCPreset::Standard => Self::default(),
            QCPreset::Lenient => Self {
                min_quality: 15,
                window_size: 4,
                min_window_quality: 15.0,
                min_length: 30,
                min_avg_quality: 20.0,
                quality_offset: 33,
            },
            QCPreset::None => Self {
                min_quality: 0,
                window_size: 1,
                min_window_quality: 0.0,
                min_length: 1,
                min_avg_quality: 0.0,
                quality_offset: 33,
            },
        }
    }
}

impl From<QCPreset> for AdapterConfig {
    fn from(preset: QCPreset) -> Self {
        use super::adapter_trimmer::*;

        match preset {
            QCPreset::Strict => Self {
                adapters: vec![
                    ILLUMINA_TRUSEQ_R1.to_string(),
                    ILLUMINA_TRUSEQ_R2.to_string(),
                    ILLUMINA_SMALL_RNA.to_string(),
                    ILLUMINA_SMALL_DNA_1.to_string(),
                    ILLUMINA_SMALL_DNA_2.to_string(),
                    ILLUMINA_SMALL_DNA_3.to_string(),
                    NEXTERA_R1.to_string(),
                ],
                min_overlap: 10,
                max_error_rate: 0.05,
                min_adapter_length: 8,
            },
            QCPreset::Standard => Self::default(),
            QCPreset::Lenient => Self {
                adapters: vec![
                    ILLUMINA_TRUSEQ_R1.to_string(),
                    ILLUMINA_TRUSEQ_R2.to_string(),
                    ILLUMINA_SMALL_DNA_1.to_string(),
                    ILLUMINA_SMALL_DNA_2.to_string(),
                    ILLUMINA_SMALL_DNA_3.to_string(),
                ],
                min_overlap: 6,
                max_error_rate: 0.15,
                min_adapter_length: 5,
            },
            QCPreset::None => Self {
                adapters: vec![],
                min_overlap: 100,
                max_error_rate: 0.0,
                min_adapter_length: 100,
            },
        }
    }
}

impl From<QCPreset> for QCPipelineConfig {
    fn from(preset: QCPreset) -> Self {
        let enabled = preset != QCPreset::None;

        Self {
            enable_quality_filter: enabled,
            enable_adapter_trimming: enabled,
            quality_config: preset.into(),
            adapter_config: preset.into(),
            verbose: false,
        }
    }
}

// Implement Display for better error messages and logging
impl std::fmt::Display for QCPreset {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            QCPreset::Strict => write!(f, "Strict (Q30, 75bp min)"),
            QCPreset::Standard => write!(f, "Standard (Q20, 50bp min)"),
            QCPreset::Lenient => write!(f, "Lenient (Q15, 30bp min)"),
            QCPreset::None => write!(f, "Disabled"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_quality_preset_conversion() {
        let config: QualityFilterConfig = QCPreset::Strict.into();
        assert_eq!(config.min_quality, 30);
        assert_eq!(config.min_length, 75);

        let config: QualityFilterConfig = QCPreset::Lenient.into();
        assert_eq!(config.min_quality, 15);
        assert_eq!(config.min_length, 30);
    }

    #[test]
    fn test_adapter_preset_conversion() {
        let config: AdapterConfig = QCPreset::Strict.into();
        assert_eq!(config.min_overlap, 10);
        assert_eq!(config.max_error_rate, 0.05);

        let config: AdapterConfig = QCPreset::None.into();
        assert!(config.adapters.is_empty());
    }

    #[test]
    fn test_pipeline_preset_conversion() {
        let config: QCPipelineConfig = QCPreset::Standard.into();
        assert!(config.enable_quality_filter);
        assert!(config.enable_adapter_trimming);

        let config: QCPipelineConfig = QCPreset::None.into();
        assert!(!config.enable_quality_filter);
        assert!(!config.enable_adapter_trimming);
    }

    #[test]
    fn test_display() {
        assert_eq!(QCPreset::Strict.to_string(), "Strict (Q30, 75bp min)");
        assert_eq!(QCPreset::None.to_string(), "Disabled");
    }
}

# Machine Learning Module Guidelines

## Purpose
Neural network models for repeat resolution and learned data structures.

## Key Components
- `gnn_repeat_resolution.rs` - Graph neural networks for assembly
- `learned_bloom_filter.rs` - ML-enhanced probabilistic data structures

## Development Rules
- Use stable numerical algorithms for training
- Validate model inputs and outputs carefully
- Handle model loading/saving errors gracefully
- Document model architectures and training procedures

## Model Design Principles
- Keep models lightweight and fast for inference
- Design for incremental learning when possible
- Use appropriate activation functions and regularization
- Consider model interpretability for biological insights

## Performance Requirements
- Optimize inference for real-time analysis
- Use SIMD/vectorization for model operations
- Consider GPU acceleration for training
- Profile memory usage during inference

## Training Considerations
- Use robust loss functions for genomic data
- Implement early stopping and regularization
- Validate on held-out biological datasets
- Document training hyperparameters and procedures

## Integration Guidelines
- Provide fallback methods when ML models fail
- Make ML features optional with feature flags
- Handle model version compatibility
- Test with both synthetic and real genomic data
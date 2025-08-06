# REST API Directory Structure and Organization

## Complete Directory Structure

```
src/
├── api/                          # API layer - HTTP handling
│   ├── mod.rs                    # API module exports
│   ├── server.rs                 # Server setup and routing
│   ├── routes/                   # Route definitions
│   │   ├── mod.rs
│   │   ├── analysis.rs           # Analysis pipeline routes
│   │   ├── samples.rs            # Sample management routes
│   │   ├── assemblies.rs         # Assembly operation routes
│   │   ├── annotations.rs        # Taxonomic classification routes
│   │   ├── jobs.rs               # Background job routes
│   │   ├── users.rs              # User management routes
│   │   ├── admin.rs              # Admin operation routes
│   │   └── health.rs             # Health check routes
│   ├── handlers/                 # Request handlers
│   │   ├── mod.rs
│   │   ├── analysis_handlers.rs
│   │   ├── sample_handlers.rs
│   │   ├── assembly_handlers.rs
│   │   ├── annotation_handlers.rs
│   │   ├── job_handlers.rs
│   │   ├── user_handlers.rs
│   │   ├── admin_handlers.rs
│   │   └── health_handlers.rs
│   ├── middleware/               # HTTP middleware
│   │   ├── mod.rs
│   │   ├── auth.rs               # JWT authentication
│   │   ├── cors.rs               # CORS handling
│   │   ├── logging.rs            # Request/response logging
│   │   ├── rate_limiting.rs      # Rate limiting
│   │   ├── validation.rs         # Input validation
│   │   ├── error_handling.rs     # Error handling middleware
│   │   └── metrics.rs            # Metrics collection
│   ├── extractors/               # Custom axum extractors
│   │   ├── mod.rs
│   │   ├── auth.rs               # Authentication extractor
│   │   ├── pagination.rs         # Pagination extractor
│   │   └── validation.rs         # Validation extractor
│   └── responses/                # Response utilities
│       ├── mod.rs
│       ├── api_response.rs       # Standard response wrapper
│       ├── error_response.rs     # Error response formatting
│       └── pagination.rs         # Pagination response helper
│
├── services/                     # Business logic layer
│   ├── mod.rs
│   ├── analysis_service.rs       # Analysis orchestration
│   ├── sample_service.rs         # Sample lifecycle management
│   ├── assembly_service.rs       # Assembly operations
│   ├── annotation_service.rs     # Taxonomic classification
│   ├── job_service.rs            # Background job management
│   ├── user_service.rs           # User management
│   ├── file_service.rs           # File handling
│   ├── notification_service.rs   # Notifications (email, webhooks)
│   └── traits/                   # Service traits for testing
│       ├── mod.rs
│       ├── analysis_service_trait.rs
│       ├── sample_service_trait.rs
│       └── user_service_trait.rs
│
├── data/                         # Data access layer
│   ├── mod.rs
│   ├── database.rs               # Database connection management
│   ├── repositories/             # Repository implementations
│   │   ├── mod.rs
│   │   ├── analysis_repository.rs
│   │   ├── sample_repository.rs
│   │   ├── assembly_repository.rs
│   │   ├── annotation_repository.rs
│   │   ├── job_repository.rs
│   │   └── user_repository.rs
│   ├── models/                   # Database models
│   │   ├── mod.rs
│   │   ├── analysis.rs
│   │   ├── sample.rs
│   │   ├── assembly.rs
│   │   ├── annotation.rs
│   │   ├── job.rs
│   │   └── user.rs
│   ├── migrations/               # Database migrations
│   │   ├── 001_initial_schema.sql
│   │   ├── 002_add_jobs_table.sql
│   │   └── 003_add_indexes.sql
│   └── schema.sql                # Complete schema definition
│
├── dto/                          # Data Transfer Objects
│   ├── mod.rs
│   ├── requests/                 # API request DTOs
│   │   ├── mod.rs
│   │   ├── analysis_requests.rs
│   │   ├── sample_requests.rs
│   │   ├── assembly_requests.rs
│   │   ├── annotation_requests.rs
│   │   ├── job_requests.rs
│   │   └── user_requests.rs
│   ├── responses/                # API response DTOs
│   │   ├── mod.rs
│   │   ├── analysis_responses.rs
│   │   ├── sample_responses.rs
│   │   ├── assembly_responses.rs
│   │   ├── annotation_responses.rs
│   │   ├── job_responses.rs
│   │   └── user_responses.rs
│   └── internal/                 # Internal service DTOs
│       ├── mod.rs
│       ├── pipeline_dto.rs
│       └── job_dto.rs
│
├── jobs/                         # Background job processing
│   ├── mod.rs
│   ├── queue.rs                  # Job queue implementation
│   ├── executor.rs               # Job execution engine
│   ├── scheduler.rs              # Cron-like job scheduler
│   ├── handlers/                 # Job-specific handlers
│   │   ├── mod.rs
│   │   ├── analysis_jobs.rs      # Analysis pipeline jobs
│   │   ├── file_processing_jobs.rs
│   │   ├── maintenance_jobs.rs
│   │   └── notification_jobs.rs
│   └── retry.rs                  # Retry logic and backoff
│
├── auth/                         # Authentication and authorization
│   ├── mod.rs
│   ├── jwt.rs                    # JWT token handling
│   ├── password.rs               # Password hashing/verification
│   ├── roles.rs                  # Role-based access control
│   └── permissions.rs            # Permission definitions
│
├── storage/                      # File storage abstraction
│   ├── mod.rs
│   ├── local.rs                  # Local filesystem storage
│   ├── s3.rs                     # AWS S3 storage
│   ├── traits.rs                 # Storage traits
│   └── manager.rs                # Storage tier management
│
├── monitoring/                   # Observability components
│   ├── mod.rs
│   ├── metrics.rs                # Prometheus metrics
│   ├── tracing.rs                # Distributed tracing
│   ├── health.rs                 # Health check implementations
│   └── alerts.rs                 # Alerting logic
│
├── config/                       # Configuration management
│   ├── mod.rs
│   ├── settings.rs               # Settings structure
│   ├── database.rs               # Database configuration
│   ├── server.rs                 # Server configuration
│   ├── auth.rs                   # Auth configuration
│   ├── storage.rs                # Storage configuration
│   └── monitoring.rs             # Monitoring configuration
│
├── validation/                   # Input validation
│   ├── mod.rs
│   ├── analysis_validation.rs
│   ├── sample_validation.rs
│   ├── user_validation.rs
│   └── common.rs                 # Common validation rules
│
├── error/                        # Error handling
│   ├── mod.rs
│   ├── api_error.rs              # API-specific errors
│   ├── service_error.rs          # Service-layer errors
│   ├── database_error.rs         # Database errors
│   └── job_error.rs              # Job processing errors
│
├── utils/                        # Utility functions
│   ├── mod.rs
│   ├── pagination.rs             # Pagination utilities
│   ├── time.rs                   # Time/date utilities
│   ├── crypto.rs                 # Cryptographic utilities
│   ├── file_utils.rs             # File handling utilities
│   └── test_utils.rs             # Testing utilities
│
├── integration/                  # External service integrations
│   ├── mod.rs
│   ├── email.rs                  # Email service integration
│   ├── webhook.rs                # Webhook delivery
│   └── notification.rs           # Notification service
│
├── main.rs                       # Application entry point
└── lib.rs                        # Library root

tests/                            # Integration tests
├── common/                       # Test utilities
│   ├── mod.rs
│   ├── fixtures.rs               # Test fixtures
│   ├── database.rs               # Test database setup
│   └── client.rs                 # Test HTTP client
├── api/                          # API endpoint tests
│   ├── analysis_tests.rs
│   ├── sample_tests.rs
│   ├── user_tests.rs
│   └── auth_tests.rs
├── services/                     # Service layer tests
│   ├── analysis_service_tests.rs
│   └── user_service_tests.rs
└── integration_tests.rs          # Full integration tests

benches/                          # Performance benchmarks
├── api_benchmarks.rs
├── database_benchmarks.rs
└── pipeline_benchmarks.rs

docs/                             # Documentation
├── api/                          # API documentation
│   ├── openapi.yaml              # OpenAPI specification
│   ├── authentication.md         # Authentication guide
│   ├── rate_limiting.md          # Rate limiting info
│   └── examples/                 # Usage examples
│       ├── curl_examples.sh
│       ├── python_client.py
│       └── javascript_client.js
├── deployment/                   # Deployment guides
│   ├── docker.md
│   ├── kubernetes.md
│   └── monitoring.md
├── development/                  # Development guides
│   ├── setup.md
│   ├── testing.md
│   ├── database.md
│   └── contributing.md
└── architecture/                 # Architecture documentation
    ├── system_overview.md
    ├── data_flow.md
    └── security.md

config/                           # Configuration files
├── default.toml                  # Default API configuration
├── development.toml              # Development overrides
├── production.toml               # Production configuration
└── test.toml                     # Test configuration

scripts/                          # Utility scripts
├── setup_dev.sh                  # Development environment setup
├── run_tests.sh                  # Test runner script
├── build_docker.sh               # Docker build script
├── deploy.sh                     # Deployment script
└── migrate.sh                    # Database migration script

migrations/                       # Database migrations (SQLx)
├── 20240101000000_initial_schema.sql
├── 20240102000000_add_jobs_table.sql
└── 20240103000000_add_indexes.sql

.github/                          # GitHub workflows
└── workflows/
    ├── ci.yml                    # Continuous integration
    ├── cd.yml                    # Continuous deployment
    └── security.yml              # Security scanning

docker/                           # Docker-related files
├── Dockerfile                    # Production image
├── Dockerfile.dev                # Development image
├── docker-compose.yml            # Development stack
└── docker-compose.prod.yml       # Production stack

k8s/                             # Kubernetes manifests
├── namespace.yaml
├── deployment.yaml
├── service.yaml
├── ingress.yaml
├── configmap.yaml
└── secrets.yaml
```

## File Organization Principles

### 1. Layer Separation
Each major architectural layer has its own directory with clear responsibilities:
- **API Layer**: HTTP handling, routing, middleware
- **Service Layer**: Business logic and orchestration  
- **Data Layer**: Database access and data models
- **Job Layer**: Background processing and scheduling

### 2. Domain-Driven Structure
Within each layer, files are organized by domain (analysis, samples, users, etc.) rather than by technical concerns.

### 3. Trait-Based Testing
Service traits in `src/services/traits/` enable easy mocking and testing by providing interfaces that can be implemented by both real services and test doubles.

### 4. Configuration Management
Environment-specific configuration files allow easy deployment across different environments while maintaining security through environment variables for secrets.

### 5. Documentation Co-location
API documentation lives alongside code to ensure it stays current, with examples and guides organized by audience (developers, operators, users).

## Key Design Decisions

### Module Organization
- **Flat hierarchy within domains**: Avoid deep nesting that makes navigation difficult
- **Clear module boundaries**: Each module has well-defined responsibilities
- **Consistent naming**: Use domain-specific naming consistently across layers

### Testing Strategy
- **Unit tests**: Co-located with source code using `#[cfg(test)]` modules
- **Integration tests**: Separate `tests/` directory for cross-service testing
- **Test utilities**: Shared fixtures and utilities in `tests/common/`

### Configuration Approach
- **Layered configuration**: Default -> Environment -> CLI arguments
- **Type-safe configuration**: Rust structs with serde for configuration parsing
- **Environment-specific overrides**: Different files for dev/staging/production

### Error Handling
- **Structured errors**: Different error types for different layers
- **Error context**: Rich error information for debugging
- **API error mapping**: Clean conversion from internal errors to HTTP responses

This structure provides a solid foundation for the REST API while maintaining clear separation of concerns and enabling easy testing and maintenance.
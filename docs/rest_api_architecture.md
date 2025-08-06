# Metagenomics REST API - System Architecture Design

## System Overview

This document outlines the comprehensive REST API architecture for the MetaForge metagenomics analysis pipeline. The API provides programmatic access to all pipeline capabilities while maintaining scalability, security, and maintainability.

## Architecture Principles

### Core Design Principles
- **Separation of Concerns**: Clear boundaries between API, business logic, and data layers
- **Scalability**: Horizontal scaling support with stateless design
- **Security**: Authentication, authorization, and input validation at all levels
- **Observability**: Comprehensive logging, metrics, and tracing
- **Performance**: Async processing with background job management
- **Reliability**: Circuit breakers, retries, and graceful degradation

### Technology Stack
- **Framework**: Axum (Rust async web framework)
- **Database**: SQLite/PostgreSQL with SQLx for async queries
- **Authentication**: JWT with refresh tokens
- **Background Jobs**: Tokio tasks with job queue
- **Serialization**: Serde JSON with custom validation
- **Documentation**: OpenAPI 3.0 with automated generation

## System Components

### 1. API Layer (`src/api/`)

#### Routes Module (`src/api/routes/`)
```rust
// Route organization by domain
mod analysis;     // Analysis pipeline endpoints
mod samples;      // Sample management
mod assemblies;   // Assembly operations  
mod annotations;  // Taxonomic classification
mod jobs;         // Background job management
mod admin;        // Administrative operations
mod health;       // Health checks and metrics
```

#### Handlers Module (`src/api/handlers/`)
```rust
// Handler functions for each route group
mod analysis_handlers;
mod sample_handlers;
mod assembly_handlers;
mod annotation_handlers;
mod job_handlers;
mod admin_handlers;
```

#### Middleware Stack (`src/api/middleware/`)
```rust
mod auth;           // JWT authentication
mod cors;           // CORS configuration
mod logging;        // Request/response logging
mod rate_limiting;  // Rate limiting per client
mod validation;     // Input validation
mod error_handling; // Centralized error handling
mod metrics;        // Prometheus metrics collection
```

### 2. Service Layer (`src/services/`)

#### Core Services
```rust
mod analysis_service;     // Analysis pipeline orchestration
mod sample_service;       // Sample lifecycle management
mod assembly_service;     // Assembly operations
mod annotation_service;   // Taxonomic classification
mod job_service;          // Background job management
mod user_service;         // User management
mod file_service;         // File upload/download handling
```

#### Service Traits for Testability
```rust
#[async_trait]
pub trait AnalysisService: Send + Sync {
    async fn create_analysis(&self, request: CreateAnalysisRequest) -> Result<Analysis>;
    async fn get_analysis(&self, id: Uuid) -> Result<Option<Analysis>>;
    async fn list_analyses(&self, filter: AnalysisFilter) -> Result<PaginatedResponse<Analysis>>;
}
```

### 3. Data Layer (`src/data/`)

#### Repository Pattern
```rust
mod repositories/
    analysis_repository.rs   // Analysis data access
    sample_repository.rs     // Sample data access
    user_repository.rs       // User data access
    job_repository.rs        // Job queue data access
```

#### Models and DTOs
```rust
mod models/              // Database models
    analysis.rs
    sample.rs
    user.rs
    job.rs

mod dto/                 // Data Transfer Objects
    requests/            // API request DTOs
    responses/           // API response DTOs
    internal/            // Internal service DTOs
```

### 4. Background Processing (`src/jobs/`)

#### Job System Architecture
```rust
mod queue;              // Job queue implementation
mod executor;           // Job execution engine
mod handlers/           // Job-specific handlers
    analysis_jobs.rs    // Analysis pipeline jobs
    cleanup_jobs.rs     // Cleanup and maintenance
mod scheduler;          // Cron-like job scheduling
```

### 5. Configuration (`src/config/`)

#### Environment-based Configuration
```rust
pub struct ApiConfig {
    pub server: ServerConfig,
    pub database: DatabaseConfig,
    pub auth: AuthConfig,
    pub pipeline: PipelineConfig,
    pub storage: StorageConfig,
    pub monitoring: MonitoringConfig,
}
```

## API Design Specifications

### 1. RESTful Endpoints

#### Analysis Pipeline Endpoints
```
POST   /api/v1/analyses                    # Create new analysis
GET    /api/v1/analyses                    # List analyses (with pagination/filtering)
GET    /api/v1/analyses/{id}               # Get specific analysis
DELETE /api/v1/analyses/{id}               # Delete analysis
GET    /api/v1/analyses/{id}/results       # Get analysis results
GET    /api/v1/analyses/{id}/status        # Get analysis status
POST   /api/v1/analyses/{id}/cancel        # Cancel running analysis
```

#### Sample Management Endpoints  
```
POST   /api/v1/samples                     # Upload sample files
GET    /api/v1/samples                     # List samples
GET    /api/v1/samples/{id}                # Get sample details
DELETE /api/v1/samples/{id}                # Delete sample
PUT    /api/v1/samples/{id}/metadata       # Update sample metadata
```

#### Assembly Endpoints
```
POST   /api/v1/assemblies                  # Create assembly from samples
GET    /api/v1/assemblies                  # List assemblies
GET    /api/v1/assemblies/{id}             # Get assembly details
GET    /api/v1/assemblies/{id}/contigs     # Get assembly contigs
GET    /api/v1/assemblies/{id}/stats       # Get assembly statistics
```

#### Background Job Endpoints
```
GET    /api/v1/jobs                        # List user's jobs
GET    /api/v1/jobs/{id}                   # Get job status/results
DELETE /api/v1/jobs/{id}                   # Cancel job
```

### 2. Authentication & Authorization

#### JWT-based Authentication
```rust
pub struct Claims {
    pub sub: Uuid,           // User ID
    pub role: UserRole,      // User role
    pub exp: usize,          // Expiration
    pub iat: usize,          // Issued at
}

#[derive(Serialize, Deserialize)]
pub enum UserRole {
    Admin,
    Researcher,
    ReadOnly,
}
```

#### Role-based Access Control
- **Admin**: Full system access, user management
- **Researcher**: Create/manage own analyses, upload samples
- **ReadOnly**: View shared analyses and results

### 3. Input Validation & Error Handling

#### Request Validation
```rust
#[derive(Deserialize, Validate)]
pub struct CreateAnalysisRequest {
    #[validate(length(min = 1, max = 100))]
    pub name: String,
    
    #[validate(custom = "validate_sample_ids")]
    pub sample_ids: Vec<Uuid>,
    
    pub analysis_mode: AnalysisMode,
    
    #[validate(range(min = 1, max = 64))]
    pub max_threads: Option<u32>,
}
```

#### Standardized Error Responses
```rust
#[derive(Serialize)]
pub struct ApiError {
    pub error_code: String,
    pub message: String,
    pub details: Option<serde_json::Value>,
    pub request_id: String,
}
```

### 4. Response Formats

#### Standard Response Wrapper
```rust
#[derive(Serialize)]
pub struct ApiResponse<T> {
    pub success: bool,
    pub data: Option<T>,
    pub error: Option<ApiError>,
    pub pagination: Option<PaginationInfo>,
}
```

#### Pagination Support
```rust
#[derive(Serialize, Deserialize)]
pub struct PaginationInfo {
    pub page: u32,
    pub per_page: u32,
    pub total_pages: u32,
    pub total_items: u64,
}
```

## Database Schema Design

### Core Tables

#### Users Table
```sql
CREATE TABLE users (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    email VARCHAR(255) UNIQUE NOT NULL,
    password_hash VARCHAR(255) NOT NULL,
    role VARCHAR(50) NOT NULL DEFAULT 'researcher',
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    last_login TIMESTAMPTZ
);
```

#### Samples Table
```sql
CREATE TABLE samples (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID REFERENCES users(id) ON DELETE CASCADE,
    name VARCHAR(255) NOT NULL,
    file_path VARCHAR(500) NOT NULL,
    file_size BIGINT NOT NULL,
    file_format VARCHAR(50) NOT NULL,
    metadata JSONB DEFAULT '{}',
    upload_status VARCHAR(50) DEFAULT 'pending',
    created_at TIMESTAMPTZ DEFAULT NOW()
);
```

#### Analyses Table
```sql
CREATE TABLE analyses (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID REFERENCES users(id) ON DELETE CASCADE,
    name VARCHAR(255) NOT NULL,
    analysis_mode VARCHAR(50) NOT NULL,
    status VARCHAR(50) DEFAULT 'pending',
    config JSONB NOT NULL,
    started_at TIMESTAMPTZ,
    completed_at TIMESTAMPTZ,
    error_message TEXT,
    created_at TIMESTAMPTZ DEFAULT NOW()
);
```

#### Analysis Samples (Junction Table)
```sql
CREATE TABLE analysis_samples (
    analysis_id UUID REFERENCES analyses(id) ON DELETE CASCADE,
    sample_id UUID REFERENCES samples(id) ON DELETE CASCADE,
    PRIMARY KEY (analysis_id, sample_id)
);
```

#### Jobs Table
```sql
CREATE TABLE jobs (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    job_type VARCHAR(100) NOT NULL,
    status VARCHAR(50) DEFAULT 'pending',
    payload JSONB NOT NULL,
    result JSONB,
    error_message TEXT,
    attempts INTEGER DEFAULT 0,
    max_attempts INTEGER DEFAULT 3,
    scheduled_at TIMESTAMPTZ DEFAULT NOW(),
    started_at TIMESTAMPTZ,
    completed_at TIMESTAMPTZ,
    created_at TIMESTAMPTZ DEFAULT NOW()
);
```

## Background Job Processing

### Job Types
```rust
#[derive(Serialize, Deserialize)]
pub enum JobType {
    AnalysisPipeline,
    FileProcessing,
    DatabaseMaintenance,
    ReportGeneration,
    DataCleanup,
}
```

### Job Execution Flow
1. **Job Creation**: API endpoints create jobs and store in database
2. **Job Scheduling**: Background worker polls for pending jobs
3. **Job Execution**: Workers execute jobs with proper error handling
4. **Progress Updates**: Jobs update status and progress in database
5. **Result Storage**: Completed jobs store results and notify clients

### Job Retry Strategy
```rust
pub struct RetryConfig {
    pub max_attempts: u32,
    pub backoff_strategy: BackoffStrategy,
    pub retry_on: Vec<ErrorType>,
}

pub enum BackoffStrategy {
    Fixed(Duration),
    Exponential { base: Duration, max: Duration },
    Linear(Duration),
}
```

## File Storage Architecture

### Multi-tier Storage Strategy
1. **Hot Storage**: Recent uploads and active analyses (Local SSD)
2. **Warm Storage**: Completed analyses (Network attached storage)
3. **Cold Storage**: Archived results (Object storage like S3)

### File Management Service
```rust
pub trait FileStorageService: Send + Sync {
    async fn upload_file(&self, file: FileUpload) -> Result<FileMetadata>;
    async fn download_file(&self, file_id: &str) -> Result<FileStream>;
    async fn delete_file(&self, file_id: &str) -> Result<()>;
    async fn move_to_tier(&self, file_id: &str, tier: StorageTier) -> Result<()>;
}
```

## Security Considerations

### Input Sanitization
- All file uploads scanned for malicious content
- Input validation on all API endpoints
- SQL injection prevention through parameterized queries
- XSS prevention through proper encoding

### Rate Limiting
```rust
pub struct RateLimitConfig {
    pub requests_per_minute: u32,
    pub requests_per_hour: u32,
    pub burst_allowance: u32,
    pub whitelist: Vec<IpAddr>,
}
```

### API Key Management
```rust
pub struct ApiKey {
    pub id: Uuid,
    pub key_hash: String,
    pub user_id: Uuid,
    pub permissions: Vec<Permission>,
    pub expires_at: Option<DateTime<Utc>>,
    pub created_at: DateTime<Utc>,
}
```

## Monitoring and Observability

### Metrics Collection
- Request/response metrics (latency, throughput, errors)
- Business metrics (analyses per day, user activity)
- System metrics (CPU, memory, disk usage)
- Custom pipeline metrics (assembly time, contig count)

### Logging Strategy
```rust
pub struct RequestLogger {
    pub log_level: LogLevel,
    pub include_request_body: bool,
    pub include_response_body: bool,
    pub mask_sensitive_data: bool,
}
```

### Health Checks
```rust
pub struct HealthCheckService {
    pub database_check: DatabaseHealthCheck,
    pub storage_check: StorageHealthCheck,
    pub pipeline_check: PipelineHealthCheck,
    pub external_services: Vec<ExternalServiceCheck>,
}
```

## Testing Strategy

### Test Categories
1. **Unit Tests**: Individual functions and methods
2. **Integration Tests**: Service layer interactions
3. **API Tests**: HTTP endpoint testing
4. **Load Tests**: Performance under load
5. **Contract Tests**: API contract validation

### Test Infrastructure
```rust
pub struct TestFixture {
    pub test_db: TestDatabase,
    pub mock_services: MockServiceContainer,
    pub test_client: TestHttpClient,
    pub test_data: TestDataGenerator,
}
```

## Deployment Architecture

### Container Strategy
```dockerfile
# Multi-stage build for optimized production image
FROM rust:1.70 AS builder
WORKDIR /app
COPY . .
RUN cargo build --release

FROM debian:bookworm-slim
RUN apt-get update && apt-get install -y ca-certificates
COPY --from=builder /app/target/release/meta-api /usr/local/bin/
EXPOSE 8080
CMD ["meta-api"]
```

### Environment Configuration
```yaml
# docker-compose.yml
services:
  api:
    image: meta-api:latest
    ports:
      - "8080:8080"
    environment:
      - DATABASE_URL=${DATABASE_URL}
      - JWT_SECRET=${JWT_SECRET}
      - LOG_LEVEL=info
    depends_on:
      - database
      - redis
  
  database:
    image: postgres:15
    environment:
      POSTGRES_DB: metagenomics
      POSTGRES_USER: ${DB_USER}
      POSTGRES_PASSWORD: ${DB_PASSWORD}
    volumes:
      - postgres_data:/var/lib/postgresql/data
```

## Performance Considerations

### Caching Strategy
- **Application-level**: In-memory caching for frequently accessed data
- **Database-level**: Query result caching
- **HTTP-level**: Response caching for static data

### Database Optimization
- Connection pooling with configurable limits
- Read replicas for query scaling
- Proper indexing strategy for common queries
- Query optimization and monitoring

### Async Processing
- Non-blocking I/O for all operations
- Background job processing for long-running tasks
- Streaming responses for large datasets
- Connection reuse and HTTP/2 support

## API Documentation

### OpenAPI Integration
- Automated schema generation from Rust types
- Interactive API documentation with Swagger UI
- Client SDK generation for multiple languages
- Contract testing based on OpenAPI specs

### Documentation Structure
```
/docs/
  ├── api/
  │   ├── openapi.yaml          # OpenAPI specification
  │   ├── authentication.md     # Auth documentation
  │   └── examples/            # Usage examples
  ├── deployment/
  │   ├── docker.md            # Container deployment
  │   └── kubernetes.md        # K8s deployment
  └── development/
      ├── setup.md             # Development setup
      └── testing.md           # Testing guide
```

This architecture provides a solid foundation for building a scalable, maintainable, and secure REST API for the metagenomics pipeline. The design emphasizes separation of concerns, testability, and operational excellence while leveraging Rust's performance and safety benefits.
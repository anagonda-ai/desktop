# Multi-stage Docker build for Plant MGC Analysis Pipeline
FROM python:3.11-slim as builder

# Set build arguments
ARG DEBIAN_FRONTEND=noninteractive

# Install system dependencies for building
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    wget \
    git \
    && rm -rf /var/lib/apt/lists/*

# Create virtual environment
RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Copy requirements and install Python dependencies
COPY pyproject.toml /tmp/
RUN pip install --upgrade pip setuptools wheel && \
    pip install "/tmp/[dev]"

# Production stage
FROM python:3.11-slim

# Set environment variables
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PATH="/opt/venv/bin:$PATH" \
    DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    ncbi-blast+ \
    hmmer \
    curl \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Create application user
RUN useradd --create-home --shell /bin/bash app

# Copy virtual environment from builder
COPY --from=builder /opt/venv /opt/venv

# Set working directory
WORKDIR /app

# Copy application code
COPY src/ /app/src/
COPY pyproject.toml /app/
COPY README.md /app/

# Install the package
RUN pip install -e .

# Create necessary directories
RUN mkdir -p /app/data /app/output /app/cache && \
    chown -R app:app /app

# Switch to application user
USER app

# Set default environment variables
ENV PLANT_MGC_DATA_DIR=/app/data \
    PLANT_MGC_OUTPUT_DIR=/app/output \
    PLANT_MGC_CACHE_DIR=/app/cache

# Expose ports (if needed for web interface)
EXPOSE 8000

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD plant-mgc --help || exit 1

# Default command
ENTRYPOINT ["plant-mgc"]
CMD ["--help"]
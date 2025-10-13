#!/bin/bash
# Download GTDB Taxonomy Database for MetaForge
#
# This script downloads the GTDB bacterial taxonomy database (Release 220+)
# and prepares it for import into MetaForge's SQLite database.
#
# Usage: ./scripts/download_taxonomy.sh

set -e  # Exit on error
set -u  # Exit on undefined variable

# Configuration
TAXONOMY_DIR="data/taxonomy"
GTDB_BASE_URL="https://data.gtdb.ecogenomic.org/releases/latest"
TAXONOMY_FILE="bac120_taxonomy.tsv"

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${BLUE}═══════════════════════════════════════════════════════${NC}"
echo -e "${BLUE}  MetaForge Taxonomy Database Downloader${NC}"
echo -e "${BLUE}═══════════════════════════════════════════════════════${NC}"
echo

# Create taxonomy directory
echo -e "${YELLOW}📁 Creating taxonomy directory...${NC}"
mkdir -p "$TAXONOMY_DIR"

# Download taxonomy file
echo -e "${YELLOW}📥 Downloading GTDB bacterial taxonomy...${NC}"
echo "   Source: $GTDB_BASE_URL/$TAXONOMY_FILE.gz"
echo "   Destination: $TAXONOMY_DIR/$TAXONOMY_FILE"
echo

if [ -f "$TAXONOMY_DIR/$TAXONOMY_FILE" ]; then
    echo -e "${YELLOW}⚠️  Taxonomy file already exists. Backup and re-download? (y/n)${NC}"
    read -r response
    if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]; then
        mv "$TAXONOMY_DIR/$TAXONOMY_FILE" "$TAXONOMY_DIR/${TAXONOMY_FILE}.backup.$(date +%Y%m%d_%H%M%S)"
        echo "   Backed up existing file"
    else
        echo "   Using existing file"
        exit 0
    fi
fi

# Download with wget or curl
if command -v wget &> /dev/null; then
    wget --progress=bar:force:noscroll \
         -O "$TAXONOMY_DIR/$TAXONOMY_FILE.gz" \
         "$GTDB_BASE_URL/$TAXONOMY_FILE.gz"
elif command -v curl &> /dev/null; then
    curl -# \
         -o "$TAXONOMY_DIR/$TAXONOMY_FILE.gz" \
         "$GTDB_BASE_URL/$TAXONOMY_FILE.gz"
else
    echo "Error: Neither wget nor curl found. Please install one of them."
    exit 1
fi

# Extract
echo
echo -e "${YELLOW}📦 Extracting taxonomy file...${NC}"
gunzip -f "$TAXONOMY_DIR/$TAXONOMY_FILE.gz"

# Verify
if [ ! -f "$TAXONOMY_DIR/$TAXONOMY_FILE" ]; then
    echo -e "${YELLOW}❌ Error: Failed to extract taxonomy file${NC}"
    exit 1
fi

FILE_SIZE=$(du -h "$TAXONOMY_DIR/$TAXONOMY_FILE" | cut -f1)
LINE_COUNT=$(wc -l < "$TAXONOMY_DIR/$TAXONOMY_FILE")

echo
echo -e "${GREEN}✅ Taxonomy download complete!${NC}"
echo
echo "   File: $TAXONOMY_DIR/$TAXONOMY_FILE"
echo "   Size: $FILE_SIZE"
echo "   Taxa: $LINE_COUNT entries"
echo

# Show sample
echo -e "${BLUE}Sample entries:${NC}"
head -n 3 "$TAXONOMY_DIR/$TAXONOMY_FILE"
echo "   ..."
echo

# Next steps
echo -e "${BLUE}═══════════════════════════════════════════════════════${NC}"
echo -e "${GREEN}Next Steps:${NC}"
echo
echo "1. Initialize database (if not already done):"
echo "   ${YELLOW}./target/release/meta-forge database init data/metagenomics.db${NC}"
echo
echo "2. Import taxonomy into database:"
echo "   ${YELLOW}./target/release/meta-forge database import-taxonomy \\"
echo "     data/metagenomics.db \\"
echo "     --taxonomy $TAXONOMY_DIR/$TAXONOMY_FILE${NC}"
echo
echo "3. Verify import:"
echo "   ${YELLOW}sqlite3 data/metagenomics.db 'SELECT COUNT(*) FROM taxonomy;'${NC}"
echo
echo -e "${BLUE}═══════════════════════════════════════════════════════${NC}"

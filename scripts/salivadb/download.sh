#!/bin/bash

BASE_URL="https://webs.iiitd.edu.in/raghava/salivadb"
FILES=(protein metabolite microbe mirna gene)

# Usage info
usage() {
    echo "Usage: $0 -d <download_dir>"
    echo "  -d  Directory to save files (required)"
    exit 1
}

# Parse arguments
while getopts "d:" opt; do
    case $opt in
        d) DOWNLOAD_DIR="$OPTARG" ;;
        *) usage ;;
    esac
done

# Validate
[ -z "$DOWNLOAD_DIR" ] && usage

# Create dir if it doesn't exist
mkdir -p "$DOWNLOAD_DIR"

for file in "${FILES[@]}"; do
    echo "Downloading ${file}.tsv..."
    curl -fSL "${BASE_URL}/${file}.tsv" -o "${DOWNLOAD_DIR}/${file}.tsv"
done

echo "Done! Files saved to: $DOWNLOAD_DIR"

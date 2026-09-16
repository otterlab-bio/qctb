#!/bin/bash
# Build script for optimized release binary

set -e

echo "Building qctb in release mode..."

cargo build --release

echo "Build complete!"
echo "Binary: target/release/qctb"
echo ""
echo "Testing binary..."
./target/release/qctb --version
./target/release/qctb --help

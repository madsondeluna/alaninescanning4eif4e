#!/bin/bash
#
# Quick test script for Rosetta Alanine Scanning
#
# This script runs a minimal test to verify the installation
# and basic functionality.

echo "======================================================================"
echo "Rosetta Alanine Scanning - Quick Test"
echo "======================================================================"

# Check if virtual environment is active
if [[ -z "${VIRTUAL_ENV}" ]]; then
    echo "  Warning: No virtual environment detected"
    echo "   Consider activating your venv: source venv/bin/activate"
    echo ""
fi

# Check for ROSETTA environment variable
if [[ -z "${ROSETTA}" ]]; then
    echo " Error: ROSETTA environment variable not set"
    echo "   Set it with: export ROSETTA=/path/to/rosetta"
    exit 1
fi

echo " ROSETTA path: $ROSETTA"
echo ""

# Create test directory
TEST_DIR="test_run"
mkdir -p "$TEST_DIR"

# Check if example PDB exists
if [[ ! -f "protein.pdb" ]]; then
    echo "  No protein.pdb found. Please provide a test PDB file."
    echo ""
    echo "You can download a test structure:"
    echo "  wget https://files.rcsb.org/download/1ABC.pdb -O protein.pdb"
    exit 1
fi

echo "Running test workflow..."
echo ""

# Step 1: Generate mutations
echo "[1/3] Generating mutations..."
rosetta-scan scan protein.pdb \
    -c A \
    --output "$TEST_DIR/mutations" \
    --format rosetta

if [[ $? -ne 0 ]]; then
    echo " Mutation generation failed"
    exit 1
fi

echo " Mutations generated"
echo ""

# Step 2: Initialize config (optional)
echo "[2/3] Creating configuration..."
rosetta-scan init-config -o "$TEST_DIR/config.yaml"

# Edit config for test run (small nstruct)
sed -i.bak 's/nstruct: 35/nstruct: 2/' "$TEST_DIR/config.yaml"

echo " Configuration created"
echo ""

# Step 3: Display help
echo "[3/3] Checking CLI commands..."
rosetta-scan --help

echo ""
echo "======================================================================"
echo " Quick test completed successfully!"
echo "======================================================================"
echo ""
echo "To run a full analysis:"
echo ""
echo "  rosetta-scan pipeline protein.pdb \\"
echo "    --chains A B \\"
echo "    --interface-only \\"
echo "    --nstruct 35 \\"
echo "    --output results/"
echo ""
echo "Test files created in: $TEST_DIR/"
echo ""

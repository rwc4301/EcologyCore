#!/bin/bash
set -e

# Go to project root directory
cd "$(dirname "$0")"

# Check if virtual environment exists
if [ ! -d "venv" ]; then
    echo "Virtual environment not found. Setting it up first..."
    ./setup_env.sh
    echo "Now you can run this script again."
    exit 0
fi

# Activate virtual environment
echo "Activating virtual environment..."
source venv/bin/activate

# Create output directory
mkdir -p examples/plots

# Run the example
echo "Running basic workflow example..."
$PYTHON examples/basic_workflow.py

echo ""
echo "Example completed successfully! Check the 'examples/plots' directory for output visualizations."
#!/bin/bash
set -e

# Go to project root directory
cd "$(dirname "$0")"

# Use Python 3 explicitly
PYTHON="python3"

# Check if Python 3 is available
if ! command -v $PYTHON &> /dev/null; then
    echo "Python 3 is required but not found. Please install Python 3 and try again."
    exit 1
fi

echo "Setting up virtual environment for EcologyCore..."

# Create virtual environment if it doesn't exist
if [ ! -d "venv" ]; then
    echo "Creating virtual environment..."
    $PYTHON -m venv venv
else
    echo "Virtual environment already exists."
fi

# Activate virtual environment
echo "Activating virtual environment..."
source venv/bin/activate

# Upgrade pip
echo "Upgrading pip..."
pip install --upgrade pip

# Install dependencies
echo "Installing dependencies..."
pip install -r requirements.txt

# Install the package in development mode
echo "Installing EcologyCore in development mode..."
pip install -e .

echo ""
echo "Setup complete! To activate the environment, run:"
echo "source venv/bin/activate"
echo ""
echo "To run the example, use:"
echo "./run_example.sh"
#!/bin/bash
set -e

# Print commands for installation options
echo "EcologyCore Python Installation Script"
echo "======================================"
echo
echo "This script provides options to install EcologyCore using conda, pip, or uv."
echo

# Function to check if a command exists
command_exists() {
  command -v "$1" >/dev/null 2>&1
}

# Check for available package managers
CONDA_AVAILABLE=false
UV_AVAILABLE=false
PIP_AVAILABLE=false

if command_exists conda; then
  CONDA_AVAILABLE=true
fi

if command_exists uv; then
  UV_AVAILABLE=true
fi

if command_exists pip; then
  PIP_AVAILABLE=true
fi

# Display installation options
echo "Available installation methods:"
if $CONDA_AVAILABLE; then
  echo "1) conda (recommended for scientific packages)"
fi
if $UV_AVAILABLE; then
  echo "2) uv (fast modern Python package manager)"
fi
if $PIP_AVAILABLE; then
  echo "3) pip (standard Python package manager)"
fi
echo "q) Quit installation"
echo

# Get user choice
read -p "Select installation method (1-3 or q): " choice

case $choice in
  1)
    if ! $CONDA_AVAILABLE; then
      echo "Error: conda is not installed."
      exit 1
    fi
    
    echo "Installing with conda..."
    echo "Creating conda environment from environment.yml..."
    conda env create -f environment.yml
    
    echo
    echo "Installation complete! To activate the environment, run:"
    echo "conda activate ecology-core"
    ;;
    
  2)
    if ! $UV_AVAILABLE; then
      echo "Error: uv is not installed."
      exit 1
    fi
    
    echo "Installing with uv..."
    echo "Creating virtual environment..."
    uv venv
    
    echo "Activating virtual environment..."
    source .venv/bin/activate
    
    echo "Installing dependencies..."
    uv pip install -r requirements.txt
    
    echo "Installing package in development mode..."
    uv pip install -e .
    
    echo
    echo "Installation complete! To activate the environment, run:"
    echo "source .venv/bin/activate"
    ;;
    
  3)
    if ! $PIP_AVAILABLE; then
      echo "Error: pip is not installed."
      exit 1
    fi
    
    echo "Installing with pip..."
    echo "Creating virtual environment..."
    python -m venv venv
    
    echo "Activating virtual environment..."
    source venv/bin/activate
    
    echo "Installing dependencies..."
    pip install -r requirements.txt
    
    echo "Installing package in development mode..."
    pip install -e .
    
    echo
    echo "Installation complete! To activate the environment, run:"
    echo "source venv/bin/activate"
    ;;
    
  q|Q)
    echo "Installation cancelled."
    exit 0
    ;;
    
  *)
    echo "Invalid option. Installation cancelled."
    exit 1
    ;;
esac
#!/bin/bash

# Start GUI script for DISPLATH
set -e

echo "🚀 Starting DISPLATH Web GUI..."

# Change to GUI directory
cd "$(dirname "$0")"

# Create logs directory
mkdir -p logs

# Check if conda environment exists
if ! conda env list | grep -q "displath-gui"; then
    echo "❌ Conda environment 'displath-gui' not found!"
    echo "Creating environment from ../environment.yml..."
    conda env create -f ../environment.yml
fi

# Activate conda environment
echo "🔧 Activating conda environment..."
source /beegfs/home/xuke/Softwares/anaconda3/etc/profile.d/conda.sh
conda activate displath-gui

# Install/update Python dependencies
echo "📦 Installing/updating Python dependencies..."
pip install -r requirements.txt

# Check Julia availability
echo "🔍 Checking Julia availability..."
if ! which julia >/dev/null 2>&1; then
    echo "⚠️  Julia not found in PATH. Simulations will not work."
    echo "Please install Julia and ensure it's in your PATH."
else
    julia_version=$(julia --version)
    echo "✅ Found Julia: $julia_version"
fi

# Check DISPLATH source files
DISPLATH_MAIN="../../src/main.jl"
if [ ! -f "$DISPLATH_MAIN" ]; then
    echo "⚠️  DISPLATH source not found at $DISPLATH_MAIN"
    echo "GUI will start but simulations may fail."
else
    echo "✅ DISPLATH source found"
fi

# Start the Flask server
echo ""
echo "🌐 Starting DISPLATH GUI server..."
echo "📍 URL: http://localhost:5000"
echo "📝 Logs: logs/"
echo "⌨️  Press Ctrl+C to stop the server"
echo ""

python3 server.py
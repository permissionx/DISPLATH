#!/bin/bash

# ARCS Installation Script
# This script sets up the ARCS environment

set -e  # Exit on any error

echo "Setting up ARCS environment..."

# Create the main ARCS repository directory
ARCS_REPO="$HOME/.arcs_repository"
echo "Creating ARCS repository directory: $ARCS_REPO"
mkdir -p "$ARCS_REPO"

# Create subdirectories
echo "Creating dte_repository directory..."
mkdir -p "$ARCS_REPO/dte_repository"

echo "Creating thetatau_repository directory..."
mkdir -p "$ARCS_REPO/thetatau_repository"

# Set up environment variable
echo "Setting up environment variable..."

# Add to .bashrc if it doesn't already exist
if ! grep -q "# ARCS Environment" "$HOME/.bashrc"; then
    echo "" >> "$HOME/.bashrc"
    echo "# ARCS Environment" >> "$HOME/.bashrc"
    echo "export ARCS_REPO=$ARCS_REPO" >> "$HOME/.bashrc"
    echo "export ARCS_HOME=/beegfs/science-share/arcs/DISPLATH" >> "$HOME/.bashrc"
    echo "ARCS environment variables added to .bashrc"
else
    echo "ARCS environment variables already exist in .bashrc"
fi

# Add Julia environment setup to .bashrc if it doesn't already exist
if ! grep -q "juliaup initialize" "$HOME/.bashrc"; then
    echo "" >> "$HOME/.bashrc"
    echo "# >>> juliaup initialize >>>" >> "$HOME/.bashrc"
    echo "" >> "$HOME/.bashrc"
    echo "# !! Contents within this block are managed by juliaup !!" >> "$HOME/.bashrc"
    echo "" >> "$HOME/.bashrc"
    echo 'case ":$PATH:" in' >> "$HOME/.bashrc"
    echo '    *:/beegfs/science-share/julia/bin:*)' >> "$HOME/.bashrc"
    echo '        ;;' >> "$HOME/.bashrc"
    echo "" >> "$HOME/.bashrc"
    echo '    *)' >> "$HOME/.bashrc"
    echo '        export PATH=/beegfs/science-share/julia/bin${PATH:+:${PATH}}' >> "$HOME/.bashrc"
    echo '        ;;' >> "$HOME/.bashrc"
    echo 'esac' >> "$HOME/.bashrc"
    echo "" >> "$HOME/.bashrc"
    echo "# <<< juliaup initialize <<<" >> "$HOME/.bashrc"
    echo "Julia environment setup added to .bashrc"
else
    echo "Julia environment setup already exists in .bashrc"
fi

# Also add to .bash_profile if it exists
if [ -f "$HOME/.bash_profile" ]; then
    if ! grep -q "# ARCS Environment" "$HOME/.bash_profile"; then
        echo "" >> "$HOME/.bash_profile"
        echo "# ARCS Environment" >> "$HOME/.bash_profile"
        echo "export ARCS_REPO=$ARCS_REPO" >> "$HOME/.bash_profile"
        echo "export ARCS_HOME=/beegfs/science-share/arcs/DISPLATH" >> "$HOME/.bash_profile"
        echo "ARCS environment variables added to .bash_profile"
    fi
    
    # Add Julia environment setup to .bash_profile if it doesn't already exist
    if ! grep -q "juliaup initialize" "$HOME/.bash_profile"; then
        echo "" >> "$HOME/.bash_profile"
        echo "# >>> juliaup initialize >>>" >> "$HOME/.bash_profile"
        echo "" >> "$HOME/.bash_profile"
        echo "# !! Contents within this block are managed by juliaup !!" >> "$HOME/.bash_profile"
        echo "" >> "$HOME/.bash_profile"
        echo 'case ":$PATH:" in' >> "$HOME/.bash_profile"
        echo '    *:/beegfs/science-share/julia/bin:*)' >> "$HOME/.bash_profile"
        echo '        ;;' >> "$HOME/.bash_profile"
        echo "" >> "$HOME/.bash_profile"
        echo '    *)' >> "$HOME/.bash_profile"
        echo '        export PATH=/beegfs/science-share/julia/bin${PATH:+:${PATH}}' >> "$HOME/.bash_profile"
        echo '        ;;' >> "$HOME/.bash_profile"
        echo 'esac' >> "$HOME/.bash_profile"
        echo "" >> "$HOME/.bash_profile"
        echo "# <<< juliaup initialize <<<" >> "$HOME/.bash_profile"
        echo "Julia environment setup added to .bash_profile"
    fi
fi

# Set environment variables for current session
export ARCS_REPO="$ARCS_REPO"
export ARCS_HOME="/beegfs/science-share/arcs/DISPLATH"

# Set Julia PATH for current session
case ":$PATH:" in
    *:/beegfs/science-share/julia/bin:*)
        ;;
    *)
        export PATH=/beegfs/science-share/julia/bin${PATH:+:${PATH}}
        ;;
esac

echo ""
echo "ARCS installation completed successfully!"
echo ""
echo "Directory structure created:"
echo "  $ARCS_REPO/"
echo "  ├── dte_repository/"
echo "  └── thetatau_repository/"
echo ""
echo "Environment variables set:"
echo "  ARCS_REPO=$ARCS_REPO"
echo "  ARCS_HOME=/beegfs/science-share/arcs/DISPLATH"
echo "  Julia PATH=/beegfs/science-share/julia/bin (added to PATH)"
echo ""
echo "Examples can be found at: /beegfs/science-share/arcs/DISPLATH/examples"
echo ""
echo "Please restart your terminal or run 'source ~/.bashrc' to load the environment variables." 
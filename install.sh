#!/usr/bin/env bash

set -euo pipefail

###############################################################################
# Choose environment type
###############################################################################

USE_CONDA=0
USE_VENV=0

if command -v conda >/dev/null 2>&1; then

    echo
    echo "Conda detected."
    echo
    echo "Choose environment type:"
    echo "  1) Conda (recommended)"
    echo "  2) venv"
    echo

    while true; do
        printf "Selection [1/2] (default=1): "
        read -r CHOICE
        CHOICE="${CHOICE:-1}"
        case "$CHOICE" in
            1)
                USE_CONDA=1
                break
                ;;
            2)
                USE_VENV=1
                break
                ;;
            *)
                echo "Please enter 1 or 2."
                ;;
        esac
    done

else

    echo
    echo "Conda was not found."
    echo

    printf "Use a Python venv instead? [y/N] "
    read -r USE_VENV_REPLY

    if [[ "$USE_VENV_REPLY" =~ ^([Yy]|[Yy][Ee][Ss])$ ]]; then

        USE_VENV=1

    else

        echo
        echo "Please install Miniforge and rerun this script."
        echo
        echo "https://conda-forge.org/miniforge/"
        echo

        exit 1
    fi
fi

###############################################################################
# Environment creation
###############################################################################

ENV_NAME="$(basename "$PWD")"

if [[ "$USE_CONDA" -eq 1 ]]; then

    echo
    echo "Creating Conda environment '$ENV_NAME'..."

    eval "$(conda shell.bash hook)"

    conda create \
        -n "$ENV_NAME" \
        "python>=3.12,<3.14" \
        pip

    ENV_PREFIX="$(
        conda env list |
        awk -v env="$ENV_NAME" '$1 == env {print $NF}'
    )"

    if [[ -z "$ENV_PREFIX" ]]; then
        echo "ERROR: Environment creation failed."
        exit 1
    fi

    conda activate "$ENV_NAME"

    if [[ "${CONDA_DEFAULT_ENV:-}" != "$ENV_NAME" ]]; then
        echo "ERROR: Failed to activate Conda environment."
        exit 1
    fi

    PYTHON="python"

else

    VENV_DIR=".venv"

    if [[ -e "$VENV_DIR" ]]; then
        echo "ERROR: $VENV_DIR already exists."
        exit 1
    fi

    PYTHON_HOST=""

    for py in python3.13 python3.12; do
        if command -v "$py" >/dev/null 2>&1; then
            PYTHON_HOST="$py"
            break
        fi
    done

    if [[ -z "$PYTHON_HOST" ]]; then
        echo "ERROR: Could not find python3.12 or python3.13."
        exit 1
    fi

    echo
    echo "Creating venv..."

    "$PYTHON_HOST" -m venv "$VENV_DIR"

    if [[ ! -x "$VENV_DIR/bin/python" ]]; then
        echo "ERROR: Failed to create venv."
        exit 1
    fi

    PYTHON="$VENV_DIR/bin/python"
fi

###############################################################################
# Install local package
###############################################################################

echo
echo "Upgrading packaging tools..."

"$PYTHON" -m pip install --upgrade \
    pip \
    setuptools \
    wheel

echo
echo "Installing current package..."

"$PYTHON" -m pip install .

###############################################################################
# Optional SciGUI installation
###############################################################################

printf "\nInstall SciGUI? ([y]/n) "
read -r INSTALL_SCIGUI
INSTALL_SCIGUI=${INSTALL_SCIGUI:-y}
if [[ "$INSTALL_SCIGUI" =~ ^([Yy]|[Yy][Ee][Ss])$ ]]; then

    SCIGUI_TMPDIR="$(mktemp -d)"

    cleanup() {
        rm -rf "$SCIGUI_TMPDIR"
    }

    trap cleanup EXIT

    echo
    echo "Downloading SciGUI..."

    curl -fL \
        "https://github.com/times-software/SciGUI/archive/refs/heads/main.tar.gz" \
        -o "$SCIGUI_TMPDIR/scigui.tar.gz"

    echo "Extracting SciGUI..."

    tar -xzf \
        "$SCIGUI_TMPDIR/scigui.tar.gz" \
        -C "$SCIGUI_TMPDIR"

    SETUP_PY="$(
        find "$SCIGUI_TMPDIR" \
            -type f \
            -name setup.py \
            -print \
            -quit
    )"

    if [[ -z "$SETUP_PY" ]]; then
        echo "ERROR: Could not locate SciGUI setup.py."
        exit 1
    fi

    SCIGUI_DIR="$(dirname "$SETUP_PY")"

    echo "Installing SciGUI..."

    (
        cd "$SCIGUI_DIR"
        "$PYTHON" -m pip install .
    )

    echo "SciGUI installed successfully."
else
    echo "Skipping SciGUI installation."
fi

###############################################################################
# Final instructions
###############################################################################

echo
echo "Setup complete."
echo

if [[ "$USE_CONDA" -eq 1 ]]; then

    echo "To activate the environment later:"
    echo
    echo "    conda activate $ENV_NAME"
    echo

else

    echo "To activate the environment later:"
    echo
    echo "    source .venv/bin/activate"
    echo

fi

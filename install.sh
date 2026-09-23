#!/usr/bin/env bash

set -euo pipefail

###############################################################################
# Choose environment type
###############################################################################

OS="$(uname -s)"
USE_CONDA=0
USE_VENV=0

if command -v conda >/dev/null 2>&1; then

    echo
    echo "Conda detected."
    echo
    USE_CONDA=1
    #echo "Choose environment type:"
    #echo "  1) Conda (recommended)"
    #echo "  2) venv"
    #echo

    #while true; do
    #    printf "Selection [1/2] (default=1): "
    #    read -r CHOICE
    #    CHOICE="${CHOICE:-1}"
    #    case "$CHOICE" in
    #        1)
    #            USE_CONDA=1
    #            break
    #            ;;
    #        2)
    #            USE_VENV=1
    #            break
    #            ;;
    #        *)
    #            echo "Please enter 1 or 2."
    #            ;;
    #    esac
    #done

else

    echo
    echo "Conda was not found."
    echo

    #printf "Use a Python venv instead? [y/N] "
    #read -r USE_VENV_REPLY
#
#    if [[ "$USE_VENV_REPLY" =~ ^([Yy]|[Yy][Ee][Ss])$ ]]; then
#
#        USE_VENV=1
#
#    else

        echo
	echo "Please install Miniforge (or other conda) and rerun this script."
        echo
        echo "https://conda-forge.org/miniforge/"
        echo

        return 1
    #fi
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
        return 1
    fi

    conda activate "$ENV_NAME"

    if [[ "${CONDA_DEFAULT_ENV:-}" != "$ENV_NAME" ]]; then
        echo "ERROR: Failed to activate Conda environment."
        return 1
    fi

    PYTHON="python"

else

    printf "Enter the name of the venv environment directory: .venv/${ENV_NAME}"
    read -r VENV_DIR
    if [[ -z "$VENV_DIR" ]]
    then
        VENV_DIR="$HOME/.venv/$ENV_NAME"
    fi

    if [[ -e "$VENV_DIR" ]]; then
        echo "ERROR: $VENV_DIR already exists."
        return 1
    fi
    mkdir -p "$(dirname "$VENV_DIR")"
    PYTHON_HOST=""

    for py in python3.13 python3.12; do
        if command -v "$py" >/dev/null 2>&1; then
            PYTHON_HOST="$py"
            break
        fi
    done

    if [[ -z "$PYTHON_HOST" ]]; then
        echo "ERROR: Could not find python3.12 or python3.13."
        return 1
    fi

    echo
    echo "Creating venv..."

    "$PYTHON_HOST" -m venv "$VENV_DIR"

    if [[ ! -x "$VENV_DIR/bin/python" ]]; then
        echo "ERROR: Failed to create venv."
        return 1
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
    corvus="corvus"
    SCIGUI_TMPDIR="$(mktemp -d)"

    cleanup() {
        rm -rf "$SCIGUI_TMPDIR"
    }

    trap cleanup EXIT

    if [[ "$OS" == "Linux" ]]; then
	    echo
	    echo "Installing wxpython."
	    conda install conda-forge::wxpython
    fi
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
        return 1
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
# Create desktop launcher. 
###############################################################################
PROJECT_DIR="${HOME}/corvus_examples"
if [[ "$OS" == "Darwin" ]]; then

    DESKTOP_LAUNCHER="$HOME/Desktop/${ENV_NAME}_GUI.command"
    TERMINAL_LAUNCHER="$HOME/Desktop/${ENV_NAME}.command"

    if [[ "$USE_CONDA" -eq 1 ]]; then

        cat > "$DESKTOP_LAUNCHER" <<EOF
#!/usr/bin/env bash

cd "$PROJECT_DIR"

eval "\$(conda shell.bash hook)"
conda activate "$ENV_NAME"

echo
echo "Activated conda environment: $ENV_NAME"
echo
corvus
exec $SHELL -i
EOF
	cat > "$TERMINAL_LAUNCHER" <<EOF
#!/usr/bin/env bash

cd "$PROJECT_DIR"

eval "\$(conda shell.bash hook)"
conda activate "$ENV_NAME"

echo
echo "Activated conda environment: $ENV_NAME"
echo
exec $SHELL -i
EOF

    else

        cat > "$DESKTOP_LAUNCHER" <<EOF
#!/usr/bin/env bash

cd "$PROJECT_DIR"

source "$VENV_DIR/bin/activate"

echo
echo "Activated virtual environment: $ENV_NAME"
echo
corvus
exec bash -i
EOF

    fi

    chmod +x "$DESKTOP_LAUNCHER"

    echo "Created launcher:"
    echo "  $DESKTOP_LAUNCHER"

elif [[ "$OS" == "Linux" ]]; then

    DESKTOP_GUI_LAUNCHER="$HOME/Desktop/${ENV_NAME}_GUI.desktop"
    DESKTOP_LAUNCHER="$HOME/Desktop/${ENV_NAME}.desktop"

    if [[ "$USE_CONDA" -eq 1 ]]; then

        cat > "$DESKTOP_GUI_LAUNCHER" <<EOF
[Desktop Entry]
Version=1.0
Type=Application
Name=$ENV_NAME
Terminal=true
Exec=bash -c 'cd "$PROJECT_DIR"; eval "\$($CONDA_EXE shell.bash hook)"; conda activate "$ENV_NAME"; corvus exec bash -i'
EOF

        cat > "$DESKTOP_LAUNCHER" <<EOF
[Desktop Entry]
Version=1.0
Type=Application
Name=$ENV_NAME
Terminal=true
Exec=bash -c 'cd "$PROJECT_DIR"; eval "\$($CONDA_EXE shell.bash hook)"; conda activate "$ENV_NAME"; corvus exec bash -i'
EOF

    else

        cat > "$DESKTOP_LAUNCHER" <<EOF
[Desktop Entry]
Version=1.0
Type=Application
Name=$ENV_NAME
Terminal=true
Exec=bash -c 'cd "$PROJECT_DIR"; source "$PROJECT_DIR/$VENV_DIR/bin/activate"; corvus; exec bash -i'
EOF

    fi

    chmod +x "$DESKTOP_LAUNCHER"

    echo "Created launcher:"
    echo "  $DESKTOP_LAUNCHER"

fi
###############################################################################
# Final instructions
###############################################################################
echo "Copying examples to $HOME/corvus_examples"
echo
ex_dir="$HOME/corvus_examples"
ans=y
if [[ -e "$ex_dir" ]]; then
    echo "WARNING: $ex_dir already exists."
    echo "Copy files anyway? [y/N]"
    read ans
fi
if [[ "$ans" =~ ^([Yy]|[Yy][Ee][Ss])$ ]]; then
    cp -r examples/ $HOME/corvus_examples
else
    echo "Will not copy example files."
fi

echo
echo "Setup complete."
echo

if [[ "$USE_CONDA" -eq 1 ]]; then

    echo "You can run corvus by opening a terminal and typing:"
    echo
    echo "    conda activate $ENV_NAME"
    echo
    echo "or double clicking the ${ENV_NAME}.command or"
    echo "${ENV_NAME}_GUI.command scripts located on your"
    echo "desktop."

else

    echo "To activate the environment later:"
    echo
    echo "    source .venv/bin/activate"
    echo

fi

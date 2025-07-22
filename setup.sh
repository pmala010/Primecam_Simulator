#!/bin/bash
set -e

# 1. Install build dependencies
sudo apt update
sudo apt install -y build-essential curl git libssl-dev zlib1g-dev \
  libbz2-dev libreadline-dev libsqlite3-dev wget llvm libncurses5-dev \
  libncursesw5-dev xz-utils tk-dev libffi-dev liblzma-dev

# Setup variables
PYENV_ROOT="$HOME/.pyenv"
PYTHON_VERSION="3.10.4"
VENV_NAME="repo-env"

# 2. Install or update pyenv
if [ ! -d "$PYENV_ROOT" ]; then
  git clone --branch v2.5.5 https://github.com/pyenv/pyenv.git "$PYENV_ROOT"
else
  cd "$PYENV_ROOT"
  git fetch --tags
  git checkout v2.5.5
  cd -
fi

# Add pyenv to path
export PATH="$PYENV_ROOT/bin:$PATH"
eval "$(pyenv init --path)"
eval "$(pyenv init -)"
eval "$(pyenv virtualenv-init -)"

# 3. Install Python
if ! pyenv versions --bare | grep -qx "$PYTHON_VERSION"; then
  pyenv install "$PYTHON_VERSION"
fi

# 4. Create virtualenv
if ! pyenv virtualenvs --bare | grep -qx "$VENV_NAME"; then
  pyenv virtualenv "$PYTHON_VERSION" "$VENV_NAME"
fi

# 5. Activate virtualenv
pyenv activate "$VENV_NAME"

# 6. Set local .python-version
echo "$VENV_NAME" > .python-version

# 7. Install pip and requirements
pip install --upgrade pip==22.0.4
if [ -f requirements.txt ]; then
  pip install -r requirements.txt
else
  echo "No requirements.txt found, skipping pip install."
fi

# 8. Install Jupyter kernel
pip install ipykernel
python -m ipykernel install --user --name="$VENV_NAME" --display-name="Python ($VENV_NAME)"

# 9. Add pyenv setup to ~/.bashrc
if ! grep -q 'pyenv init' ~/.bashrc; then
  cat >> ~/.bashrc << 'EOF'

# Pyenv setup
export PYENV_ROOT="$HOME/.pyenv"
export PATH="$PYENV_ROOT/bin:$PATH"
eval "$(pyenv init --path)"
eval "$(pyenv init -)"
eval "$(pyenv virtualenv-init -)"
EOF
fi

if ! grep -q "pyenv activate $VENV_NAME" ~/.bashrc; then
  echo "pyenv activate $VENV_NAME" >> ~/.bashrc
fi

echo "Setup complete."
echo "Restart your terminal or run: source ~/.bashrc"
echo "Then launch Jupyter and select: Python ($VENV_NAME)"

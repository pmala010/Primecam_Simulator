#!/bin/bash
set -e

# 1. Install dependencies for pyenv and build Python
sudo apt update
sudo apt install -y build-essential curl git libssl-dev zlib1g-dev \
  libbz2-dev libreadline-dev libsqlite3-dev wget llvm libncurses5-dev \
  libncursesw5-dev xz-utils tk-dev libffi-dev liblzma-dev

# 2. Install pyenv 2.5.5 if not installed or wrong version
PYENV_ROOT="$HOME/.pyenv"
if [ ! -d "$PYENV_ROOT" ]; then
  echo "Installing pyenv 2.5.5..."
  git clone --branch v2.5.5 https://github.com/pyenv/pyenv.git "$PYENV_ROOT"
else
  CURRENT_VERSION=$(cd "$PYENV_ROOT" && git describe --tags --abbrev=0)
  if [ "$CURRENT_VERSION" != "v2.5.5" ]; then
    echo "Updating pyenv to version 2.5.5..."
    cd "$PYENV_ROOT"
    git fetch --tags
    git checkout v2.5.5
    cd -
  else
    echo "pyenv 2.5.5 is already installed."
  fi
fi

# 3. Setup pyenv environment for this script
export PYENV_ROOT="$HOME/.pyenv"
export PATH="$PYENV_ROOT/bin:$PATH"
eval "$(pyenv init --path)"
eval "$(pyenv init -)"
eval "$(pyenv virtualenv-init -)"

# 4. Install Python 3.10.4 if missing
if ! pyenv versions --bare | grep -q "^3.10.4\$"; then
  echo "Installing Python 3.10.4..."
  pyenv install 3.10.4
else
  echo "Python 3.10.4 already installed."
fi

# 5. Create virtualenv repo-env if missing
if ! pyenv virtualenvs --bare | grep -q "^repo-env\$"; then
  echo "Creating virtualenv repo-env with Python 3.10.4..."
  pyenv virtualenv 3.10.4 repo-env
else
  echo "Virtualenv repo-env already exists."
fi

# 6. Activate the virtualenv
pyenv activate repo-env

# 7. Upgrade pip to 22.0.4
pip install --upgrade pip==22.0.4

# 8. Install requirements.txt if exists
if [ -f requirements.txt ]; then
  pip install -r requirements.txt
else
  echo "No requirements.txt found, skipping pip install."
fi

# 9. Add pyenv and virtualenv activation to ~/.bashrc if not present
if ! grep -q 'pyenv init' ~/.bashrc; then
  echo "Adding pyenv initialization to ~/.bashrc..."
  cat >> ~/.bashrc << 'EOF'

# Pyenv setup for repo-env
export PYENV_ROOT="$HOME/.pyenv"
export PATH="$PYENV_ROOT/bin:$PATH"
eval "$(pyenv init --path)"
eval "$(pyenv init -)"
eval "$(pyenv virtualenv-init -)"
pyenv activate repo-env

EOF
else
  echo "pyenv initialization already present in ~/.bashrc."
fi

echo "Setup complete! Please restart your terminal or run 'source ~/.bashrc' to start using the environment."

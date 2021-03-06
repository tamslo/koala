#!/bin/bash

# Install node version manager, to use node right after script execution,
# run `source ~/.bashrc`
curl -o- https://raw.githubusercontent.com/creationix/nvm/v0.33.11/install.sh | bash
export NVM_DIR="$HOME/.nvm"
[ -s "$NVM_DIR/nvm.sh" ] && \. "$NVM_DIR/nvm.sh"  # This loads nvm
[ -s "$NVM_DIR/bash_completion" ] && \. "$NVM_DIR/bash_completion"  # This loads nvm bash_completion
nvm install node
nvm use node
nvm alias default node

# Add code folder and virtualenv
mkdir ~/code
virtualenv ~/code/venv -p python3

# Init empty git repository for deployment
git init --bare ~/git

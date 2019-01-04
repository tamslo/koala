#!/bin/bash

sudo apt update
sudo apt install -y git
sudo apt install -y python3
sudo apt install -y virtualenv

# Install node version manager
sudo apt install -y build-essential libssl-dev
curl -o- https://raw.githubusercontent.com/creationix/nvm/v0.33.11/install.sh | bash
export NVM_DIR="$HOME/.nvm"
[ -s "$NVM_DIR/nvm.sh" ] && \. "$NVM_DIR/nvm.sh"  # This loads nvm
[ -s "$NVM_DIR/bash_completion" ] && \. "$NVM_DIR/bash_completion"  # This loads nvm bash_completion
nvm install node
nvm use node
nvm alias default node
source ~/.bashrc

# Install docker
sudo apt install -y \
    apt-transport-https \
    ca-certificates \
    curl \
    software-properties-common
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
sudo add-apt-repository \
   "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
   $(lsb_release -cs) \
   stable"
sudo apt update
sudo apt install -y docker-ce

# Add code folder and virtualenv
sudo mkdir /home/deploy/code
sudo virtualenv /home/deploy/code/venv -p python3

# Init empty git repository for deployment
sudo git init --bare /home/deploy/git

# Make app accessible at default port
sudo iptables -t nat -A PREROUTING -i eth0 -p tcp --dport 80 -j REDIRECT --to-port 5000

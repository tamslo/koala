#!/bin/bash

# Install node version manager manually with deploy user after this script
# with the install_nvm.sh script

sudo apt update
sudo apt install -y git
sudo apt install -y python3
sudo apt install -y virtualenv

# For node version manager
sudo apt install -y build-essential libssl-dev

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

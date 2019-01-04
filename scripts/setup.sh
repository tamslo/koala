#!/bin/bash

# Run setup_deploy.sh script after this one
# Make sure, deploy user exists on the machine (create otherwise)

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
sudo groupadd docker
sudo usermod -aG docker deploy

# Make app accessible at default port
sudo iptables -t nat -A PREROUTING -i eth0 -p tcp --dport 80 -j REDIRECT --to-port 5000

# Install Python requirements
pip install -r requirements.txt --no-cache-dir

# Install Node requirements
cd client
npm install
cd ..

bash scripts/build-docker.sh

# Download data (skip data download for now, comment back in when needed)
# python scripts/download-data.py
# python scripts/download-calls.py

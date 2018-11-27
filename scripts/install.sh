# Install Python requirements
pip install -r requirements.txt --no-cache-dir

# Install Node requirements
cd client
npm install
cd ..

bash scripts/build-docker.sh

# Download data
bash scripts/download-data.sh
python scripts/download-data.py

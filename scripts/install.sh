# Install Python requirements
pip install -r requirements.txt --no-cache-dir

# Install Node requirements
cd client
npm install
cd ..

bash scripts/build-docker.sh

# Download data
python scripts/download-data.py
bash scripts/get-annotations.sh
python scripts/download-calls.py

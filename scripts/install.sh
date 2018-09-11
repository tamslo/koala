# Install Python requirements
pip install -r requirements.txt --no-cache-dir

# Install Node requirements
cd client
npm install
cd ..

# Build services as Docker images
cd services
for dockerfile in ./Dockerfile.*; do
  image_name=${dockerfile#"./Dockerfile."}
  docker build -t $image_name -f $dockerfile .
done
cd ..

# Download data
python scripts/download-data.py

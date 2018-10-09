# Install Python requirements
# pip install -r requirements.txt --no-cache-dir

# Install Node requirements
# cd client
# npm install
# cd ..

# Build services as Docker images
for service in services/*; do
  image_name=${service#"services/"}
  docker build -t $image_name $service
done

# Download data
# python scripts/download-data.py

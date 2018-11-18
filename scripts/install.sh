pip install -r requirements.txt --no-cache-dir
bash scripts/build-docker.sh
python scripts/download-data.py

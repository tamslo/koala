docker-compose up -d
cd client
npm run build
cd ..
echo "Starting server with $(python --version)"
nohup python app.py >> info.log 2>> error.log &

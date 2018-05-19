docker-compose up -d
cd client && npm run build
cd .. && nohub python app.py >> info.log 2>> error.log &

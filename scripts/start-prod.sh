# docker-compose up -d
cd client && npm run build
cd .. && nohub python3 app.py >> info.log 2>> error.log &

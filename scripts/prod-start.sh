docker-compose up -d
cd client
npm run build
cd ..
nohup python app.py >> info.log 2>> error.log &

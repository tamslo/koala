killall python3

cd client
npm run build
cd ..
../scripts/share-constants.sh
echo "Starting server with $(python --version)"
nohup python app.py >> info.log 2>> error.log &

killall python

cd client
npm run build
cd ..
echo "Starting server with $(python --version)"
nohup python app.py >> koala.log 2>&1&

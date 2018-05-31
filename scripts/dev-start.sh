function finish {
  cd ..
}
trap finish EXIT
cd client
npm start

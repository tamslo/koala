# Build services as Docker images
for service in services/*; do
  if [ -d $service ]; then
    image_name=${service#"services/"}
    docker build -t $image_name $service
  fi
done

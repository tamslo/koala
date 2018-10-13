import os
import yaml

services = []
for service in os.listdir("services"):
    service_directory = "services/{}/".format(service)
    if os.path.isdir(service_directory):
        service_yaml_path = service_directory + "service.yml"
        if os.path.exists(service_yaml_path):
            with open(service_yaml_path, "r") as service_yaml:
                services.append(yaml.load(service_yaml))

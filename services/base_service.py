import yaml

class BaseService:
    def __init__(self, config_path):
        with open(config_path, "r") as config_file:
            config = yaml.load(config_file)
            for key in config:
                setattr(self, key, config[key])

    def frontend_information(self):
        return { "id": self.id, "name": self.name, "type": self.type }

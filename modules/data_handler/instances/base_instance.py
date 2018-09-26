import os
import json
import yaml
from time import localtime
import modules.file_utils as file_utils

class BaseInstance:
    def __init__(self, content, path):
        self.content = content
        self.path = path
        with open("constants.yml", "r") as constants_file:
            self.constants = yaml.load(constants_file)
        if not os.path.exists(self.path):
            self.store()

    def setup(self):
        self.content["created"] = localtime()
        self.store()

    def store(self):
        file_utils.write(json.dumps(self.content), self.path)

    def delete(self):
        file_utils.delete(self.path)

    def update(self, content):
        self.content = content
        self.store()

    def get(self, property):
        return self.content[property]

import os
import json
import uuid
from collections import OrderedDict
import modules.file_utils as file_utils

class InstanceHandler:
    def __init__(self, directory, Instance):
        self.directory = directory
        self.Instance = Instance
        file_utils.create_directory(self.directory)

    def __instance_path(self, id):
        return self.directory + id + ".json"

    def create(self, content):
        if not "id" in content:
            content["id"] = str(uuid.uuid4())
        instance = self.Instance(content, self.__instance_path(content["id"]))
        instance.setup()
        return instance

    def update(self, content):
        instance = self.select(content["id"])
        instance.update(content)
        return instance

    def all(self):
        instances = {}
        for path in os.listdir(self.directory):
            if not os.path.isdir(os.path.join(self.directory, path)):
                id = path.split(".json")[0]
                instance = self.select(id)
                instances[instance.get("id")] = instance

        return OrderedDict(
            sorted(
                instances.items(),
                # tuple[0] is experiment_id, tuple[1] is experiment
                key = lambda tuple: tuple[1].get("created")
            )
        )

    def select(self, id):
        path = self.__instance_path(id)
        if os.path.exists(path):
            with open(path, "r") as file:
                content = json.load(file, object_pairs_hook=OrderedDict)
            return self.Instance(content, self.__instance_path(content["id"]))
        else:
            return None

    def delete(self, id):
        instance = self.select(id)
        instance.delete()
        return instance

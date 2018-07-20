import os
import json
import modules.file_utils as file_utils
from .base_handler import BaseHandler
from .base_instance import BaseInstance

class Datasets(BaseHandler):
    def __init__(self, directory):
        super().__init__(directory)
        self.Instance = BaseInstance
        with open("modules/constants/constants.json", "r") as constants_file:
            self.constants = json.load(constants_file)["dataset"]

    def __store_data(self, content, files):
        for file_key in content["data"]:
            file_path = self.directory + content["id"] + "/" + file_key + ".fastq"
            if content["method"] == self.constants["URL"]:
                url = content["data"][file_key]["name"]
                file_utils.download(url, file_path)
            else:
                file = files[file_key]
                name = file.filename
                file.save(file_path)
                content["data"][file_key]["name"] = name
            content["data"][file_key]["path"] = file_path
        return content

    def create(self, content, files=None):
        directory = self.directory + content["id"] + "/"
        os.mkdir(directory)
        content = self.__store_data(content, files)
        return super().create(content)

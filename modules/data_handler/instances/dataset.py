import os
import yaml
import modules.file_utils as file_utils
from .base_instance import BaseInstance

class Dataset(BaseInstance):
    def __init__(self, content, path):
        self.directory = path.split(".json")[0] + "/"
        if "files" in content:
            self.files = content["files"]
            content.pop("files", None)
        super().__init__(content, path)

    def setup(self):
        try:
            with open("constants.yml", "r") as constants_file:
                self.constants = yaml.load(constants_file)["dataset"]
            os.mkdir(self.directory)
            self.__store_data()
            super().setup()
            return super().store()
        except Exception as error:
            file_utils.delete(self.directory)
            file_utils.delete(self.path)
            self.content["error"] = True
            return self

    def __store_data(self):
        for file_key in self.content["data"]:
            file_path = self.directory + file_key + ".fastq"
            if self.content["method"] == self.constants["URL"]:
                url = self.content["data"][file_key]["name"]
                file_utils.download(url, file_path)
            else:
                file = self.files[file_key]
                name = file.filename
                file.save(file_path)
                self.content["data"][file_key]["name"] = name
            self.content["data"][file_key]["path"] = file_path

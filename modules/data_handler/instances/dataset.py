import os
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
        super().setup()
        os.mkdir(self.directory)
        try:
            self.content["error"] = False
            self.__store_data()
        except Exception as error:
            file_utils.delete(self.directory)
            file_utils.delete(self.path)
            self.content["error"] = True
            raise error

    def __store_data(self):
        def is_zipped(name):
            return name.endswith(".gz")

        def maybe_zipped_path(name, path):
            if is_zipped(name):
                return path + ".gz"
            else:
                return path

        def maybe_unzip(path):
            if is_zipped(path):
                file_utils.unzip(path)

        for file_key in self.content["data"]:
            unzipped_file_path = self.directory + file_key + ".fastq"
            if self.content["method"] == self.constants["dataset"]["URL"]:
                url = self.content["data"][file_key]["name"]
                file_path = maybe_zipped_path(url, unzipped_file_path)
                file_utils.download(url, file_path)
                maybe_unzip(file_path)
            else:
                file = self.files[file_key]
                name = file.filename
                file_path = maybe_zipped_path(name, unzipped_file_path)
                file.save(file_path)
                maybe_unzip(file_path)
                self.content["data"][file_key]["name"] = name
            self.content["data"][file_key]["path"] = unzipped_file_path

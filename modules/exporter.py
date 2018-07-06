import os, zipfile, shutil
import modules.file_utils as file_utils

class Exporter:
    def __init__(self, data_directory):
        self.directory = data_directory + "tmp/"
        self.__setup()

    def __setup(self):
        if not os.path.isdir(self.directory):
            os.makedirs(self.directory)

    def file_name(self, path):
        return path.split("/")[-1]

    def zip(self, experiment):
        archive_name = experiment["name"] + ".zip"
        archive_path = self.directory + archive_name
        archive = zipfile.ZipFile(archive_path, "w")
        for key, action in experiment["pipeline"].items():
            if "file" in action:
                path = action["file"]
                file_name = path.split("/")[-1]
                archive.write(path, file_name)
        archive.close()
        return archive_path, archive_name

    def clean_up(self):
        file_utils.delete(self.directory)

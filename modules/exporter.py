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

    def zip(self, directory, archive_name):
        archive_path = self.directory + archive_name
        archive = zipfile.ZipFile(archive_path, "w")
        for root, dirs, files in os.walk(directory):
            for file in files:
                archive.write(
                    os.path.join(root, file),
                    os.path.join(root.replace(directory, ""), file)
                )
        archive.close()
        return archive_path

    def zip_experiment(self, experiment):
        archive_name = experiment["name"] + ".zip"
        archive_path = self.directory + archive_name
        archive = zipfile.ZipFile(archive_path, "w", zipfile.ZIP_DEFLATED)
        for key, action in experiment["pipeline"].items():
            if "file" in action:
                path = action["file"]
                file_name = path.split("/")[-1]
                archive.write(path, file_name)
        archive.close()
        return archive_path, archive_name

    def clean_up(self):
        file_utils.delete(self.directory)

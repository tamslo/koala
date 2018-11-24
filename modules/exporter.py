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
        if path.endswith("/"):
            split_index = -2
        else:
            split_index = -1
        return path.split("/")[split_index]

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

    def zip_experiment_directory(self, path):
        archive_name = self.file_name(path) + ".zip"
        archive_path = self.directory + archive_name
        archive = zipfile.ZipFile(archive_path, "w", zipfile.ZIP_DEFLATED)
        for file_name in os.listdir(path):
            file_path = os.path.join(path, file_name)
            if not os.path.isdir(file_path):
                archive.write(file_path, file_name)
        archive.close()
        return archive_path, archive_name

    def zip_experiment(self, experiment):
        archive_name = experiment.get("name") + ".zip"
        archive_path = self.directory + archive_name
        archive = zipfile.ZipFile(archive_path, "w", zipfile.ZIP_DEFLATED)
        for key, action in experiment.get("pipeline").items():
            if "directory" in action:
                path = action["directory"]
                for file_name in os.listdir(path):
                    file_path = os.path.join(path, file_name)
                    if not os.path.isdir(file_path):
                        archive.write(file_path, os.path.join(action["id"], file_name))
        archive.close()
        return archive_path, archive_name

    def clean_up(self):
        file_utils.delete(self.directory)

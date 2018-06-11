import os, zipfile

class Downloader:
    def __init__(self, download_directory):
        self.directory = download_directory
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
        for key, path in experiment["files"].items():
            file_name = path.split("/")[-1]
            archive.write(path, file_name)
        archive.close()
        return archive_path, archive_name

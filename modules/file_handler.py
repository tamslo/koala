import os, json, urllib

class FileHandler:
    def __init__(self, directory=None):
        self.directory = directory
        self.__setup()

    def __setup(self):
        if self.directory and not os.path.isdir(self.directory):
            os.makedirs(self.directory)

    def write(self, content, path):
        with open(path, "w") as file:
            file.write(content)

    def delete(self, path):
        try:
            os.remove(path)
        except OSError:
            pass

    def download(self, url, destination):
        urllib.request.urlretrieve(
            url,
            destination
        )

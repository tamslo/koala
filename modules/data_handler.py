import os, json

class DataHandler:
    def __init__(self, directory):
        self.directory = directory
        self.__setup()

    def __setup(self):
        if not os.path.isdir(self.directory):
            os.makedirs(self.directory)

    def write(self, content, path):
        with open(path, "w") as file:
            file.write(content)

    def delete(self, path):
        os.remove(path)

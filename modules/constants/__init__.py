import os, json, shutil

class Constants:
    def __init__(self, root_directory):
        constants_path = os.path.join(os.path.dirname(__file__), "constants.json")
        self.__copy_to_frontend(constants_path, root_directory)
        with open(constants_path, "r") as constants_file:
            constants = json.load(constants_file)
        self.actions = constants["actions"]
        self.dataset = constants["dataset"]

    def __copy_to_frontend(self, constants_path, root_directory):
        frontend_path = os.path.join(root_directory, "client/src/constants.json")
        shutil.copyfile(constants_path, frontend_path)

    def actions(self):
        return self.actions

    def dataset(self):
        return self.dataset

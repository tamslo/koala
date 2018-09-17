import os, shutil, time
import modules.file_utils as file_utils

class Cache:
    def __init__(self, directory, error_directory):
        self.directory = directory
        self.error_directory = error_directory
        file_utils.create_directory(self.directory)
        file_utils.create_directory(self.error_directory)

    def __cache_path(self, experiment, action):
        cache_path = (
            self.directory +
            experiment.get("dataset") + "/" +
            experiment.get("reference") + "/"
        )
        action_id = experiment.action_id(action)
        for step_name, step in experiment.get("pipeline").items():
            step_id = step["id"]
            cache_path += (step_id + "/")
            if step_id == action_id:
                break
        return cache_path

    def __error_path(self, experiment, action):
        action_id = experiment.action_id(action)
        timestamp = int(time.time())
        return self.error_directory + action_id + "-" + str(timestamp)

    def lookup(self, experiment, action):
        action_path = self.__cache_path(experiment, action)
        if os.path.isdir(action_path):
            return os.listdir(action_path)[0]

    def create_path(self, experiment, action):
        action_path = self.__cache_path(experiment, action)
        os.makedirs(action_path)
        return action_path

    def clean_up(self, action, experiment):
        action_path = self.__cache_path(action, experiment)
        error_path = self.__error_path(action, experiment)
        shutil.copytree(action_path, error_path)
        file_utils.delete(action_path)

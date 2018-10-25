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
        timestamp = time.strftime("%Y%m%dT%H%M%S", time.localtime())
        return self.error_directory + action_id + "-" + timestamp

    def lookup(self, experiment, action):
        action_path = self.__cache_path(experiment, action)
        if os.path.isdir(action_path):
            return action_path

    def create_path(self, experiment, action):
        action_path = self.__cache_path(experiment, action)
        os.makedirs(action_path)
        return action_path

    def clean_up(self, experiment, action):
        action_path = self.__cache_path(experiment, action)
        error_path = self.__error_path(experiment, action)
        shutil.copytree(action_path, error_path)
        file_utils.delete(action_path)

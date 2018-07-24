import os
import modules.file_utils as file_utils

class Cache:
    def __init__(self, directory):
        self.directory = directory
        file_utils.create_directory(self.directory)

    def __cache_path(self, experiment, action):
        dataset_id = experiment.get("dataset")
        reference_id = experiment.get("reference")
        dataset_directory = self.directory + dataset_id + "/"
        action_id = experiment.action_id(action)
        action_directory = dataset_directory + reference_id + "/" + action_id
        return action_directory

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
        file_utils.delete(action_path)

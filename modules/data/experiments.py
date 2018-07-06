import os, json
from time import localtime, strftime
from collections import OrderedDict
import modules.file_utils as file_utils

class Experiments:
    def __init__(self, data_directory):
        self.directory = data_directory + "experiments/"
        file_utils.create_directory(self.directory)

    def __write(self, experiment):
        content = json.dumps(experiment)
        path = self.experiment_path(experiment["id"])
        file_utils.write(content, path)
        return experiment

    def __set(self, experiment, name, value):
        experiment[name] = value
        return self.__write(experiment)

    def log_action_status(self, experiment, action, status):
        time = localtime()
        experiment["pipeline"][action][status] = time
        return self.__write(experiment)

    def create(self, experiment):
        time = localtime()
        experiment["created"] = time
        experiment["error"] = False
        experiment["done"] = False
        experiment["interrupted"] = False
        experiment["report"] = None
        return self.__write(experiment)

    def select(self, experiment_id):
        with open(self.experiment_path(experiment_id), "r") as file:
            return json.load(file, object_pairs_hook=OrderedDict)

    def delete(self, experiment_id):
        experiment = self.select(experiment_id)
        file_utils.delete(self.experiment_path(experiment_id))
        return experiment

    def update(self, experiment):
        return self.__write(experiment)

    def start_action(self, experiment, action, cached = False):
        experiment = self.log_action_status(experiment, action, "started")
        if cached:
            experiment["pipeline"][action]["cached"] = True
            experiment = self.complete_action(experiment, action)
        return self.__write(experiment)

    def add_download(self, experiment, action, path):
        experiment["pipeline"][action]["file"] = path
        return self.__write(experiment)

    def complete_action(self, experiment, action):
        return self.log_action_status(experiment, action, "completed")

    def mark_error(self, experiment_id, action, error):
        experiment = self.select(experiment_id)
        experiment = self.log_action_status(experiment, action, "error")
        return self.__set(experiment, "error", str(error))

    def mark_interrupted(self, experiment_id, action):
        experiment = self.select(experiment_id)
        experiment = self.log_action_status(experiment, action, "interrupted")
        return self.__set(experiment, "interrupted", True)

    def mark_done(self, experiment_id):
        time = localtime()
        experiment = self.select(experiment_id)
        return self.__set(experiment, "done", time)

    def all(self):
        experiments = {}
        for path in os.listdir(self.directory):
            experiment = self.select(path.split(".json")[0])
            experiments[experiment["id"]] = experiment
        return OrderedDict(
            sorted(
                experiments.items(),
                # tuple[0] is experiment_id, tuple[1] is experiment
                key = lambda tuple: tuple[1]["created"]
            )
        )
        return experiments

    def experiment_path(self, experiment_id):
        return self.directory + experiment_id + ".json"

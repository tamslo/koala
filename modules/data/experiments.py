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

    def __last_log_entry(self, experiment, action):
        return [entry for entry in experiment["log"] if entry["action"] == action][-1]

    def __edit_log_entry(self, experiment, action, field):
        time = localtime()
        entry = self.__last_log_entry(experiment, action)
        entry[field] = time
        return self.__write(experiment)

    def create(self, experiment):
        experiment["error"] = False
        experiment["done"] = False
        experiment["interrupted"] = False
        experiment["report"] = None
        experiment["log"] = []
        experiment["files"] = {}
        return self.add_log_entry(experiment, "create", one_step = True)

    def select(self, experiment_id):
        with open(self.experiment_path(experiment_id), "r") as file:
            return json.load(file, object_pairs_hook=OrderedDict)

    def delete(self, experiment_id):
        experiment = self.select(experiment_id)
        file_utils.delete(self.experiment_path(experiment_id))
        return experiment

    def update(self, experiment):
        experiment = self.__write(experiment)
        return self.add_log_entry(experiment, "update", one_step = True)

    def add_log_entry(self, experiment, action, one_step = False):
        time = localtime()
        entry = ({
            "action": action,
            "started": time,
            "completed": False,
            "error": False
        })
        if one_step:
            entry["completed"] = time
        experiment["log"].append(entry)
        return self.__write(experiment)

    def add_download(self, experiment, key, path):
        experiment["files"][key] = path
        return self.__write(experiment)

    def log_complete(self, experiment, action):
        return self.__edit_log_entry(experiment, action, "completed")

    def log_error(self, experiment):
        action = experiment["log"][-1]["action"]
        return self.__edit_log_entry(experiment, action, "error")

    def mark_error(self, experiment_id, error):
        experiment = self.select(experiment_id)
        experiment = self.log_error(experiment)
        return self.__set(experiment, "error", str(error))

    def mark_interrupted(self, experiment_id):
        experiment = self.select(experiment_id)
        return self.__set(experiment, "interrupted", True)

    def mark_done(self, experiment_id):
        experiment = self.select(experiment_id)
        experiment = self.add_log_entry(experiment, "done", one_step = True)
        return self.__set(experiment, "done", True)

    def all(self):
        experiments = {}
        for path in os.listdir(self.directory):
            experiment = self.select(path.split(".json")[0])
            experiments[experiment["id"]] = experiment
        return OrderedDict(
            sorted(
                experiments.items(),
                key = lambda tuple: self.__last_log_entry(tuple[1], "create")["started"]
            )
        )
        return experiments

    def experiment_path(self, experiment_id):
        return self.directory + experiment_id + ".json"

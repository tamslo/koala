import os, json, uuid
from time import localtime, strftime
from collections import OrderedDict

class Experiments:
    def __init__(self, data_directory):
        self.directory = data_directory + "experiments/"
        self.__setup()

    def __setup(self):
        if not os.path.isdir(self.directory):
            os.makedirs(self.directory)

    def __write(self, experiment):
        with open(self.experiment_path(experiment["id"]), "w") as file:
            file.write(json.dumps(experiment))
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
        experiment["id"] = str(uuid.uuid4())
        experiment["error"] = False
        experiment["done"] = False
        experiment["interrupted"] = False
        experiment["report"] = None
        experiment["log"] = []
        return self.add_log_entry(experiment, "create", one_step = True)

    def select(self, experiment_id):
        with open(self.experiment_path(experiment_id), "r") as file:
            return json.load(file, object_pairs_hook=OrderedDict)

    def delete(self, experiment_id):
        experiment = self.select(experiment_id)
        os.remove(self.experiment_path(experiment_id))
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

    def log_complete(self, experiment, name):
        return self.__edit_log_entry(experiment, name, "completed")

    def log_error(self, experiment):
        name = list(experiment["log"].keys())[-1]
        return self.__edit_log_entry(experiment, name, "error")

    def mark_error(self, experiment_id, error):
        experiment = self.select(experiment_id)
        experiment = self.log_error(experiment)
        return self.__set(experiment, "error", str(error))

    def mark_interrupted(self, experiment_id):
        experiment = self.select(experiment_id)
        experiment = self.add_log_entry(experiment, "interrupted", one_step = True)
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

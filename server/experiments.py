import os, json, uuid, sys
from time import gmtime, strftime

class Experiments:
    def __init__(self, data_directory):
        self.directory = data_directory + "experiments/"
        self.__setup()

    def __setup(self):
        if not os.path.isdir(self.directory):
            os.makedirs(self.directory)

    def __status(self, name):
        return {
            "time": strftime("%a, %d %b %Y %H:%M", gmtime()),
            "name": name
        }

    def create(self, experiment):
        experiment_id = str(uuid.uuid4())
        experiment["id"] = experiment_id
        experiment["error"] = False
        experiment["done"] = False
        experiment["report"] = None
        experiment["statuses"] = [self.__status("Created experiment")]
        experiment = json.dumps(experiment)

        with open(self.experiment_path(experiment_id), "w") as file:
            file.write(experiment)
            return experiment

    def select(self, experiment_id):
        with open(self.experiment_path(experiment_id), "r") as file:
            return json.load(file)

    def set(self, experiment_id, name, value):
        experiment = self.select(experiment_id)
        experiment[name] = value
        experiment = json.dumps(experiment)
        with open(self.experiment_path(experiment_id), "w") as file:
            file.write(experiment)
            return experiment

    def add_status(self, experiment_id, status):
        experiment = self.select(experiment_id)
        statuses = experiment["statuses"]
        statuses.append(self.__status(status))
        return self.set(experiment_id, "statuses", statuses)

    def mark_error(self, experiment_id, error):
        experiment = self.select(experiment_id)
        experiment = self.add_status(experiment_id, "Error")
        return self.set(experiment_id, "error", str(error))

    def mark_done(self, experiment_id):
        experiment = self.select(experiment_id)
        experiment = self.add_status(experiment_id, "Done")
        return self.set(experiment_id, "done", True)

    def all(self):
        experiments = {}
        for path in os.listdir(self.directory):
            with open(self.directory + path, "r") as file:
                experiment = json.load(file)
            experiment_id = path.split("json")[0]
            experiments[experiment_id] = experiment
        return experiments

    def experiment_path(self, experiment_id):
        return self.directory + experiment_id + ".json"

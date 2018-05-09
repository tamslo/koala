import os, json, uuid, sys

class Experiments:
    def __init__(self, data_directory):
        self.directory = data_directory + "experiments/"
        self.__setup()

    def __setup(self):
        if not os.path.isdir(self.directory):
            os.makedirs(self.directory)

    def create(self, params):
        experiment_id = str(uuid.uuid4())
        params["id"] = experiment_id
        experiment = json.dumps(params)
        with open(self.directory + experiment_id + ".json", "w") as file:
            file.write(experiment)
            return experiment

    def select(self, experiment_id):
        with open(self.directory + experiment_id + ".json", "r") as file:
            return json.load(file)

    def all(self):
        experiments = {}
        for path in os.listdir(self.directory):
            with open(self.directory + path, "r") as file:
                experiment = json.load(file)
            experiment_id = path.split("json")[0]
            experiments[experiment_id] = experiment
        return experiments

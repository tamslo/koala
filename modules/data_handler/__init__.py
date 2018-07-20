from .experiments import Experiments
from .datasets import Datasets
from .cache import Cache

class DataHandler:
    def __init__(self, data_directory):
        self.experiments_directory = data_directory + "experiments/"
        self.datasets_directory = data_directory + "datasets/"

        self.experiments = Experiments(self.experiments_directory)
        self.datasets = Datasets(self.datasets_directory)
        self.cache = Cache(self.datasets_directory)

    def clean_up(self):
        # In case of an interruption, clean up experiments and datasets.
        # If an error occurred, the components already handled it.
        for experiment_id, experiment in self.experiments.all().items():
            if not experiment.get("done") and not experiment.get("error"):
                for action, pipeline_step in experiment.get("pipeline").items():
                    started = "started" in pipeline_step and pipeline_step["started"]
                    completed = "completed" in pipeline_step and pipeline_step["completed"]
                    if started and not completed:
                        experiment.mark_interrupted(action)
                        self.experiments.store(experiment)
                        self.cache.clean_up(action, experiment)

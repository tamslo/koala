import os
import modules.file_utils as file_utils
from .instance_handler import InstanceHandler
from .instances.experiment import Experiment
from .instances.dataset import Dataset
from .cache import Cache

class DataHandler:
    def __init__(self, data_directory):
        self.experiments_directory = data_directory + "experiments/"
        self.datasets_directory = data_directory + "datasets/"
        self.error_directory = data_directory + "errored/"
        self.reference_directory = data_directory + "references/"

        self.experiments = InstanceHandler(self.experiments_directory, Experiment)
        self.datasets = InstanceHandler(self.datasets_directory, Dataset)
        self.cache = Cache(self.datasets_directory, self.error_directory)

    def reference_path(self, experiment):
        reference_id = experiment.get("reference")
        return self.reference_directory + "{}.fa".format(reference_id)

    def genome_index_path(self, experiment, aligner):
        reference_id = experiment.get("reference")
        return self.reference_directory + "{}_{}_index".format(reference_id, aligner)

    def clean_up(self):
        # In case of an interruption, clean up experiments, references, and
        # datasets.
        # If an error occurred, the components already handled it.
        for reference in os.listdir(self.reference_directory):
            if reference.endswith(".running"):
                file_utils.delete(os.path.join(self.reference_directory, reference))
        for experiment_id, experiment in self.experiments.all().items():
            if not experiment.get("done") and not experiment.get("error"):
                for action, pipeline_step in experiment.get("pipeline").items():
                    started = "started" in pipeline_step and pipeline_step["started"]
                    completed = "completed" in pipeline_step and pipeline_step["completed"]
                    if started and not completed:
                        experiment.mark_interrupted(action)
                        self.cache.clean_up(experiment, action)

import docker, yaml, traceback
from modules.file_handler import FileHandler

class Runner:
    def __init__(self, datasets, experiments, data_directory, constants):
        self.datasets = datasets
        self.experiments = experiments
        self.action_names = constants.actions
        self.actions = {
            self.action_names["DATASET"]: self.__get_dataset,
            self.action_names["ALIGNMENT"]: self.__align
        }
        self.file_handler = FileHandler()
        self.docker_client = docker.from_env()
        with open("config.yml", "r") as config_file:
            config = yaml.load(config_file)
            self.absolute_data_path = config["absolute_repo_path"] + data_directory

    def execute(self, action, experiment_id):
        try:
            experiment = self.experiments.select(experiment_id)
            file_path = self.datasets.lookup(experiment, action)
            if file_path:
                experiment = self.experiments.add_log_entry(
                    experiment,
                    "{}_loaded".format(action),
                    one_step = True
                )
            else:
                experiment = self.experiments.add_log_entry(experiment, action)
                file_path = self.actions[action](experiment)
                experiment = self.experiments.log_complete(experiment, action)
            experiment = self.experiments.add_download(experiment, action, file_path)
        except Exception as error:
            print(error, flush=True)
            traceback.print_exc()
            experiment = self.experiments.mark_error(experiment_id, error)
            self.datasets.clean_up(action, experiment)
        return experiment

    def __get_dataset(self, experiment):
        dataset_id = experiment[self.action_names["DATASET"]]
        dataset = self.datasets.select(dataset_id)
        dataset_folder = self.datasets.dataset_folder(dataset_id)
        for file in dataset["content"]:
            if dataset["method"] == self.constants["dataset"]["URL"]:
                destination = dataset["content"][file]["path"]
                url = dataset["content"][file]["name"]
                file_handler.download(url, destination)
        return dataset_folder

    def __align(self, experiment):
        file_ending = "bam"
        alignment_path = self.datasets.create_path(
            experiment,
            self.action_names["ALIGNMENT"]
        )
        aligner = experiment[self.action_names["ALIGNMENT"]]
        self.docker_client.containers.run(
            aligner,
            "touch {}/dummy.bam".format(alignment_path),
            volumes={
                self.absolute_data_path: {
                    "bind": "/data",
                    "mode": "rw"
                }
            },
            auto_remove=True
        )
        return alignment_path

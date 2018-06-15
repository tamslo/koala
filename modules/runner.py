import urllib, docker, yaml

class Runner:
    def __init__(self, datasets, experiments, data_directory):
        self.datasets = datasets
        self.experiments = experiments
        self.actions = {
            "dataset": self.__get_dataset,
            "alignment": self.__align
        }
        self.docker_client = docker.from_env()
        with open("config.yml", "r") as config_file:
            config = yaml.load(config_file)
            self.absolute_data_path = config["absolute_repo_path"] + data_directory

    def execute(self, action, experiment_id):
        try:
            experiment = self.experiments.select(experiment_id)
            file_path = self.datasets.lookup(experiment, action)
            print("Executing action " + action, flush=True)
            print(file_path, flush=True)
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
            experiment = self.experiments.mark_error(experiment_id, error)
            self.datasets.clean_up(action, experiment)
        return experiment

    def __get_dataset(self, experiment):
        dataset = self.datasets.select(experiment["dataset"])
        dataset_folder = self.datasets.dataset_folder(dataset["id"])
        for file in dataset["content"]:
            path = dataset["content"][file]["path"]
            url = dataset["content"][file]["origin"]
            urllib.request.urlretrieve(
                url,
                path
            )
        return dataset_folder

    def __align(self, experiment):
        file_ending = "bam"
        alignment_path = self.datasets.create_path(experiment, "alignment")
        self.docker_client.containers.run(
            "star",
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

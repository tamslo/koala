import urllib, docker, yaml

class Runner:
    def __init__(self, cache, experiments, data_directory):
        self.cache = cache
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
        experiment = self.experiments.select(experiment_id)
        experiment_data = self.cache.get_experiment_data(experiment)
        try:
            if experiment_data[action]:
                file_path = experiment_data[action]
                experiment = self.experiments.add_log_entry(
                    experiment,
                    "{}_cached".format(action),
                    one_step = True
                )
            else:
                experiment = self.experiments.add_log_entry(experiment, action)
                file_path = self.actions[action](experiment)
                experiment = self.experiments.log_complete(experiment, action)
            experiment = self.experiments.add_download(experiment, action, file_path)
        except Exception as error:
            self.cache.clean_up(action, experiment)
            experiment = self.experiments.mark_error(experiment_id, error)
        return experiment

    def __get_dataset(self, experiment):
        dataset_path, headers = urllib.request.urlretrieve(
            experiment["dataset"],
            self.cache.create_dataset(experiment["dataset"])
        )
        return dataset_path

    def __align(self, experiment):
        alignment_directory = self.cache.create_path(experiment, "alignment")
        self.docker_client.containers.run(
            "star",
            "touch {}/alignment.bam; ".format(alignment_directory) +
            "echo 'TEST' > {}/alignment.bam".format(alignment_directory),
            volumes={
                self.absolute_data_path: {
                    "bind": "/data",
                    "mode": "rw"
                }
            },
            auto_remove=True
        )
        return alignment_directory + "/" + "alignment.bam"

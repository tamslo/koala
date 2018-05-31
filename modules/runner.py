import urllib

class Runner:
    def __init__(self, cache, experiments):
        self.cache = cache
        self.experiments = experiments
        self.step_actions = {
            "dataset": self.__get_dataset,
            "alignment": self.__align
        }

    def execute(self, step, experiment_id):
        experiment = self.experiments.select(experiment_id)
        experiment_data = self.cache.get_experiment_data(experiment)
        try:
            if experiment_data[step]:
                file_path = experiment_data[step]
                experiment = self.experiments.add_log_entry(
                    experiment,
                    "{}_cached".format(step),
                    one_step = True
                )
            else:
                print("Creating start log", flush=True)
                experiment = self.experiments.add_log_entry(experiment, step)
                print("Running action", flush=True)
                file_path = self.step_actions[step](experiment)
                print("Creating end log", flush=True)
                experiment = self.experiments.log_complete(experiment, step)
            experiment = self.experiments.add_download(experiment, step, file_path)
        except Exception as error:
            print("Error in step " + step + ": " + str(error))
            self.cache.clean_up(step, experiment)
            experiment = self.experiments.mark_error(id, error)
        return experiment

    def __get_dataset(self, experiment):
        print("Getting dataset", flush=True)
        try:
            dataset_path, headers = urllib.request.urlretrieve(
                experiment["dataset"],
                self.cache.create_dataset(experiment["dataset"])
            )
        except Exception as error:
            print(error, flush=True)
        return dataset_path

    def __align(self, experiment):
        return ""

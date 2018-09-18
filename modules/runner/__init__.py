import traceback
import modules.file_utils as file_utils
from .docker import Docker
from .alignment import align

class Runner:
    def __init__(self, data_handler, data_directory, constants):
        self.data_handler = data_handler
        self.action_names = constants.actions
        self.actions = {
            self.action_names["ALIGNMENT"]: align
        }
        self.docker_client = Docker(data_directory)

    def execute(self, experiment_id):
        current_action = ""
        experiment = self.data_handler.experiments.select(experiment_id)
        try:
            for action in experiment.get("pipeline"):
                current_action = action
                experiment = self.__execute_step(action, experiment)
            experiment.mark_done(experiment_id)
        except Exception as error:
            print("[Error in {}] {}".format(current_action, error), flush=True)
            traceback.print_exc()
            experiment.mark_error(current_action, error)
            self.data_handler.cache.clean_up(experiment, current_action)
        return experiment

    def __execute_step(self, action, experiment):
        file_path = self.data_handler.cache.lookup(experiment, action)
        if file_path:
            experiment.start_action(
                action,
                cached = True
            )
        else:
            experiment.start_action(action)
            file_path = self.actions[action](
                self.docker_client,
                self.data_handler,
                experiment,
                self.action_names
            )
            experiment.complete_action(action)
        experiment.add_download(action, file_path)
        return experiment

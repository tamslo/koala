import os
import re
import time
import traceback
import modules.file_utils as file_utils
from .docker import Docker
from services import get_services

class Runner:
    def __init__(self, data_handler, data_directory, constants):
        self.data_handler = data_handler
        self.experiment_statuses = constants["experiment"]
        self.docker_client = Docker(data_directory)
        self.tasks = []
        self.current_task = None
        self.services = get_services()

    def add_task(self, experiment_id):
        experiment = self.data_handler.experiments.select(experiment_id)
        hampering_statuses = [
            self.experiment_statuses["DONE"],
            self.experiment_statuses["RUNNING"]
        ]
        if (
            not experiment.get("status") in hampering_statuses and
            not experiment.get("id") in self.tasks
        ):
            self.tasks.append(experiment_id)

    def remove_task(self, experiment_id):
        experiment_id in self.tasks and self.tasks.remove(experiment_id)

    def run(self):
        if len(self.tasks) > 0:
            self.__execute_next_task()

    def __execute_next_task(self):
        self.current_task = self.tasks.pop(0)
        current_action = "" # For error handling
        experiment = self.data_handler.experiments.select(self.current_task)
        if experiment != None:
            try:
                experiment.mark_running()
                for action in experiment.get("pipeline"):
                    current_action = action
                    experiment = self.__execute_step(action, experiment)
                experiment.mark_done()
            except Exception as error:
                print("[Error in {}] {}".format(current_action, error), flush=True)
                traceback.print_exc()
                experiment.mark_error(current_action, error)
                self.data_handler.cache.clean_up(experiment, current_action)
            return experiment

    def __execute_step(self, action, experiment):
        action_handler_id = experiment.get("pipeline")[action]["id"]
        action_handler = next(
            (service for service in self.services if service.id == action_handler_id),
            None
        )
        file_path = self.data_handler.cache.lookup(experiment, action)
        if file_path: # Step is cached
            experiment.start_action(
                action,
                cached = True
            )
        else: # Run step
            experiment.start_action(action)
            file_path = self.data_handler.cache.create_path(
                experiment,
                action
            )
            file_utils.create_directory(file_path)
            action_handler.run({
                "docker_client": self.docker_client,
                "data_handler": self.data_handler,
                "experiment": experiment,
                "destination": file_path
            })
            self.__run_evaluation_if_specified(file_path, experiment, action)
            experiment.complete_action(action)
        experiment.add_download(action, file_path)
        return experiment

    def __run_evaluation_if_specified(self, destination, experiment, step):
        dataset = self.data_handler.datasets.select(experiment.get("dataset"))
        evaluation = dataset.get("evaluation")
        if evaluation:
            evaluation_type = evaluation["type"]
            evaluation_handler = next(
                (service for service in self.services if service.id == evaluation_type),
                None
            )
            step_id = re.sub(r"_\d$", "", step) # in case a number was appended
            if (step_id in evaluation_handler.steps):
                evaluation_handler.run({
                    "docker_client": self.docker_client,
                    "destination": destination,
                    "dataset": dataset
                })

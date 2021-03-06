import os
import re
import time
import traceback
import logging
from apscheduler.schedulers.background import BackgroundScheduler
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
        self.scheduler = BackgroundScheduler()
        self.scheduler.add_job(
            func=self.check_run,
            trigger="interval",
            seconds=5,
            timezone="Europe/Berlin"
        )
        self.scheduler.start()
        logging.getLogger('apscheduler').setLevel("ERROR")

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

    def check_run(self):
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
            # self.current_task = None
            return experiment

    def __execute_step(self, action, experiment):
        action_handler_id = experiment.get("pipeline")[action]["id"]
        action_handler = next(
            (service for service in self.services if service.id == action_handler_id),
            None
        )
        file_path = self.data_handler.cache.lookup(experiment, action)
        if file_path: # Step is cached
            # Check if evaluation was run, maybe run evaluation
            evaluation_handler = self.__evaluation_handler(experiment, action)
            has_evaluation = next(
                (file for file in os.listdir(file_path) if file.startswith("Evaluation")),
                False
            )
            if evaluation_handler and not has_evaluation:
                experiment.start_action(action)
                experiment.mark_cached(action)
                self.__run_evaluation(file_path, experiment, evaluation_handler, action_handler)
                experiment.complete_action(action)
            else:
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
            self.__run_evaluation_if_specified(file_path, experiment, action, action_handler)
            experiment.complete_action(action)
        experiment.add_download(action, file_path)
        return experiment

    def __evaluation_handler(self, experiment, step):
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
                return evaluation_handler

    def __run_evaluation(self, destination, experiment, evaluation_handler, action_handler):
        dataset = self.data_handler.datasets.select(experiment.get("dataset"))
        evaluation_handler.run({
            "docker_client": self.docker_client,
            "destination": destination,
            "experiment": experiment,
            "dataset": dataset,
            "action_handler": action_handler
        })


    def __run_evaluation_if_specified(self, destination, experiment, step, action_handler):
        evaluation_handler = self.__evaluation_handler(experiment, step)
        if evaluation_handler:
            self.__run_evaluation(destination, experiment, evaluation_handler, action_handler)

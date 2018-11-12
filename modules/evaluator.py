#!/usr/bin/env python3

import os
import shutil
import json
import time

class Evaluator:
    def __init__(self, data_directory):
        self.directory = data_directory + "evaluation/"
        self.evaluation_accumulators = {
            "beers": self.accumulate_beers,
            "giab": self.accumulate_giab
        }
        self.__setup()

    def __setup(self):
        if not os.path.isdir(self.directory):
            os.makedirs(self.directory)

    def assemble_step_id(self, current_step_id, pipeline):
        step_id = ""
        for step_name, step_content in pipeline.items():
            this_step_id = step_content["id"]
            step_id += "_" + this_step_id
            if this_step_id == current_step_id:
                break
        return step_id

    def collect_results(self, data_handler):
        current_direcory = self.directory + time.strftime(
            "%Y%m%dT%H%M%S", time.localtime()
        ) + "/"
        runstats_directory = current_direcory + "runstats/"
        # Extract information from evaluation files
        for experiment_id, experiment in data_handler.experiments.all().items():
            if experiment.get("status") != "done":
                continue

            dataset = data_handler.datasets.select(experiment.get("dataset"))
            evaluation_type = dataset.get("evaluation") and dataset.get("evaluation")["type"]
            evaluation_directory = evaluation_type and current_direcory + evaluation_type + "/"

            pipeline = experiment.get("pipeline")
            for step_name, step_content in pipeline.items():
                directory = step_content["directory"]
                file_name = "{}_{}_{}.txt".format(
                    dataset.get("id"),
                    experiment.get("reference"),
                    self.assemble_step_id(step_content["id"], pipeline)
                )

                runstats_path = runstats_directory + file_name
                if not os.path.isdir(runstats_directory):
                    os.makedirs(runstats_directory)

                if not os.path.exists(runstats_path):
                    if "Runstats.txt" in os.listdir(directory):
                        shutil.copy(
                            directory + "Runstats.txt",
                            runstats_path
                        )

                if not evaluation_type:
                    continue

                if not os.path.isdir(evaluation_directory):
                    os.makedirs(evaluation_directory)

                evaluation_path = evaluation_directory + file_name
                if not os.path.exists(evaluation_path):
                    if "Evaluation.txt" in os.listdir(directory):
                        shutil.copy(
                            directory + "Evaluation.txt",
                            evaluation_path
                        )
        self.accumulate_runstats(runstats_directory)
        if evaluation_type in self.evaluation_accumulators:
            self.evaluation_accumulators[evaluation_type](evaluation_directory)
        return current_direcory

    def accumulate_runstats(self, directory):
        return None

    def accumulate_beers(self, directory):
        return None

    def accumulate_giab(self, directory):
        return None

#!/usr/bin/env python3

import os
import shutil
import json
import time

def accumulate_runstats(directory):
    return None

def accumulate_beers(directory):
    return None

def accumulate_giab(directory):
    return None

def assemble_step_id( current_step_id, pipeline):
    step_id = ""
    for step_name, step_content in pipeline.items():
        this_step_id = step_content["id"]
        step_id += "_" + this_step_id
        if this_step_id == current_step_id:
            break
    return step_id

def collect_evaluation_results(data_directory, data_handler):
    evaluation_directory = data_directory + "evaluation/"
    if not os.path.isdir(evaluation_directory):
        os.makedirs(evaluation_directory)

    evaluation_accumulators = {
        "beers": accumulate_beers,
        "giab": accumulate_giab
    }

    current_direcory = evaluation_directory + time.strftime(
        "%Y%m%dT%H%M%S", time.localtime()
    ) + "/"
    if not os.path.isdir(current_direcory):
        os.makedirs(current_direcory)
    runstats_directory = current_direcory + "runstats/"
    evaluation_types = {}

    # Extract information from evaluation files
    for experiment_id, experiment in data_handler.experiments.all().items():
        if experiment.get("status") != "done":
            continue

        dataset = data_handler.datasets.select(experiment.get("dataset"))
        evaluation_type = dataset.get("evaluation") and dataset.get("evaluation")["type"]
        evaluation_directory = evaluation_type and current_direcory + evaluation_type + "/"
        if evaluation_type and evaluation_type not in evaluation_types:
            evaluation_types[evaluation_type] = evaluation_directory

        pipeline = experiment.get("pipeline")
        for step_name, step_content in pipeline.items():
            directory = step_content["directory"]
            file_name = "{}_{}{}".format(
                dataset.get("id"),
                experiment.get("reference"),
                assemble_step_id(step_content["id"], pipeline)
            )

            runstats_path = runstats_directory + file_name + ".txt"
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

            evaluation_results = []
            for file in os.listdir(directory):
                if file.startswith("Evaluation"):
                    evaluation_results.append(file)

            if len(evaluation_results) == 0:
                continue
            elif len(evaluation_results) == 1:
                evaluation_path = evaluation_directory + file_name + ".txt"
                if not os.path.exists(evaluation_path):
                    shutil.copy(
                        directory + evaluation_results[0],
                        evaluation_path
                    )
            else:
                evaluation_path = evaluation_directory + file_name + "/"
                if not os.path.exists(evaluation_path):
                    os.makedirs(evaluation_path)
                    for file in evaluation_results:
                        shutil.copy(
                            directory + file,
                            evaluation_path + file
                        )

    accumulate_runstats(runstats_directory)
    for evaluation_type, evaluation_directory in evaluation_types.items():
        if evaluation_type in evaluation_accumulators:
            evaluation_accumulators[evaluation_type](evaluation_directory)
    return current_direcory

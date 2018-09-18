from .star import star
from .novoalign import novoalign
from .test import test

aligner_actions = {
    "test": test,
    "star": star,
    "novoalign": novoalign
}

def align(docker_client, data_handler, experiment, action_names):
    alignment_path = data_handler.cache.create_path(
        experiment,
        action_names["ALIGNMENT"]
    )
    aligner = experiment.get("pipeline")[action_names["ALIGNMENT"]]["id"]
    return aligner_actions[aligner](docker_client, alignment_path, data_handler, experiment)

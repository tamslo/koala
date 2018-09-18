from .star import star
from .novoalign import novoalign
from .test import test

aligner_actions = {
    "test": test,
    "star": star,
    "novoalign": novoalign,
}

def align(docker_client, aligner, destination, reference_id, dataset):
    return aligner_actions[aligner](docker_client, destination, reference_id, dataset)

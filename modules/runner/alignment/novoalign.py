import yaml, os
import modules.file_utils as file_utils

def novoalign(docker_client, destination, data_handler, experiment):
    file_utils.create_directory(destination)
    genome_index_path = data_handler.genome_index_path(experiment, "novoalign")

    try:
        build_genome_index(
            docker_client,
            genome_index_path,
            data_handler,
            experiment
        )
    except:
        file_utils.delete(genome_index_path)
        raise

    dataset = data_handler.datasets.select(experiment.get("dataset"))
    run(docker_client, dataset, genome_index_path, destination)
    return destination;

def build_genome_index(docker_client, genome_index_path, data_handler, experiment):
    if not os.path.exists(genome_index_path):
        reference_id = experiment.get("reference")
        reference_path = data_handler.reference_path(experiment)
        command = "novoindex -n {} {} {}".format(reference_id, genome_index_path, reference_path)
        docker_client.run(
            "novoalign",
            command
        )

def run(docker_client, dataset, genome_index_path, destination):
    command = "novoalign -f "
    for direction, specification in dataset.get("data").items():
        command += " /{}".format(specification["path"])
    command += " -d {}".format(genome_index_path)
    docker_client.run(
        "novoalign",
        command
    )

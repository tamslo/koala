import yaml, os
import modules.file_utils as file_utils

def novoalign(docker_client, destination, data_handler, experiment):
    file_utils.create_directory(destination)
    reference_id = experiment.get("reference")
    dataset = data_handler.datasets.select(experiment.get("dataset"))
    reference_directory = "data/references"
    reference_path = "{}/{}.fa".format(reference_directory, reference_id)
    genome_index_path = "{}/{}_novo_index".format(reference_directory, reference_id)

    try:
        build_genome_index(
            docker_client,
            genome_index_path,
            reference_id,
            reference_path
        )
    except:
        file_utils.delete(genome_index_path)
        raise

    run(docker_client, dataset, genome_index_path, destination)
    return destination;

def build_genome_index(docker_client, genome_index_path, reference_id, reference_path):
    if not os.path.exists(genome_index_path):
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

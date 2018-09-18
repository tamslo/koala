import yaml, os
import modules.file_utils as file_utils

def star(docker_client, destination, data_handler, experiment):
    with open("config.yml", "r") as config_file:
        config = yaml.load(config_file)
        num_threads = int(config["cores"])

    genome_index_path = data_handler.genome_index_path(experiment, "star")
    star_preamble = "STAR --runThreadN {} --genomeDir /{} --outFileNamePrefix {}".format(
        num_threads,
        "/" + genome_index_path,
        destination
    )

    try:
        reference_path = data_handler.reference_path(experiment)
        build_genome_index(
            docker_client,
            genome_index_path,
            star_preamble,
            reference_path
        )
    except:
        file_utils.delete(genome_index_path)
        raise

    dataset = data_handler.datasets.select(experiment.get("dataset"))
    run(docker_client, star_preamble, dataset, destination)
    return destination;

def build_genome_index(docker_client, genome_index_path, preamble, reference_path):
    if not os.path.isdir(genome_index_path):
        file_utils.create_directory(genome_index_path)
        command = preamble + " --runMode genomeGenerate"
        command += " --genomeFastaFiles /{}".format(reference_path)
        docker_client.run(
            "star",
            command
        )

def run(docker_client, preamble, dataset, destination):
    file_utils.create_directory(destination)
    command = preamble + " --readFilesIn"
    for direction, specification in dataset.get("data").items():
        command += " /{}".format(specification["path"])
    docker_client.run(
        "star",
        command
    )

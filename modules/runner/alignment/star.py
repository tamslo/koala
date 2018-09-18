import yaml, os
import modules.file_utils as file_utils

def star(docker_client, destination, data_handler, experiment):
    reference_id = experiment.get("reference")
    dataset = data_handler.datasets.select(experiment.get("dataset"))
    reference_directory = "data/references"
    reference_path = "{}/{}.fa".format(reference_directory, reference_id)
    genome_index_path = "{}/{}_star_index".format(reference_directory, reference_id)
    with open("config.yml", "r") as config_file:
        config = yaml.load(config_file)
        num_threads = int(config["cores"])

    star_preamble = "STAR --runThreadN {} --genomeDir /{} --outFileNamePrefix {}".format(
        num_threads,
        "/" + genome_index_path,
        destination
    )

    try:
        build_genome_index(
            docker_client,
            genome_index_path,
            star_preamble,
            reference_path
        )
    except:
        file_utils.delete(genome_index_path)
        raise

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

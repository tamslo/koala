import yaml, os
import modules.file_utils as file_utils

def __preamble(destination, genome_index_path):
    with open("config.yml", "r") as config_file:
        config = yaml.load(config_file)
        num_threads = int(config["cores"])
    return "STAR --runThreadN {} --genomeDir {} --outFileNamePrefix {}".format(
        num_threads,
        "/" + genome_index_path,
        destination
    )

def build_genome_index(container_name, docker_client, genome_index_path, data_handler, experiment, destination):
    file_utils.create_directory(genome_index_path)
    reference_path = data_handler.reference_path(experiment)
    command = __preamble(destination, genome_index_path) + " --runMode genomeGenerate"
    command += " --genomeFastaFiles /{}".format(reference_path)
    # For debugging and manual execution
    with open(os.path.join(destination, "Commands.txt"), "a") as command_file:
        command_file.write("{}\n".format(command))
    docker_client.run(
        container_name,
        command
    )

def run(container_name, docker_client, dataset, genome_index_path, destination):
    command = __preamble(destination, genome_index_path) + " --readFilesIn"
    for direction, specification in dataset.get("data").items():
        command += " /{}".format(specification["path"])
    # For debugging and manual execution
    with open(os.path.join(destination, "Commands.txt"), "a") as command_file:
        command_file.write("{}\n".format(command))
    docker_client.run(
        container_name,
        command
    )

import yaml, os, shutil
import modules.file_utils as file_utils

def star(docker_client, destination, data_handler, experiment):
    file_utils.create_directory(destination)
    genome_index_path = data_handler.genome_index_path(experiment, "star")
    temp_genome_index_path = genome_index_path + ".running"

    if not os.path.isdir(genome_index_path):
        try:
            build_genome_index(
                docker_client,
                temp_genome_index_path,
                data_handler,
                experiment,
                destination
            )
        except:
            file_utils.delete(temp_genome_index_path)
            raise
        os.rename(temp_genome_index_path, genome_index_path)

    dataset = data_handler.datasets.select(experiment.get("dataset"))
    run(docker_client, dataset, genome_index_path, destination)
    return destination;

def __preamble(destination, genome_index_path):
    with open("config.yml", "r") as config_file:
        config = yaml.load(config_file)
        num_threads = int(config["cores"])
    return "STAR --runThreadN {} --genomeDir {} --outFileNamePrefix {}".format(
        num_threads,
        "/" + genome_index_path,
        destination
    )

def build_genome_index(docker_client, genome_index_path, data_handler, experiment, destination):
    file_utils.create_directory(genome_index_path)
    reference_path = data_handler.reference_path(experiment)
    command = __preamble(destination, genome_index_path) + " --runMode genomeGenerate"
    command += " --genomeFastaFiles /{}".format(reference_path)
    with open(os.path.join(destination, "Commands.txt"), "a") as command_file:
        command_file.write("{}\n".format(command))
    docker_client.run(
        "star",
        command
    )

def run(docker_client, dataset, genome_index_path, destination):
    command = __preamble(destination, genome_index_path) + " --readFilesIn"
    for direction, specification in dataset.get("data").items():
        command += " /{}".format(specification["path"])
    with open(os.path.join(destination, "Commands.txt"), "a") as command_file:
        command_file.write("{}\n".format(command))
    docker_client.run(
        "star",
        command
    )

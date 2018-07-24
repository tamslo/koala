import yaml, os
import modules.file_utils as file_utils

def star(docker_client, destination, experiment):
    genome_index_path = destination + "/index"
    with open("config.yml", "r") as config_file:
        config = yaml.load(config_file)
        num_threads = int(config["cores"])

    star_preamble = "STAR --runThreadN {} --genomeDir /{}".format(
        num_threads,
        genome_index_path
    )
    build_genome_index(
        docker_client,
        genome_index_path,
        star_preamble,
        experiment.get("reference")
    )
    return destination;

def build_genome_index(docker_client, genome_index_path, preamble, reference):
    if not os.path.isdir(genome_index_path):
        os.makedirs(genome_index_path)
        reference_path = "/data/references/{}.fa".format(reference)
        command = preamble + " --runMode genomeGenerate"
        command += " --genomeFastaFiles {}".format(reference_path)
        docker_client.run(
            "star",
            command
        )

def run(docker_client, destination, options):
    # --readFilesIn /data/datasets/... [/data/datasets/...]
    # --outFileNamePrefix /data/datasets/...
    return None

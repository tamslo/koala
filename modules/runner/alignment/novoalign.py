import os
import modules.file_utils as file_utils

def novoalign(docker_client, destination, data_handler, experiment):
    file_utils.create_directory(destination)
    genome_index_path = data_handler.genome_index_path(experiment, "novoalign")
    temp_genome_index_path = genome_index_path + ".running"

    if not os.path.exists(genome_index_path):
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

def build_genome_index(docker_client, genome_index_path, data_handler, experiment, destination):
    reference_id = experiment.get("reference")
    reference_path = data_handler.reference_path(experiment)
    command = "novoindex -n {} {} {}".format(reference_id, genome_index_path, reference_path)
    with open(os.path.join(destination, "Commands.txt"), "a") as command_file:
        command_file.write("{}\n".format(command))
    docker_client.run(
        "novoalign",
        command
    )

def run(docker_client, dataset, genome_index_path, destination):
    command = "novoalign -o SAM -f "
    for direction, specification in dataset.get("data").items():
        command += " /{}".format(specification["path"])
    command += " -d /{}".format(genome_index_path)
    with open(os.path.join(destination, "Commands.txt"), "a") as command_file:
        command_file.write("{}\n".format(command))

    output = docker_client.run(
        "novoalign",
        command,
        log_config={"type": "json-file"},
        detach=False,
        stdout=True,
        stderr=True
    )
    print(output, flush=True)
    out_file_path = os.path.join(destination, "Out.sam")
    # open(out_file_path, "w").close()
    log_file_path = os.path.join(destination, "Out.log")
    # open(log_file_path, "w").close()
    raise Exception("Testing NovoAlign, still needs to be implemented")

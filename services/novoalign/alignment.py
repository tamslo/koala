import os

def build_genome_index(container_name, docker_client, genome_index_path, data_handler, experiment, destination):
    reference_id = experiment.get("reference")
    reference_path = data_handler.reference_path(experiment)
    command = "novoindex -n {} {} {}".format(reference_id, genome_index_path, reference_path)
    # For debugging and manual execution
    with open(os.path.join(destination, "Commands.txt"), "a") as command_file:
        command_file.write("{}\n".format(command))
    docker_client.run(
        container_name,
        command
    )

def run(container_name, docker_client, dataset, genome_index_path, destination):
    command = "novoalign -o SAM -f "
    for direction, specification in dataset.get("data").items():
        command += " /{}".format(specification["path"])
    command += " -d /{}".format(genome_index_path)
    # For debugging and manual execution
    with open(os.path.join(destination, "Commands.txt"), "a") as command_file:
        command_file.write("{}\n".format(command))

    output = docker_client.run(
        container_name,
        command,
        log_config={"type": "json-file"},
        stderr=True
    )
    print(output, flush=True)
    out_file_path = os.path.join(destination, "Out.sam")
    # open(out_file_path, "w").close()
    log_file_path = os.path.join(destination, "Out.log")
    # open(log_file_path, "w").close()
    raise Exception("Testing NovoAlign, still needs to be implemented")

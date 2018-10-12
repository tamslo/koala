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

    container = docker_client.run(
        container_name,
        command,
        stderr=True,
        detach=True,
        auto_remove=False
    )

    while(container.status != "exited"):
        container.reload()

    out_file_path = os.path.join(destination, "Out.sam")
    with open(out_file_path, "wb") as out_file:
        output = container.logs(stdout=True, stderr=False)
        out_file.write(output)

    log_file_path = os.path.join(destination, "Out.log")
    with open(log_file_path, "wb") as log_file:
        output = container.logs(stderr=True, stdout=False)
        print("NovoAlign STDERR output is: {}".format(str(output)))
        log_file.write(output)

    container.kill()
    container.remove()

import os, time

def build_genome_index(container_name, docker_client, genome_index_path, data_handler, experiment, destination):
    reference_path = data_handler.reference_path(experiment)
    command = "touch {}".format(genome_index_path)
    docker_client.run(
        container_name,
        command
    )

def run(container_name, docker_client, dataset, genome_index_path, destination):
    command = "bash ./echo.sh"
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
        log_file.write(output)

    container.remove()

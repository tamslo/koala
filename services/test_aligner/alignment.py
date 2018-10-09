import os

def build_genome_index(container_name, docker_client, genome_index_path, data_handler, experiment, destination):
    reference_path = data_handler.reference_path(experiment)
    command = "touch {}".format(genome_index_path)
    docker_client.run(
        container_name,
        command
    )

def run(container_name, docker_client, dataset, genome_index_path, destination):
    print("Container name: {}".format(container_name))
    command = "bash ./echo.sh"
    output = docker_client.run(
        container_name,
        command,
        log_config={"type": "json-file"},
        stderr=True,
        stdout=False
    )
    print(output, flush=True)
    out_file_path = os.path.join(destination, "Out.sam")
    with open(out_file_path, "w") as out_file:
        out_file.write(str(output))
    log_file_path = os.path.join(destination, "Out.log")
    with open(log_file_path, "w") as log_file:
        log_file.write(str(output))
    raise Exception("Enabling re-run without deleting everything")

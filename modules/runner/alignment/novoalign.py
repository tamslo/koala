def novoalign(docker_client, destination, reference_id, dataset):
    path = destination + "/novoalign.bam"
    docker_client.run(
        "novoalign",
        "touch " + path
    )
    return path

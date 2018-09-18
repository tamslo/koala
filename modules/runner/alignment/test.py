def test(docker_client, destination, reference_id, dataset):
    path = destination + "/dummy.bam"
    docker_client.run(
        "test",
        "touch " + path
    )
    return path

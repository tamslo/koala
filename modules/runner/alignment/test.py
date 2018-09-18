def test(docker_client, destination, data_handler, experiment):
    path = destination + "/dummy.bam"
    docker_client.run(
        "test",
        "touch " + path
    )
    return path

def test(docker_client, destination, experiment):
    path = destination + "/dummy.bam"
    docker_client.run(
        "test",
        "touch " + path
    )
    return path

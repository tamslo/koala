def test(docker_client, destination, experiment):
    docker_client.run(
        "test",
        "touch {}/dummy.bam".format(destination)
    )

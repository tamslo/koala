def test(docker_client, aligner, destination):
    docker_client.run(
        "test",
        "touch {}/dummy.bam".format(destination)
    )

def star(docker_client, aligner, destination):
    docker_client.run(
        "star",
        "touch {}/star-dummy.bam".format(destination)
    )

def align(docker_client, aligner, destination):
    aligner_actions = {
        "test": test,
        "star": star
    }
    aligner_actions[aligner](docker_client, aligner, destination)

import yaml, docker

class Docker:
    def __init__(self, data_directory):
        self.docker_client = docker.from_env()
        with open("config.yml", "r") as config_file:
            config = yaml.load(config_file)
            self.absolute_data_path = config["absolute_repo_path"] + data_directory

    def run(self, container, command, auto_remove=True, **kwargs):
        return self.docker_client.containers.run(
            container,
            command,
            volumes={
                self.absolute_data_path: {
                    "bind": "/data",
                    "mode": "rw"
                }
            },
            auto_remove=auto_remove,
            **kwargs
        )

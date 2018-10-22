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

    def run_and_write_logs(self, container, command, out_file_path, log_file_path=None):
        stderr = True if log_file_path != None else False
        container = self.run(
            container,
            command,
            stderr=stderr,
            detach=True,
            auto_remove=False
        )

        out_file = open(out_file_path, "wb")
        for line in container.logs(stdout=True, stderr=False, stream=True):
            out_file.write(line)
        out_file.close()

        if stderr:
            log_file = open(log_file_path, "wb")
            for line in container.logs(stdout=False, stderr=True, stream=True):
                log_file.write(line)
            log_file.close()

        container.reload()
        if container.status != "exited":
            container.stop()
        container.remove()

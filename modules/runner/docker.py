import yaml, docker, os

class Docker:
    def __init__(self, data_directory):
        self.docker_client = docker.from_env()
        with open("config.yml", "r") as config_file:
            config = yaml.load(config_file)
            self.absolute_data_path = config["absolute_repo_path"] + data_directory

    def run(self, container, command, auto_remove=True, **kwargs):
        if command == None or command == "":
            print("[WARNING] Executing an empty command in {}".format(container), flush=True)

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

    def run_and_write_logs(self, container, command, stdout_file_path, stderr_file_path=None):
        stderr = stderr_file_path != None
        container = self.run(
            container,
            command,
            stderr=stderr,
            detach=True,
            auto_remove=False
        )

        out_file = open(stdout_file_path, "wb")
        for line in container.logs(stdout=True, stderr=False, stream=True):
            out_file.write(line)
        out_file.close()

        if stderr:
            log_file = open(stderr_file_path, "wb")
            for line in container.logs(stdout=False, stderr=True, stream=True):
                log_file.write(line)
            log_file.close()

        container.reload()
        if container.status != "exited":
            container.stop()
        container.remove()

        # If log files are empty, delete them
        log_file_path = stdout_file_path if not stderr else stderr_file_path
        if os.stat(log_file_path).st_size == 0:
            os.remove(log_file_path)

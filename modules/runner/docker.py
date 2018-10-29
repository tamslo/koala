import yaml, docker, os
import modules.file_utils as file_utils

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

    def run_and_write_logs(self, container, command, stdout_file_path=None, stderr_file_path=None):
        stderr = stderr_file_path != None
        container = self.run(
            container,
            command,
            stderr=stderr,
            detach=True,
            auto_remove=False
        )

        if stdout_file_path != None:
            out_file = open(stdout_file_path, "wb")
            for line in container.logs(stdout=True, stderr=False, stream=True):
                out_file.write(line)
            out_file.close()

        if stderr_file_path != None:
            log_file = open(stderr_file_path, "wb")
            for line in container.logs(stdout=False, stderr=True, stream=True):
                log_file.write(line)
            log_file.close()

        container.reload()
        if container.status != "exited":
            container.stop()
        container.remove()

        # If log file is empty, delete it. If only one file is written it is
        # expected to be the log file, if both files are written, it is expected
        # to be stderr_file_path
        log_file_path = stderr_file_path or stdout_file_path
        if not file_utils.file_has_content(log_file_path):
            file_utils.delete(log_file_path)

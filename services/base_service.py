import os
import yaml

class BaseService:
    def __init__(self, config_path, service_id):
        with open(config_path, "r") as config_file:
            config = yaml.load(config_file)
            self.id = service_id
            for key in config:
                setattr(self, key, config[key])

    # self.command needs to be implemented by subclasses or self.run overwritten
    def run(self, parameters):
        if hasattr(self, "command"):
            return self.run_docker(parameters, self.command(parameters))

    def frontend_information(self):
        return { "id": self.id, "name": self.name, "type": self.type }

    # For debugging and manual execution
    def __log_command(self, parameters, command):
        destination = parameters["destination"]
        with open(os.path.join(destination, "Commands.txt"), "a") as command_file:
            command_file.write("{}\n".format(command))

    def __image_name(self):
        if hasattr(self, "image"):
            return self.image
        else:
            return self.id

    def run_docker(self, parameters, command, log_to_output=False, rename_output=False):
        self.__log_command(parameters, command)
        docker_client = parameters["docker_client"]
        destination = parameters["destination"]
        stdout_file_path = os.path.join(destination, "Out.log")
        stderr_file_path = None

        if log_to_output:
            stderr_file_path = stdout_file_path # Write logs from stderr
            stdout_file_path = os.path.join(destination, parameters["out_file_name"])

        docker_client.run_and_write_logs(
            self.__image_name(),
            command,
            stdout_file_path,
            stderr_file_path
        )

        if rename_output:
            os.rename(
                os.path.join(parameters["destination"], self.output_file_name),
                os.path.join(parameters["destination"], parameters["out_file_name"])
            )

import os
import yaml

class BaseService:
    def __init__(self, config_path, service_id, image_name):
        with open(config_path, "r") as config_file:
            config = yaml.load(config_file)
            self.id = service_id
            self.image = image_name
            for key in config:
                setattr(self, key, config[key])

    # Gets parameters and builds command. Possibly adds output_parameters and
    # runs docker `self.run_docker(command, parameters, [, output_parameters])`.
    def run(self, parameters):
        raise Exception("Method base_service.run needs to be implemented by subclasses")

    def frontend_information(self):
        return { "id": self.id, "name": self.name, "type": self.type }

    def __log_path(self, parameters):
        destination = parameters["destination"]
        return os.path.join(destination, "Runstats.txt")

    def __image_name(self, parameters):
        if "docker_image" in parameters:
            return parameters["docker_image"]
        else:
            return self.image

    def run_docker(self, command, parameters, output_parameters={}):
        docker_client = parameters["docker_client"]
        destination = parameters["destination"]

        # Set defaults for output_parameters
        def set_default(dict, key, default):
            return dict[key] if key in dict else default
        log_is_output = set_default(output_parameters, "log_is_output", False)
        log_file_path = set_default(
            output_parameters,
            "log_file_path",
            os.path.join(destination, "Out.log")
        )

        # Default. The command automatically writes to file, write to log what
        # is printed in STDOUT and STDERR.
        stdout_file_path = log_file_path
        stderr_file_path = log_file_path

        # The command does not write to file, write to output file what is
        # printed in STDOUT and to log what is printed in STDERR.
        if log_is_output:
            stdout_file_path = output_parameters["out_file_path"]
            stderr_file_path = log_file_path

        docker_client.run_and_write_logs(
            self.__image_name(parameters),
            command,
            stdout_file_path,
            stderr_file_path,
            self.__log_path(parameters)
        )

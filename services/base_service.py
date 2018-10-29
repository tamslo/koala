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

    # For debugging and manual execution
    def __log_command(self, command, parameters):
        destination = parameters["destination"]
        with open(os.path.join(destination, "Commands.txt"), "a") as command_file:
            command_file.write("{}\n".format(command))

    def __image_name(self, parameters):
        if "docker_image" in parameters:
            return parameters["docker_image"]
        else:
            return self.image

    def run_docker(self, command, parameters, output_parameters={}):
        self.__log_command(command, parameters)
        docker_client = parameters["docker_client"]
        destination = parameters["destination"]

        # Set defaults for output_parameters
        def set_default(dict, key, default):
            return dict[key] if key in dict else default
        log_is_output = set_default(output_parameters, "log_is_output", False)
        log_file_path = set_default(output_parameters, "log_file_path", None)

        # Default. The command automatically writes to file, write to log what
        # is printed in STDOUT.
        stdout_file_path = log_file_path or os.path.join(destination, "Out.log")
        stderr_file_path = None

        # The command does not write to file, write to output file what is
        # printed in STDOUT and to log what is printed in STDERR.
        if log_is_output:
            stderr_file_path = stdout_file_path # Write logs from stderr
            stdout_file_path = os.path.join(destination, parameters["out_file_name"])
        else:
            # The command automatically writes to file, write to log what is
            # printed in STDERR.
            if "log_from_stderr" in parameters and parameters["log_from_stderr"]:
                stderr_file_path = stdout_file_path
                stdout_file_path = None

        docker_client.run_and_write_logs(
            self.__image_name(parameters),
            command,
            stdout_file_path,
            stderr_file_path
        )            

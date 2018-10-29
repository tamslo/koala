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
    def __log_command(self, parameters, command):
        destination = parameters["destination"]
        with open(os.path.join(destination, "Commands.txt"), "a") as command_file:
            command_file.write("{}\n".format(command))

    def __image_name(self, parameters):
        if "docker_image" in parameters:
            return parameters["docker_image"]
        else:
            return self.image

    def run_docker(self, command, parameters, log_is_output=False, rename_output=False, log_file_path=None):
        self.__log_command(parameters, command)
        docker_client = parameters["docker_client"]
        destination = parameters["destination"]
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

        out_file_path = None
        if "out_file_name" in parameters:
            out_file_path = os.path.join(parameters["destination"], parameters["out_file_name"])
            if rename_output:
                os.rename(
                    os.path.join(parameters["destination"], self.output_file_name),
                    out_file_path
                )

            file_exists = os.path.exists(out_file_path)
            file_is_empty = file_exists and os.stat(out_file_path).st_size == 0
            if file_is_empty or not file_exists:
                raise Exception("Output not written")

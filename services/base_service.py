import os
import yaml

class BaseService:
    def __init__(self, config_path, service_id):
        with open(config_path, "r") as config_file:
            config = yaml.load(config_file)
            self.id = service_id
            for key in config:
                setattr(self, key, config[key])

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

    def run_main(self, parameters, command):
        if self.creates_output_files:
            self.run_simple(parameters, command)
        else:
            self.run_with_file_creation(parameters, command)

    def run_simple(self, parameters, command):
        self.__log_command(parameters, command)
        parameters["docker_client"].run(
            self.__image_name(),
            command
        )
        if "out_file_name" in parameters:
            os.rename(
                os.path.join(parameters["destination"], self.output_file_name),
                os.path.join(parameters["destination"], parameters["out_file_name"])
            )

    def run_with_file_creation(self, parameters, command):
        self.__log_command(parameters, command)
        docker_client = parameters["docker_client"]
        destination = parameters["destination"]

        container = docker_client.run(
            self.__image_name(),
            command,
            stderr=True,
            detach=True,
            auto_remove=False
        )

        out_file_path = os.path.join(destination, parameters["out_file_name"])
        out_file = open(out_file_path, "wb")
        for line in container.logs(stdout=True, stderr=False, stream=True):
            out_file.write(line)
        out_file.close()

        log_file_path = os.path.join(destination, "Out.log")
        log_file = open(log_file_path, "wb")
        for line in container.logs(stdout=False, stderr=True, stream=True):
            log_file.write(line)
        log_file.close()

        container.reload()
        if container.status != "exited":
            container.stop()
        container.remove()

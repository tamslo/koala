import os
import modules.file_utils as file_utils
from .base_service import BaseService

class BaseAligner(BaseService):
    # NOTE: TO BE IMPLEMENTED BY SPECIFIC ALIGNER
    # Returns: Command to build genome index as string
    def build_index_command(self, parameters):
        return None

    # NOTE: TO BE IMPLEMENTED BY SPECIFIC ALIGNER
    # Returns: Command to run alignment as string
    def alignment_command(self, parameters):
        return None

    # For debugging and manual execution
    def __log_command(self, parameters, command):
        destination = parameters["destination"]
        with open(os.path.join(destination, "Commands.txt"), "a") as command_file:
            command_file.write("{}\n".format(command))

    def __run_simple(self, parameters, command):
        self.__log_command(parameters, command)
        parameters["docker_client"].run(
            self.id,
            command
        )

    def __run_with_file_creation(self, parameters, command):
        docker_client = parameters["docker_client"]
        destination = parameters["destination"]
        container = docker_client.run(
            self.id,
            command,
            stderr=True,
            detach=True,
            auto_remove=False
        )

        # while(container.status != "exited"):
        #     container.reload()

        out_file_path = os.path.join(destination, "Out.sam")
        out_file = open(out_file_path, "wb")
        for line in container.logs(stdout=True, stderr=False, stream=True):
            print("STDIN (.sam)", flush=True)
            print(line, flush=True)
            out_file.write(line)
        out_file.close()

        log_file_path = os.path.join(destination, "Out.log")
        log_file = open(log_file_path, "wb")
        for line in container.logs(stdout=False, stderr=True, stream=True):
            print("STDERR (.log)", flush=True)
            print(line, flush=True)
            log_file.write(line)
        log_file.close()

        container.remove()

    def build_genome_index(self, parameters):
        command = self.build_index_command(parameters)
        self.__log_command(parameters, command)

        if self.reference_is_directory:
            file_utils.create_directory(parameters["genome_index_path"])

        self.__run_simple(parameters, command)

    def align(self, parameters):
        command = self.alignment_command(parameters)
        self.__log_command(parameters, command)

        if self.creates_output_files:
            self.__run_simple(parameters, command)
        else:
            self.__run_with_file_creation(parameters, command)

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

    def build_genome_index(self, parameters):
        if self.reference_is_directory:
            file_utils.create_directory(parameters["genome_index_path"])
        command = self.build_index_command(parameters)
        self.run(parameters, command)

    def align(self, parameters):
        write_logs = not self.creates_output_files
        command = self.alignment_command(parameters)
        self.run(
            parameters,
            command,
            write_logs=write_logs,
            rename_output=self.creates_output_files
        )

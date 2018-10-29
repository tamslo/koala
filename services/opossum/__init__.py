import modules.file_utils as file_utils
from ..base_service import BaseService

class Opossum(BaseService):
    def run(self, parameters):
        out_file_path = parameters["destination"] + "Out.bam"
        command = "python Opossum.py --BamFile=/{} --OutFile=/{}".format(
            parameters["experiment"].get_input_directory(self.id) + "Out.bam",
            out_file_path
        )
        self.run_docker(command, parameters)
        file_utils.validate_file_content(out_file_path)

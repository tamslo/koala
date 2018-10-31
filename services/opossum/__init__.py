import modules.file_utils as file_utils
from ..base_service import BaseService

class Opossum(BaseService):
    def run(self, parameters):
        experiment = parameters["experiment"]
        preceding_service = experiment.get_preceding_service(self.id)
        in_file_path = experiment.get_input_directory(self.id) + "Out.bam"
        out_file_path = parameters["destination"] + "Out.bam"
        soft_clips_exist = preceding_service.soft_clips_exist
        command = "python Opossum.py --BamFile=/{} --OutFile=/{} --SoftClipsExist={}".format(
            in_file_path,
            out_file_path,
            soft_clips_exist
        )
        self.run_docker(command, parameters)
        file_utils.validate_file_content(out_file_path)

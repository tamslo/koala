import os
import modules.file_utils as file_utils
from ..base_service import BaseService

class HaplotypeCaller(BaseService):
    def run(self, parameters):
        experiment = parameters["experiment"]
        destination = parameters["destination"]
        data_handler = parameters["data_handler"]
        docker_client = parameters["docker_client"]

        reference_path = data_handler.reference_path(experiment)
        reference_index_path = reference_path + ".fai"
        reference_dict_path = data_handler.reference_path(
            experiment,
            alternate_file_ending=".dict"
        )

        # Generate index if not there
        if not os.path.exists(reference_index_path):
            command = "samtools faidx /{}".format(reference_path)
            output_parameters = {
                "log_file_path": destination + "Index.log",
                "log_from_stderr": True
            }
            self.run_docker(command, parameters, output_parameters)

        # Generate dict if not there
        if not os.path.exists(reference_dict_path):
            command = "gatk CreateSequenceDictionary -R /{} -O /{}".format(
                reference_path,
                reference_dict_path
            )
            output_parameters = {
                "log_file_path": destination + "Dict.log",
                "log_from_stderr": True
            }
            self.run_docker(command, parameters, output_parameters)

        # Run variant calling
        out_file_path = destination + "Out.vcf"
        command = "gatk HaplotypeCaller -I /{} -O /{} -R /{} -dontUseSoftClippedBases -stand_call_conf 20".format(
            experiment.get_input_directory(self.id) + "Out.bam",
            out_file_path,
            data_handler.reference_path(experiment)
        )
        output_parameters = { "log_from_stderr": True }
        self.run_docker(command, parameters, output_parameters)
        file_utils.validate_file_content(out_file_path)

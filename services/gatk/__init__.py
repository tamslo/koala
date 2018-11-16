import os
import modules.file_utils as file_utils
from ..base_service import BaseService

class GatkFilters(BaseService):
    def run(self, parameters):
        destination = parameters["destination"]
        experiment = parameters["experiment"]
        data_handler = parameters["data_handler"]

        in_file_path = experiment.get_input_directory(self.id) + "Out.bam"
        reference_path = data_handler.reference_path(experiment)

        # Remove duplicates
        deduplicated_path = destination + "Deduplicated.bam"
        metrics_path = destination + "Deduplicate.metrics"
        command = "gatk MarkDuplicates -I /{} -O /{} -M /{} " \
            "--VALIDATION_STRINGENCY=SILENT".format(
                in_file_path,
                deduplicated_path,
                metrics_path
            )
        output_parameters = {
            "log_file_path": destination + "Deduplicate.log",
            "log_from_stderr": True
        }
        self.run_docker(command, parameters, output_parameters)
        file_utils.validate_file_content(deduplicated_path)

        # Remove introns
        out_file_path = destination + "Out.bam"
        command = "gatk SplitNCigarReads -R /{} -I /{} -O /{}".format(
                reference_path,
                deduplicated_path,
                out_file_path
            )
        output_parameters = { "log_file_path": destination + "SplitN.log" }
        self.run_docker(command, parameters, output_parameters)
        file_utils.validate_file_content(out_file_path)
        file_utils.delete(deduplicated_path)

class HaplotypeCaller(BaseService):
    def run(self, parameters):
        experiment = parameters["experiment"]
        destination = parameters["destination"]
        data_handler = parameters["data_handler"]
        docker_client = parameters["docker_client"]

        reference_path = data_handler.reference_path(experiment)
        # Run variant calling
        out_file_path = destination + "Out.vcf"
        command = "gatk HaplotypeCaller -I /{} -O /{} -R /{} " \
            "--dont-use-soft-clipped-bases " \
            "--standard-min-confidence-threshold-for-calling 20".format(
            experiment.get_input_directory(self.id) + "Out.bam",
            out_file_path,
            data_handler.reference_path(experiment)
        )
        output_parameters = { "log_from_stderr": True }
        self.run_docker(command, parameters, output_parameters)
        file_utils.validate_file_content(out_file_path)

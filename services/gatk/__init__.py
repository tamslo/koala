import os
import modules.file_utils as file_utils
from ..base_service import BaseService

class GatkVariantFiltration(BaseService):
    def run(self, parameters):
        destination = parameters["destination"]
        experiment = parameters["experiment"]
        data_handler = parameters["data_handler"]

        in_file_path = experiment.get_input_directory(self.id) + "Out.vcf"
        reference_path = data_handler.reference_path(experiment)
        out_file_path = destination + "Out.vcf"

        command = "gatk VariantFiltration -R /{} -V /{} -O /{} " \
            "-window 35 -cluster 3 --filter-name FS -filter 'FS > 30.0' " \
            " --filter-name QD -filter 'QD < 2.0'".format(
                reference_path,
                in_file_path,
                out_file_path
            )
        self.run_docker(command, parameters)
        file_utils.validate_file_content(out_file_path)

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
        output_parameters = { "log_file_path": destination + "Deduplicate.log" }
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
    # Can be overwritten by specific components
    def add_filters(self, command):
        return command

    def run(self, parameters, in_file_path=None):
        experiment = parameters["experiment"]
        destination = parameters["destination"]
        data_handler = parameters["data_handler"]
        docker_client = parameters["docker_client"]

        reference_path = data_handler.reference_path(experiment)
        # Run variant calling
        in_file_path = in_file_path or experiment.get_input_directory(self.id) + "Out.bam"
        out_file_path = destination + "Out.vcf"
        command = "gatk HaplotypeCaller -I /{} -O /{} -R /{} " \
            "--dont-use-soft-clipped-bases " \
            "--standard-min-confidence-threshold-for-calling 20".format(
            in_file_path,
            out_file_path,
            data_handler.reference_path(experiment)
        )
        command = self.add_filters(command)
        self.run_docker(command, parameters)
        file_utils.validate_file_content(out_file_path)

class HaplotypeCallerFiltered(HaplotypeCaller):
    def add_filters(self, command):
        return command + " -L {}".format(",".join(self.chromosomes))

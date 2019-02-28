import os
import modules.file_utils as file_utils
from ..base_service import BaseService

class GiabEvaluator(BaseService):
    def transcriptome_regions_path(self, alignment_path, parameters):
        transcriptome_regions_path = alignment_path + "aligned_coverage_regions.bed"
        if not os.path.exists(transcriptome_regions_path):
            bam_path = alignment_path + "Out.bam"
            coverage_path = alignment_path + "Out.base_coverage"
            min_coverage = 2

            # Create coverage file
            command = "bedtools genomecov -d -ibam /{}".format(bam_path)
            output_parameters = {
                "log_is_output": True,
                "out_file_path": coverage_path,
                "log_file_path": parameters["destination"] + "Coverage.log"
            }
            self.run_docker(command, parameters, output_parameters)
            file_utils.validate_file_content(coverage_path)

            # Create BED from coverage file
            command = "python base_coverage_to_bed.py /{} {} /{}".format(
                coverage_path,
                str(min_coverage),
                transcriptome_regions_path
            )
            self.run_docker(command, parameters, log_file_name="CoverageToBed.log")
            file_utils.validate_file_content(transcriptome_regions_path)

        return transcriptome_regions_path

    def bedtools(self, function, a_file_path, b_file_path, out_file_path, parameters, options=""):
        destination = parameters["destination"]
        log_file_path = destination + function.capitalize() + ".log"
        command = "bedtools {} " \
            "-a /{} " \
            "-b /{} {}".format(function, a_file_path, b_file_path, options)
        output_parameters = {
            "log_is_output": True,
            "out_file_path": out_file_path,
            "log_file_path": log_file_path
        }
        self.run_docker(command, parameters, output_parameters)

    def run(self, parameters):
        experiment = parameters["experiment"]
        reference_id = experiment.get("reference")
        destination = parameters["destination"]
        vcf_file_path = destination + "Out.vcf"
        alignment_path = experiment.get("pipeline")["alignment"]["directory"]
        confidence_regions_path = alignment_path + "confidence_calls.bed".format(reference_id)

        # Intersect confidence regions with transcriptome regions if not already done
        if not os.path.exists(confidence_regions_path):
            confidence_genome_regions_path = "data/giab/{}/confidence_calls.bed".format(reference_id)
            transcriptome_regions_path = self.transcriptome_regions_path(alignment_path, parameters)
            self.bedtools(
                "intersect",
                confidence_genome_regions_path,
                transcriptome_regions_path,
                confidence_regions_path,
                parameters
            )
            file_utils.validate_file_content(confidence_regions_path)


        # Filter data if necessary
        action_handler = parameters["action_handler"]
        additional_commands = ""
        if hasattr(action_handler, "chromosomes"):
            # Escape spaces for bash
            space_escape = "%%"
            additional_commands = "--location{}{}".format(
                space_escape,
                ",".join(action_handler.chromosomes)
            )

        command = "./hap.py /data/giab/{0}/confidence_calls.vcf /{1}Out.vcf " \
            "-f /{2} " \
            "-o /{1}Evaluation " \
            "-r /data/references/{0}.fa " \
            "--location {3}".format(
                reference_id,
                destination,
                confidence_regions_path,
                additional_commands
            )
        output_parameters = { "log_file_path": destination + "Evaluation.log" }
        self.run_docker(command, parameters, output_parameters)

        for file_name in os.listdir(destination):
            if file_name.startswith("Evaluation"):
                file_path = destination + file_name
                if not file_utils.file_has_content(file_path):
                    file_utils.delete(file_path)

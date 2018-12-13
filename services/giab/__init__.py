import os
import modules.file_utils as file_utils
from ..base_service import BaseService

class GiabEvaluator(BaseService):
    def run(self, parameters):
        experiment = parameters["experiment"]
        reference_id = experiment.get("reference")
        destination = parameters["destination"]
        vcf_file_path = destination + "Out.vcf"

        # Intersect confidence regions with coding regions if not already done
        confidence_regions_path = "data/giab/{}/confidence_calls_coding.bed".format(reference_id)
        if not os.path.exists(confidence_regions_path):
            command = "bedtools intersect " \
                "-a /data/giab/{0}/confidence_calls.bed " \
                "-b /data/annotations/{0}_coding_exons.bed".format(reference_id)
            output_parameters = {
                "log_is_output": True,
                "out_file_path": confidence_regions_path,
                "log_file_path": destination + "Intersect.log"
            }
            self.run_docker(command, parameters, output_parameters)

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
                ",".join(action_handler.chromosomes)
            )
        output_parameters = { "log_file_path": destination + "Evaluation.log" }
        self.run_docker(command, parameters, output_parameters)

        for file_name in os.listdir(destination):
            if file_name.startswith("Evaluation"):
                file_path = destination + file_name
                if not file_utils.file_has_content(file_path):
                    file_utils.delete(file_path)

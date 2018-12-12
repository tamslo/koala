import os
import modules.file_utils as file_utils
from ..base_service import BaseService

class GiabEvaluator(BaseService):
    def run(self, parameters):
        experiment = parameters["experiment"]
        reference_id = experiment.get("reference")
        destination = parameters["destination"]
        vcf_file_path = destination + "Out.vcf"

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
            "-f /data/giab/{0}/confidence_calls.bed " \
            "-o /{1}Evaluation " \
            "-r /data/references/{0}.fa " \
            "--location {2}".format(
                reference_id,
                destination,
                ",".join(action_handler.chromosomes)
            )
        output_parameters = { "log_file_path": destination + "Evaluation.log" }
        self.run_docker(command, parameters, output_parameters)

        for file_name in os.listdir(destination):
            if file_name.startswith("Evaluation"):
                file_path = destination + file_name
                if not file_utils.file_has_content(file_path):
                    file_utils.delete(file_path)

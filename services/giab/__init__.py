import os
import modules.file_utils as file_utils
from ..base_service import BaseService

class GiabEvaluator(BaseService):
    def run(self, parameters):
        experiment = parameters["experiment"]
        reference_id = experiment.get("reference")
        destination = parameters["destination"]
        path_prefix = destination
        if path_prefix.endswith("/"):
            path_prefix = path_prefix[:-1] # Trim trailing slash

        # Filter data if necessary
        action_handler = parameters["action_handler"]
        additional_commands = ""
        if hasattr(action_handler, "chromosomes"):
            additional_commands = "-l {}".format(
                " ".join(action_handler.chromosomes)
            )

        command = "bash evaluate_variants.sh /{} {} {} {}".format(
            path_prefix,
            "Out.vcf",
            reference_id,
            additional_commands
        )
        output_parameters = { "log_file_path": destination + "Evaluation.log" }
        self.run_docker(command, parameters, output_parameters)

        for file_name in os.listdir(destination):
            if file_name.startswith("Evaluation"):
                file_path = destination + file_name
                if not file_utils.file_has_content(file_path):
                    file_utils.delete(file_path)

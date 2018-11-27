import modules.file_utils as file_utils
from ..base_service import BaseService

class BeersEvaluator(BaseService):
    def run(self, parameters):
        dataset = parameters["dataset"]
        destination = parameters["destination"]
        command = "bash evaluate_alignment.sh {} {} /{}".format(
            dataset.get("readLength"),
            destination,
            dataset.get("evaluation")["truth_file"]["path"]
        )
        output_parameters = { "log_file_path": destination + "Evaluation.log" }
        self.run_docker(command, parameters, output_parameters)

        for file_name in ["Evaluation.multi.txt", "Evaluation.txt"]:
            file_path = destination + file_name
            if not file_utils.file_has_content(file_path):
                file_utils.delete(file_path)

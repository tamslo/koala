from ..base_service import BaseService

class GiabEvaluator(BaseService):
    def run(self, parameters):
        dataset = parameters["dataset"]
        destination = parameters["destination"]
        path_prefix = destination
        if path_prefix.endswith("/"):
            path_prefix = path_prefix[:-1] # Trim trailing slash
        command = "bash evaluate_variants.sh /{} {} {}".format(
            destination,
            "Out.vcf",
            dataset.get("reference")
        )
        output_parameters = {
            "log_from_stderr": True,
            "log_file_path": destination + "Evaluation.log"
        }
        self.run_docker(command, parameters, output_parameters)

        # for file_name in ["Evaluation.multi.txt", "Evaluation.txt"]:
        #     file_path = destination + file_name
        #     if not file_utils.file_has_content(file_path):
        #         file_utils.delete(file_path)

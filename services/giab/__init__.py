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
        filters = []
        filter_postfix = ""
        if hasattr(action_handler, "filters"):
            filters = action_handler.filters
            giab_path_prefix = "/giab/{}".format(reference_id)

            # Filter BED and VCF file
            command = "awk '"
            for filter, index in enumerate(filters):
                filter_postfix += "_" + filter
                if reference_id == "hg19":
                    filter = filter.replace("chr", "")
                if index != 0:
                    command += " || "
                command += "/^{}/".format(filter)
            command += "' {}/confidence_calls{}".format(giab_path_prefix, file_ending)

            for file_ending in [".bed", ".vcf"]:
                output_parameters = {
                    "log_is_output": True,
                    "out_file_path": "{}/confidence_calls{}{}".format(
                        giab_path_prefix,
                        filter_postfix,
                        file_ending
                    ),
                    "log_file_path": "Filter{}.log".format(file_ending)
                }
                self.run_docker(command, parameters, output_parameters)

        command = "bash evaluate_variants.sh /{} {} {} {}".format(
            path_prefix,
            "Out.vcf",
            reference_id,
            filter_postfix
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

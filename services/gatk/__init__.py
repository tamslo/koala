import os, time, datetime
import modules.file_utils as file_utils
from ..base_service import BaseService

class HaplotypeCaller(BaseService):
    def run(self, parameters):
        experiment = parameters["experiment"]
        destination = parameters["destination"]
        data_handler = parameters["data_handler"]
        docker_client = parameters["docker_client"]
        runtime_log_path = parameters["runtime_log_path"]

        reference_path = data_handler.reference_path(experiment)
        reference_index_path = reference_path + ".fai"
        reference_dict_path = data_handler.reference_path(
            experiment,
            alternate_file_ending=".dict"
        )

        # Generate index if not there
        if not os.path.exists(reference_index_path):
            start_time = time.time()
            command = "samtools faidx /{}".format(reference_path)
            output_parameters = {
                "log_file_path": destination + "Index.log",
                "log_from_stderr": True
            }
            self.run_docker(command, parameters, output_parameters)
            with open(runtime_log_path, "a") as runtime_log:
                runtime = str(datetime.timedelta(seconds=time.time() - start_time))
                runtime_log.write("Index generation: {}\n".format(runtime))
        else:
            with open(runtime_log_path, "a") as runtime_log:
                runtime_log.write("Index already present\n")

        # Generate dict if not there
        if not os.path.exists(reference_dict_path):
            start_time = time.time()
            command = "gatk CreateSequenceDictionary -R /{} -O /{}".format(
                reference_path,
                reference_dict_path
            )
            output_parameters = {
                "log_file_path": destination + "Dict.log",
                "log_from_stderr": True
            }
            self.run_docker(command, parameters, output_parameters)
            with open(runtime_log_path, "a") as runtime_log:
                runtime = str(datetime.timedelta(seconds=time.time() - start_time))
                runtime_log.write("Dict generation: {}\n".format(runtime))
        else:
            with open(runtime_log_path, "a") as runtime_log:
                runtime_log.write("Dict already present\n")

        # Run variant calling
        out_file_path = destination + "Out.vcf"
        command = "gatk HaplotypeCaller -I /{} -O /{} -R /{}".format(
            experiment.get_input_directory(self.id) + "Out.bam",
            out_file_path,
            data_handler.reference_path(experiment)
        )
        output_parameters = { "log_from_stderr": True }
        self.run_docker(command, parameters, output_parameters)
        file_utils.validate_file_content(out_file_path)

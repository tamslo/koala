import os, time, datetime
from ..base_service import BaseService

class HaplotypeCaller(BaseService):
    def run(self, parameters):
        experiment = parameters["experiment"]
        destination = parameters["destination"]
        data_handler = parameters["data_handler"]
        docker_client = parameters["docker_client"]
        runtime_log_path = parameters["runtime_log_path"]
        parameters["log_from_stderr"] = True

        reference_path = data_handler.reference_path(experiment)
        reference_index_path = reference_path + ".fai"
        reference_dict_path = data_handler.reference_path(
            experiment,
            alternate_file_ending=".dict"
        )

        if not os.path.exists(reference_index_path):
            # TODO generate index
            start_time = time.time()
            command = "samtools faidx /{}".format(reference_path)
            log_file_path = destination + "Index.log"
            self.run_docker(command, parameters, log_file_path=log_file_path)
            with open(runtime_log_path, "a") as runtime_log:
                runtime = str(datetime.timedelta(seconds=time.time() - start_time))
                runtime_log.write("Index generation: {}\n".format(runtime))
        else:
            with open(runtime_log_path, "a") as runtime_log:
                runtime_log.write("Index already present\n")

        if not os.path.exists(reference_dict_path):
            # TODO generate index
            start_time = time.time()
            command = "gatk CreateSequenceDictionary -R /{} -O /{}".format(
                reference_path,
                reference_dict_path
            )
            log_file_path = destination + "Dict.log"
            self.run_docker(command, parameters, log_file_path=log_file_path)
            with open(runtime_log_path, "a") as runtime_log:
                runtime = str(datetime.timedelta(seconds=time.time() - start_time))
                runtime_log.write("Dict generation: {}\n".format(runtime))
        else:
            with open(runtime_log_path, "a") as runtime_log:
                runtime_log.write("Index already present\n")

        super().run(parameters)

    def command(self, parameters):
        experiment = parameters["experiment"]
        destination = parameters["destination"]
        data_handler = parameters["data_handler"]
        # Tool does not throw an error. Test in BaseService whether the output
        # file is written by setting the out_file_name parameter.
        parameters["out_file_name"] = "Out.vcf"

        return "gatk HaplotypeCaller -I /{} -O /{} -R /{}".format(
            experiment.get_input_directory(self.id) + "Out.bam",
            destination + parameters["out_file_name"],
            data_handler.reference_path(experiment)
        )

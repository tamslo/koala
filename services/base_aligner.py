import os
import time
import datetime
import modules.file_utils as file_utils
from .base_service import BaseService

class BaseAligner(BaseService):
    # Returns: Command to build genome index as string
    def build_index_command(self, parameters):
        raise Exception("Method base_aligner.build_index_command needs to be implemented by subclasses")

    # Returns: Command to run alignment as string
    def alignment_command(self, parameters):
        raise Exception("Method base_aligner.alignment_command needs to be implemented by subclasses")

    def run(self, parameters):
        docker_client = parameters["docker_client"]
        data_handler = parameters["data_handler"]
        experiment = parameters["experiment"]
        destination = parameters["destination"]
        runtime_log_path = parameters["runtime_log_path"]

        # Define genome index path and temp path (will be renamed if successful)
        genome_index_path = data_handler.genome_index_path(experiment, self.id)
        temp_genome_index_path = genome_index_path + ".running"

        # If neccessary, build genome index
        if not os.path.exists(genome_index_path):
            try:
                index_parameters = {
                    "docker_client": docker_client,
                    "destination": destination,
                    "genome_index_path": temp_genome_index_path,
                    "reference_id": experiment.get("reference"),
                    "reference_path": data_handler.reference_path(experiment)
                }
                build_genome_index_start = time.time()
                self.build_genome_index(index_parameters)
                with open(runtime_log_path, "a") as runtime_log:
                    runtime = str(datetime.timedelta(
                        seconds=time.time() - build_genome_index_start)
                    )
                    runtime_log.write("Index generation: {}\n".format(runtime))
            except:
                file_utils.delete(temp_genome_index_path)
                raise
            os.rename(temp_genome_index_path, genome_index_path)
        else:
            with open(runtime_log_path, "a") as runtime_log:
                runtime_log.write("Index already present\n")

        # Run alignment
        out_file_name = "Out.sam"
        alignment_parameters = {
            "docker_client": docker_client,
            "destination": destination,
            "genome_index_path": genome_index_path,
            "dataset": data_handler.datasets.select(experiment.get("dataset")),
            "out_file_name": out_file_name
        }
        run_start = time.time()
        self.align(alignment_parameters)
        with open(runtime_log_path, "a") as runtime_log:
            runtime = str(datetime.timedelta(seconds=time.time() - run_start))
            runtime_log.write("Alignment: {}\n".format(runtime))

        # Create sorted BAM file from SAM file

        conversion_start = time.time()
        sam_path = destination + out_file_name

        # Somehow, samtools needs a present file to write to but logs to
        # STDOUT anyways, this file will be deleted later
        dummy_file_path = destination + "tmp.file"
        open(dummy_file_path, "w").close()

        conversion_parameters = {
            "docker_client": docker_client,
            "docker_image": "gatk",
            "destination": destination,
            "out_file_name": "Out.bam"
        }
        conversion_output_parameters = {
            "log_is_output": True,
            "log_file_path": destination + "Samtools.log"
        }
        command = "samtools sort -o /{} /{}".format(dummy_file_path, sam_path)

        self.run_docker(
            command,
            conversion_parameters,
            conversion_output_parameters
        )
        os.remove(dummy_file_path)
        file_utils.validate_file_content(destination + conversion_parameters["out_file_name"])
        with open(runtime_log_path, "a") as runtime_log:
            runtime = str(datetime.timedelta(seconds=time.time() - conversion_start))
            runtime_log.write("Convert to sorted BAM: {}\n".format(runtime))

    def build_genome_index(self, parameters):
        if self.reference_is_directory:
            file_utils.create_directory(parameters["genome_index_path"])
        command = self.build_index_command(parameters)
        self.run_docker(command, parameters)

    def align(self, parameters):
        command = self.alignment_command(parameters)
        output_parameters = {
            "log_is_output": not self.creates_output_files,
            "rename_output": self.creates_output_files
        }
        self.run_docker(
            command,
            parameters,
            output_parameters
        )

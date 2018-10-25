import os
import time
import datetime
import modules.file_utils as file_utils
from .base_service import BaseService

class BaseAligner(BaseService):
    # NOTE: TO BE IMPLEMENTED BY SPECIFIC ALIGNER
    # Returns: Command to build genome index as string
    def build_index_command(self, parameters):
        return None

    # NOTE: TO BE IMPLEMENTED BY SPECIFIC ALIGNER
    # Returns: Command to run alignment as string
    def alignment_command(self, parameters):
        return None

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

        # Create BAM file from SAM file
        sam_path = destination + out_file_name
        bam_path = destination + "Out.bam"
        docker_client.run_and_write_logs(
            "gatk",
            "samtools view -Sb /{}".format(sam_path),
            bam_path
        )

    def build_genome_index(self, parameters):
        if self.reference_is_directory:
            file_utils.create_directory(parameters["genome_index_path"])
        command = self.build_index_command(parameters)
        self.run_docker(parameters, command)

    def align(self, parameters):
        write_logs = not self.creates_output_files
        command = self.alignment_command(parameters)
        self.run_docker(
            parameters,
            command,
            write_logs=write_logs,
            rename_output=self.creates_output_files
        )

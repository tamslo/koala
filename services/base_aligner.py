import os
import time
import datetime
import modules.file_utils as file_utils
from .base_service import BaseService

class BaseAligner(BaseService):
    # Optional preparation of building the index, implement in specific class
    def prepare_indexing(self, parameters):
        return None

    # Optional post-processing of alignment, implement in specific class
    def conclude_alignment(self, parameters, out_file_path):
        return None

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

        sam_file_path = destination + "Out.sam"
        bam_file_path = destination + "Out.bam"

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
        alignment_parameters = {
            "docker_client": docker_client,
            "destination": destination,
            "genome_index_path": genome_index_path,
            "dataset": data_handler.datasets.select(experiment.get("dataset"))
        }
        run_start = time.time()
        self.align(alignment_parameters, sam_file_path)
        with open(runtime_log_path, "a") as runtime_log:
            runtime = str(datetime.timedelta(seconds=time.time() - run_start))
            runtime_log.write("Alignment: {}\n".format(runtime))

        # Create sorted BAM file from SAM file
        conversion_parameters = {
            "docker_client": docker_client,
            "docker_image": "gatk",
            "destination": destination
        }
        conversion_start = time.time()
        self.convert(conversion_parameters, sam_file_path, bam_file_path)
        with open(runtime_log_path, "a") as runtime_log:
            runtime = str(datetime.timedelta(seconds=time.time() - conversion_start))
            runtime_log.write("Convert to sorted BAM: {}\n".format(runtime))

    def build_genome_index(self, parameters):
        self.prepare_indexing(parameters)
        command = self.build_index_command(parameters)
        self.run_docker(command, parameters)

    def align(self, parameters, sam_file_path):
        command = self.alignment_command(parameters)
        output_parameters = {
            "log_is_output": not self.creates_output,
            "out_file_path": sam_file_path
        }
        self.run_docker(
            command,
            parameters,
            output_parameters
        )
        self.conclude_alignment(parameters, sam_file_path)

    def convert(self, parameters, sam_file_path, bam_file_path):
        """
        Convert SAM file to sorted BAM file. Also append read groups to ensure
        compatibility with Picard and GATK tools.
        """

        destination = parameters["destination"]
        command = "gatk AddOrReplaceReadGroups -I /{} -O /{} -SO coordinate " \
            "-ID foo -LB bar -PL illumina -SM Sample1 -PU foo.bar " \
            "--CREATE_INDEX true".format(sam_file_path, bam_file_path)
        output_parameters = { "log_file_path": destination + "Samtools.log" }
        self.run_docker(command, parameters, output_parameters)
        file_utils.validate_file_content(bam_file_path)

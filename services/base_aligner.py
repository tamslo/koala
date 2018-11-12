import os
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
                self.build_genome_index(index_parameters)
            except:
                file_utils.delete(temp_genome_index_path)
                raise
            os.rename(temp_genome_index_path, genome_index_path)

        # Run alignment
        alignment_parameters = {
            "docker_client": docker_client,
            "destination": destination,
            "genome_index_path": genome_index_path,
            "dataset": data_handler.datasets.select(experiment.get("dataset"))
        }
        self.align(alignment_parameters, sam_file_path)

        # Create sorted BAM file from SAM file
        post_processing_parameters = {
            "docker_client": docker_client,
            "docker_image": "gatk",
            "destination": destination,
            "reference_path": data_handler.reference_path(experiment)
        }
        self.post_process(post_processing_parameters, sam_file_path, bam_file_path)

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

    def post_process(self, parameters, sam_file_path, bam_file_path):
        """
        Convert SAM file to sorted BAM file. Also append read groups to ensure
        compatibility with Picard and GATK tools, mark duplicates, remove soft
        clips, and standardize the mapping quality.
        Command parameters from https://gatkforums.broadinstitute.org/gatk/discussion/3891/calling-variants-in-rnaseq
        """

        destination = parameters["destination"]

        sorted_path = destination + "Sorted.bam"
        command = "gatk AddOrReplaceReadGroups -I /{} -O /{} -SO coordinate " \
            "-ID foo -LB bar -PL illumina -SM Sample1 -PU foo.bar".format(
                sam_file_path,
                sorted_path
        )
        output_parameters = {
            "log_file_path": destination + "Conversion.log",
            "log_from_stderr": True
        }
        self.run_docker(command, parameters, output_parameters)
        file_utils.validate_file_content(sorted_path)

        deduplicated_path = destination + "Deduplicated.bam"
        metrics_path = destination + "Deduplicate.metrics"
        command = "gatk MarkDuplicates -I /{} -O /{} --CREATE_INDEX=true " \
            "--VALIDATION_STRINGENCY=SILENT -M /{}".format(
                sorted_path,
                deduplicated_path,
                metrics_path
            )
        output_parameters = {
            "log_file_path": destination + "Deduplicate.log",
            "log_from_stderr": True
        }
        self.run_docker(command, parameters, output_parameters)
        file_utils.validate_file_content(deduplicated_path)
        file_utils.delete(sorted_path)

        command = "gatk SplitNCigarReads -R /{} -I /{} -o /{} " \
            "-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 " \
            "-U ALLOW_N_CIGAR_READS".format(
                parameters["reference_path"],
                deduplicated_path,
                bam_file_path
            )
        output_parameters = { "log_file_path": destination + "SplitN.log" }
        self.run_docker(command, parameters, output_parameters)
        file_utils.validate_file_content(bam_file_path)
        file_utils.delete(deduplicated_path)

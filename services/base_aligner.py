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
        dataset = data_handler.datasets.select(experiment.get("dataset"))

        sam_file_path = destination + "Out.sam"
        bam_file_path = destination + "Out.bam"

        # Define genome index path and temp path (will be renamed if successful)
        parameters["reference_id"] = experiment.get("reference")
        genome_index_path = data_handler.genome_index_path(experiment, self.id)
        temp_genome_index_path = genome_index_path + ".running"

        # If neccessary, build genome index
        if not os.path.exists(genome_index_path):
            try:
                index_parameters = {
                    "docker_client": docker_client,
                    "destination": destination,
                    "genome_index_path": temp_genome_index_path,
                    "reference_path": data_handler.reference_path(experiment),
                    "dataset": dataset,
                    "reference_base_path": data_handler.reference_directory,
                    "reference_id": parameters["reference_id"]
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
            "dataset": dataset,
            "reference_id": parameters["reference_id"],
            "reference_base_path": data_handler.reference_directory
        }
        self.align(alignment_parameters, sam_file_path)

        # Create sorted BAM file from SAM file
        post_processing_parameters = {
            "docker_client": docker_client,
            "docker_image": "gatk",
            "destination": destination,
            "data_handler": data_handler,
            "experiment": experiment,
            "dataset": dataset
        }
        self.post_process(post_processing_parameters, sam_file_path, bam_file_path)

    def build_genome_index(self, parameters):
        self.prepare_indexing(parameters)
        command = self.build_index_command(parameters)
        self.run_docker(command, parameters, log_file_name="Index.log")

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
        file_utils.validate_file_content(sam_file_path)

    def post_process(self, parameters, sam_file_path, bam_file_path):
        destination = parameters["destination"]
        dataset = parameters["dataset"]

        # Convert to BAM, add read groups and sort
        command = "gatk AddOrReplaceReadGroups -I /{} -O /{} -SO coordinate " \
            "-ID foo -LB bar -PL illumina -SM Sample1 -PU foo.bar " \
            "--TMP_DIR {} " \
            "--CREATE_INDEX".format(
                sam_file_path,
                bam_file_path,
                destination
        )
        output_parameters = { "log_file_path": destination + "Conversion.log" }
        self.run_docker(command, parameters, output_parameters)
        file_utils.validate_file_content(bam_file_path)

        # Delete SAM file if not needed in evaluation (which is for BEERS sets)
        evaluation = dataset.get("evaluation")
        if evaluation == None or evaluation["type"] != "beers":
            file_utils.delete(sam_file_path)

        # Create reference indices
        data_handler = parameters["data_handler"]
        experiment = parameters["experiment"]
        reference_path = data_handler.reference_path(experiment)
        reference_index_path = data_handler.reference_path(
            experiment,
            alternate_file_ending=".fa.fai"
        )
        reference_dict_path = data_handler.reference_path(
            experiment,
            alternate_file_ending=".dict"
        )

        # Generate index of reference if not there
        if not os.path.exists(reference_index_path):
            command = "samtools faidx /{}".format(reference_path)
            output_parameters = { "log_file_path": destination + "Index.log" }
            self.run_docker(command, parameters, output_parameters)

        # Generate dict or reference if not there
        if not os.path.exists(reference_dict_path):
            command = "gatk CreateSequenceDictionary -R /{} -O /{}".format(
                reference_path,
                reference_dict_path
            )
            output_parameters = { "log_file_path": destination + "Dict.log" }
            self.run_docker(command, parameters, output_parameters)

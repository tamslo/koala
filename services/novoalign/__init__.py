import os
import yaml
from ..base_aligner import BaseAligner
import modules.file_utils as file_utils

class NovoAlign(BaseAligner):
    def __fasta_input(self, parameters):
        dataset = parameters["dataset"]
        evaluation = dataset.get("evaluation")
        return evaluation and evaluation["type"] == "beers"

    def __annotation_path(self, parameters):
        reference_id = parameters["reference_id"]
        annotation_base_path = parameters["annotation_base_path"]
        annotation_path = annotation_base_path + reference_id + ".gtf"
        if os.path.exists(annotation_path):
            return annotation_path

    def __annotated_index_path(self, parameters, file_ending=".idx.fa"):
        reference_base_path = parameters["reference_base_path"]
        reference_id = parameters["reference_id"]
        return reference_base_path + reference_id + file_ending

    def prepare_indexing(self, parameters):
        destination = parameters["destination"]
        reference_id = parameters["reference_id"]
        annotation_base_path = parameters["annotation_base_path"]
        reference_path = parameters["reference_path"]
        annotation_path = self.__annotation_path(parameters)
        if annotation_path and not os.path.exists(self.__annotated_index_path(parameters)):
            command = "rsem-prepare-reference --gtf {} {} {}".format(
                annotation_path,
                reference_path,
                self.__annotated_index_path(parameters, file_ending="")
            )
            self.run_docker(command, parameters, log_file_name="Annotation.log")

    def build_index_command(self, parameters):
        genome_index_path = parameters["genome_index_path"]
        reference_id = parameters["reference_id"]
        reference_path = parameters["reference_path"]
        if self.__annotation_path(parameters):
            reference_path = self.__annotated_index_path(parameters)
        command = "novoindex -n {} /{} /{}".format(
            reference_id,
            genome_index_path,
            reference_path
        )
        return command

    def alignment_command(self, parameters):
        # Testing
        command = "cp /data/errored/novoalign-20181213T223754/Out.sam.tmp " \
            "/{}Out.sam".format(parameters["destination"])
        return command

        # Actual code
        dataset = parameters["dataset"]
        genome_index_path = parameters["genome_index_path"]
        command = "novoalign -o SAM -f"
        file_utils.validate_file_content(genome_index_path)
        for direction, specification in dataset.get("data").items():
            command += " /{}".format(specification["path"])
        command += " -d /{}".format(genome_index_path)
        command += " -r All 10" # report max. 10 alignments per read
        command += " -v 0 70 70 '[>]([^:]*)'" # group junction and exon sequences together
        if self.__fasta_input(parameters):
            command += " -F FA"
        return command

    def conclude_alignment(self, parameters, out_file_path):
        # file_utils.validate_file_content(out_file_path)
        if self.__annotation_path(parameters):
            destination = parameters["destination"]
            with open("config.yml", "r") as config_file:
                config = yaml.load(config_file)
                num_threads = int(config["cores"])
            intermediate_sam_path = out_file_path + ".tmp"
            intermediate_bam_path = destination + "Tmp.bam"
            fixed_bam_path = destination + "Fixed.bam"
            os.rename(out_file_path, intermediate_sam_path)

            # SAM to BAM
            command = "samtools view -bS -f 2 /{}".format(intermediate_sam_path)
            output_parameters = {
                "log_is_output": True,
                "log_file_path": destination + "IntermediateConversion.log",
                "out_file_path": intermediate_bam_path
            }
            parameters["docker_image"] = "gatk"
            self.run_docker(command, parameters, output_parameters)
            parameters.pop("docker_image", None)
            file_utils.validate_file_content(intermediate_bam_path)

            # Fix coordinates
            command = "rsem-tbam2gbam /{} /{} /{} -p {}".format(
                self.__annotated_index_path(parameters, file_ending=""),
                intermediate_bam_path,
                fixed_bam_path,
                num_threads
            )
            self.run_docker(command, parameters, log_file_name="FixCoordinates.log")
            file_utils.validate_file_content(fixed_bam_path)

            # BAM to SAM again
            command = "samtools view -h /{}".format(fixed_bam_path)
            output_parameters = {
                "log_is_output": True,
                "log_file_path": destination + "IntermediateReconversion.log",
                "out_file_path": out_file_path
            }
            parameters["docker_image"] = "gatk"
            self.run_docker(command, parameters, output_parameters)
            parameters.pop("docker_image", None)

            file_utils.validate_file_content(out_file_path)
            file_utils.delete(intermediate_sam_path)
            file_utils.delete(intermediate_bam_path)
            file_utils.delete(fixed_bam_path)

class NovoAlignIndelSensitive(NovoAlign):
    def alignment_command(self, parameters):
        has_license = os.path.exists("services/novoalign/assets/novoalign.lic")
        if not has_license:
            raise Exception("NovoAlign in indel mode can only be run with license")

        command = super().alignment_command(parameters)
        command += " -x 3"
        command += " --matchreward 3"
        command += " --softclip 50,30" # four color, use 35,0 for two color
        return command

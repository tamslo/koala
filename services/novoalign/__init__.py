import os
import yaml
from ..base_aligner import BaseAligner
import modules.file_utils as file_utils

class NovoAlign(BaseAligner):
    def __fasta_input(self, parameters):
        dataset = parameters["dataset"]
        evaluation = dataset.get("evaluation")
        return evaluation and evaluation["type"] == "beers"

    def build_index_command(self, parameters):
        genome_index_path = parameters["genome_index_path"]
        reference_id = parameters["reference_id"]
        reference_path = parameters["reference_path"]
        command = "novoindex -n {} /{} /{}".format(
            reference_id,
            genome_index_path,
            reference_path
        )
        return command

    def alignment_command(self, parameters):
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

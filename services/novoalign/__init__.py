import os
from ..base_aligner import BaseAligner

class NovoAlign(BaseAligner):
    def build_index_command(self, parameters):
        return "./build_genome.sh {} {} {} {}".format(
            parameters["reference_base_path"],
            parameters["genome_index_path"],
            parameters["reference_id"],
            parameters["dataset"].get("readLength")
        )

    def alignment_command(self, parameters):
        dataset = parameters["dataset"]
        genome_index_path = parameters["genome_index_path"]
        command = "novoalign -o SAM -f"
        for direction, specification in dataset.get("data").items():
            command += " /{}".format(specification["path"])
        command += " -d /{}".format(genome_index_path)
        command += " -r All 10" # report max. 10 alignments per read
        return command

    def conclude_alignment(self, parameters, sam_file_path):
        command = "./fix_coordinates.sh {}".format(sam_file_path)
        output_parameters = {
            "log_file_path": parameters["destination"] + "Coordinates.log",
            "log_from_stderr": True
        }
        self.run_docker(command, parameters, output_parameters)

    def genome_index_amendment(self, dataset):
        return "_" + dataset["readLength"]


class NovoAlignIndelSensitive(NovoAlign):
    def alignment_command(self, parameters):
        has_license = os.path.exists("services/novoalign/assets/novoalign.lic")
        if not has_license:
            raise Exception("NovoAlign in indel mode can only be run with license")

        command = super().alignment_command(parameters)
        command += " -x 3"
        command += " --matchreward 3"
        command += " --softclip 50,30" # four color
        # command += " --softclip 35,0" # two color
        return command

class NovoAlignFasta(NovoAlign):
    def alignment_command(self, parameters):
        command = super().alignment_command(parameters)
        command += " -F FA"
        return command

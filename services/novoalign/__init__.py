import os
from ..base_aligner import BaseAligner

class NovoAlign(BaseAligner):
    def build_index_command(self, parameters):
        return "novoindex -n {} {} {}".format(
            parameters["reference_id"],
            parameters["genome_index_path"],
            parameters["reference_path"]
        )

    def alignment_command(self, parameters):
        dataset = parameters["dataset"]
        genome_index_path = parameters["genome_index_path"]
        command = "novoalign -o SAM -f"
        for direction, specification in dataset.get("data").items():
            command += " /{}".format(specification["path"])
        command += " -d /{}".format(genome_index_path)
        return command

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

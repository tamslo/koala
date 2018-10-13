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
        command = "novoalign -o SAM -f "
        for direction, specification in dataset.get("data").items():
            command += " /{}".format(specification["path"])
        command += " -d /{}".format(genome_index_path)
        return command

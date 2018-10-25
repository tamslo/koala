import yaml
from ..base_aligner import BaseAligner

class Star(BaseAligner):
    def __preamble(self, parameters):
        destination = parameters["destination"]
        genome_index_path = parameters["genome_index_path"]
        with open("config.yml", "r") as config_file:
            config = yaml.load(config_file)
            num_threads = int(config["cores"])
        return "STAR --runThreadN {} --genomeDir {} --outFileNamePrefix {}".format(
            num_threads,
            "/" + genome_index_path,
            destination
        )

    def build_index_command(self, parameters):
        command = self.__preamble(parameters) + " --runMode genomeGenerate"
        command += " --genomeFastaFiles /{}".format(parameters["reference_path"])
        return command

    def alignment_command(self, parameters):
        dataset = parameters["dataset"]
        command = self.__preamble(parameters) + " --readFilesIn"
        for direction, specification in dataset.get("data").items():
            command += " /{}".format(specification["path"])
        return command
import os, yaml
import modules.file_utils as file_utils
from ..base_aligner import BaseAligner

class Star(BaseAligner):
    def __preamble(self, parameters, destination_postfix=""):
        destination = parameters["destination"]
        genome_index_path = parameters["genome_index_path"]
        with open("config.yml", "r") as config_file:
            config = yaml.load(config_file)
            num_threads = int(config["cores"])
        return "STAR --runThreadN {} --genomeDir {} --outFileNamePrefix {}{}".format(
            num_threads,
            "/" + genome_index_path,
            destination,
            destination_postfix
        )

    def prepare_indexing(self, parameters):
        file_utils.create_directory(parameters["genome_index_path"])

    def conclude_alignment(self, parameters, sam_file_path):
        os.rename(
            os.path.join(parameters["destination"], "Aligned.out.sam"),
            sam_file_path
        )

    def build_index_command(self, parameters):
        command = self.__preamble(parameters, "Genome.Index.") + " --runMode genomeGenerate"
        command += " --genomeFastaFiles /{}".format(parameters["reference_path"])
        return command

    def alignment_command(self, parameters):
        dataset = parameters["dataset"]
        command = self.__preamble(parameters) + " --readFilesIn"
        for direction, specification in dataset.get("data").items():
            command += " /{}".format(specification["path"])
        command += " --outSAMmapqUnique 60" # default is 255 but GATK can't interpret this
        command += " --twopassMode Basic"
        command += " --twopass1readsN 1000000000"
        command += " --sjdbOverhang {}".format(int(dataset.get("readLength")) - 1)
        return command

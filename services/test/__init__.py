import os
import modules.file_utils as file_utils
from ..base_aligner import BaseAligner

class TestAlignerWritesLog(BaseAligner):
    def build_index_command(self, parameters):
        return "touch {}".format(parameters["genome_index_path"])

    def alignment_command(self, parameters):
        return "bash ./echo.sh"

class TestAlignerWritesFile(BaseAligner):
    def prepare_indexing(self, parameters):
        file_utils.create_directory(parameters["genome_index_path"])

    def conclude_alignment(self, parameters, sam_file_path):
        os.rename(
            os.path.join(parameters["destination"], "Aligned.out.sam"),
            sam_file_path
        )

    def build_index_command(self, parameters):
        return "touch {}/test".format(parameters["genome_index_path"])

    def alignment_command(self, parameters):
        return "cp /example.sam {}{}".format(parameters["destination"], "Aligned.out.sam")

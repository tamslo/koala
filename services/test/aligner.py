from ..base_aligner import BaseAligner

class TestAlignerWritesLog(BaseAligner):
    def build_index_command(self, parameters):
        return "touch {}".format(parameters["genome_index_path"])

    def alignment_command(self, parameters):
        return "bash ./echo.sh"

class TestAlignerWritesFile(BaseAligner):
    def build_index_command(self, parameters):
        return "touch {}".format(parameters["genome_index_path"])

    def alignment_command(self, parameters):
        print("touch {}/{}".format(parameters["destination"], self.output_file_name), )
        return "touch {}/{}".format(parameters["destination"], self.output_file_name)

from ..base_aligner import BaseAligner

class TestAligner(BaseAligner):
    def build_index_command(self, parameters):
        return "touch {}".format(genome_index_path)

    def alignment_command(self, parameters):
        return "bash ./echo.sh"

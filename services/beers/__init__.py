from ..base_service import BaseService

class BeersEvaluator(BaseService):
    def run(self, parameters):
        read_length = parameters["read_length"]
        alignment_directory = parameters["alignment_directory"]
        truth_file = parameters["truth_file"]
        command = "bash evaluate_alignment.sh {} {} {}".format(
            read_length,
            alignment_directory,
            truth_file
        )
        raise Exception("Needs to be tested")
        self.run_docker(parameters, command)

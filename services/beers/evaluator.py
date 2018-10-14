from ..base_service import BaseService

class BeersEvaluator(BaseService):
    def evaluate_alignment(self, parameters):
        docker_client = parameters["docker_client"]
        read_length = parameters["read_length"]
        alignment_directory = parameters["alignment_directory"]
        truth_file = parameters["truth_file"]
        docker_client.run(
            "beers",
            "bash evaluate_alignment.sh {} {} {}".format(
                read_length,
                alignment_directory,
                truth_file
            )
        )

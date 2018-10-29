from ..base_service import BaseService

class Opossum(BaseService):
    def run(self, parameters):
        # Tool does not throw an error. Test in BaseService whether the output
        # file is written by setting the out_file_name parameter.
        parameters["out_file_name"] = "Out.bam"

        command = "python Opossum.py --BamFile=/{} --OutFile=/{}".format(
            parameters["experiment"].get_input_directory(self.id) + "Out.bam",
            parameters["destination"] + parameters["out_file_name"]
        )
        self.run_docker(command, parameters)

from ..base_service import BaseService

class Opossum(BaseService):
    def command(self, parameters):
        return "python Opossum.py --BamFile=/{} --OutFile=/{}".format(
            parameters["experiment"].get_input_directory(self.id) + "Out.bam",
            parameters["destination"] + "Out.bam"
        )

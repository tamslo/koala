from ..base_service import BaseService

class Opossum(BaseService):
    def command(self, parameters):
        command = "python Opossum.py --BamFile={} --OutFile={}".format(
            experiment.get_input_directory(self.id) + "Out.bam",
            destination + "Out.bam"
        )
        raise Exception("TODO: Build Opossum command")

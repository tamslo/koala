import os, urllib
from modules.file_handler import FileHandler

class ReferenceGenomes:
    """
    Downloads specified tools and specified reference genomes.
    """

    fasta_file_ending = ".fa"
    rsync_uri = "rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/"
    tools = ["twoBitToFa"]
    started_tasks = []
    finished_tasks = []

    def __init__(self, data_directory):
        self.directory = data_directory + "references/"
        self.file_handler = FileHandler(self.directory)
        self.reference_genomes = {
            "hg19": self.__get_human_genome,
            "hg38": self.__get_human_genome,
            "pfal": self.__get_p_falciparum
        }
        try:
            self.__get_tools()
            self.__get_genomes()
        finally:
            self.clean_up()

    def __log_task_start(self, item, path):
        self.started_tasks.append(path)
        print("[ReferenceGenomes] Downloading {}...".format(item), flush=True)

    def __log_task_end(self, item, path):
        self.finished_tasks.append(path)
        print("[ReferenceGenomes] Downloaded {}".format(item), flush=True)

    def __get_tools(self):
        for tool_name in self.tools:
            tool_path = self.directory + tool_name
            if not os.path.exists(tool_path):
                self.__log_task_start(tool_name, tool_path)
                tool_uri = self.rsync_uri + tool_name
                os.system("rsync -aPq {} {}".format(tool_uri, tool_path))
                self.__log_task_end(tool_name, tool_path)

    def __genome_path(self, genome_name):
        return self.directory + genome_name + self.fasta_file_ending

    def __get_genomes(self):
        for genome_name, genome_getter in self.reference_genomes.items():
            genome_path = self.__genome_path(genome_name)
            if not os.path.exists(genome_path):
                self.__log_task_start(genome_name, genome_path)
                genome_getter(genome_name, genome_path)
                self.__log_task_end(genome_name, genome_path)

    def __get_human_genome(self, genome_name, genome_path):
        url = "http://hgdownload.soe.ucsc.edu/goldenPath/"
        url += "{0}/bigZips/{0}.2bit".format(genome_name)
        two_bit_path = genome_path + ".2bit"
        self.started_tasks.append(two_bit_path)
        self.file_handler.download(url, two_bit_path)
        self.finished_tasks.append(two_bit_path)
        # Convert .2bit file to .fa
        # TODO log start and end of tasks
        # os.system("cd {} && twoBitToFa {} {}".format(
        #     self.directory,
        #     two_bit_path,
        #     genome_path
        # ))
        # self.file_handler.delete(two_bit_path)

    def __get_p_falciparum(self, genome_name, genome_path):
        # url = http://bp1.s3.amazonaws.com/malaria.tar.bz2
        # download_path = self.directory + "malaria.tar.bz2"
        # self.file_handler.download(url, )
        # Extract right file from tarball
        # file_name = genome_sequence_pfal.fa
        # TODO extract
        return None

    def select(self):
        return None

    def clean_up(self):
        for path in self.started_tasks:
            if not path in self.finished_tasks:
                self.file_handler.delete(path)

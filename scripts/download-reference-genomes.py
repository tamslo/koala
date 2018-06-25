import sys, os

# Make import from parent directory possible
sys.path.append(
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import modules.file_utils as file_utils

def get_human_genome(genome_name, file_path):
    url = "http://hgdownload.soe.ucsc.edu/goldenPath/"
    url += "{0}/bigZips/{0}.2bit".format(genome_name)
    two_bit_path = file_path + ".2bit"
    started_tasks.append(two_bit_path)
    file_utils.download(url, two_bit_path)
    finished_tasks.append(two_bit_path)
    # Convert .2bit file to .fa
    os.system("touch {0} && chmod +x {0}".format(file_path))
    os.system("cd {} && chmod +x twoBitToFa && ./twoBitToFa {} {}".format(
        directory,
        two_bit_path,
        file_path
    ))
    # file_utils.delete(two_bit_path)

def get_p_falciparum(genome_name, file_path):
    # url = http://bp1.s3.amazonaws.com/malaria.tar.bz2
    # download_path = directory + "malaria.tar.bz2"
    # file_utils.download(url, )
    # Extract right file from tarball
    # file_name = genome_sequence_pfal.fa
    # TODO extract
    return None

# Add new reference genomes with options here
reference_genomes = {
    "hg19": {"getter": get_human_genome, "name": "Human (hg19)"},
    # "hg38": {"getter": get_human_genome, "name": "Human (hg38)"},
    # "pfal": {"getter": get_p_falciparum, "name": "Malaria"}
}

# Add new tools needed to download reference genomes here
tools = ["twoBitToFa"]

# Constants
directory = "data/references/"
fasta_file_ending = ".fa"
rsync_uri = "rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/"

started_tasks = []
finished_tasks = []

def log_task_start(item, path):
    started_tasks.append(path)
    print("Downloading {}...".format(item), flush=True)

def log_task_end(item, path):
    finished_tasks.append(path)
    print("Downloaded {}".format(item), flush=True)

def get_tools():
    for tool_name in tools:
        tool_path = directory + tool_name
        if not os.path.exists(tool_path):
            log_task_start(tool_name, tool_path)
            tool_uri = rsync_uri + tool_name
            os.system("rsync -aPq {} {}".format(tool_uri, tool_path))
            log_task_end(tool_name, tool_path)

def genome_path(genome_name):
    return directory + genome_name + fasta_file_ending

def get_genomes():
    for genome_name, genome_options in reference_genomes.items():
        file_path = genome_path(genome_name)
        if not os.path.exists(file_path):
            log_task_start(genome_name, file_path)
            genome_options["getter"](genome_name, file_path)
            log_task_end(genome_name, file_path)

print("Downloading reference genomes", flush=True)
file_utils.create_directory(directory)
try:
    get_tools()
    get_genomes()
finally:
    for path in started_tasks:
        if not path in finished_tasks:
            file_utils.delete(path)

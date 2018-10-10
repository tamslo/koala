import sys, os, shutil, json

# Make import from parent directory possible
sys.path.append(
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import modules.file_utils as file_utils

# Directory Paths
reference_directory = "data/references/"
datasets_directory = "data/datasets/"
dataset_jsons_directory = "scripts/prepared_datasets/"

def log_task_start(item, path):
    started_tasks.append(path)
    print("Downloading {}...".format(item), flush=True)

def log_task_end(item, path):
    finished_tasks.append(path)
    print("Downloaded {}".format(item), flush=True)

def log_data_present(item):
    print("{} already present".format(item), flush=True)

####################
# REFERENCE GENOMES
####################

# Add new tools needed to download reference genomes here
tools = ["twoBitToFa"]

# Constants
fasta_file_ending = ".fa"
rsync_uri = "rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/"

started_tasks = []
finished_tasks = []

def get_human_genome(genome_id, file_path):
    url = "http://hgdownload.soe.ucsc.edu/goldenPath/"
    url += "{0}/bigZips/{0}.2bit".format(genome_id)
    two_bit_path = file_path + ".2bit"
    started_tasks.append(two_bit_path)
    file_utils.download(url, two_bit_path)
    finished_tasks.append(two_bit_path)
    # Convert .2bit file to .fa
    print("Extracting {} from 2bit file...".format(genome_id), flush=True)
    os.system("chmod +x {0}twoBitToFa && {0}twoBitToFa {1} {2}".format(
        reference_directory,
        two_bit_path,
        file_path
    ))
    file_utils.delete(two_bit_path)

def get_p_falciparum(genome_id, file_path):
    url = "http://bp1.s3.amazonaws.com/malaria.tar.bz2"
    download_path = reference_directory + "malaria.tar.bz2"
    file_utils.download(url, download_path)
    print("Unzipping {}...".format(genome_id), flush=True)
    unzipped_directory = file_utils.unzip(download_path)
    os.rename(unzipped_directory + "/genome_sequence_pfal.fa", file_path)
    file_utils.delete(download_path)
    file_utils.delete(unzipped_directory)

# Add new reference genomes with options here
genomes = {
    "hg19": {
        "getter": get_human_genome,
        "name": "Human (hg19)",
        "source": "http://hgdownload.cse.ucsc.edu/downloads.html#human"
    },
    "hg38": {
        "getter": get_human_genome,
        "name": "Human (hg38)",
        "source": "http://hgdownload.cse.ucsc.edu/downloads.html#human"
    },
    "pfal": {
        "getter": get_p_falciparum,
        "name": "Malaria",
        "source": "http://bioinf.itmat.upenn.edu/BEERS/bp1/datasets.php"
    }
}

def get_tools():
    for tool_name in tools:
        tool_path = reference_directory + tool_name
        if not os.path.exists(tool_path):
            log_task_start(tool_name, tool_path)
            tool_uri = rsync_uri + tool_name
            os.system("rsync -aPq {} {}".format(tool_uri, tool_path))
            log_task_end(tool_name, tool_path)
        else:
            log_data_present(tool_name)

def genome_path(genome_id):
    return reference_directory + genome_id + fasta_file_ending

def get_genomes():
    genome_infos_path = os.path.join(reference_directory, "references.json")
    genome_infos = []
    if os.path.exists(genome_infos_path):
        with open(genome_infos_path, "r") as genome_infos_file:
            genome_infos = json.load(genome_infos_file)

    for genome_id, genome_specification in genomes.items():
        file_path = genome_path(genome_id)
        info_path = file_path.split(fasta_file_ending)[0] + ".yml"
        genome_getter = genome_specification["getter"]
        if not os.path.exists(file_path):
            log_task_start(genome_id, file_path)
            genome_getter(genome_id, file_path)
            genome_infos.append({
                "id": genome_id,
                "name": genome_specification["name"],
                "source": genome_specification["source"]
            })
            with open(genome_infos_path, "w") as genome_infos_file:
                genome_infos_file.write(json.dumps(genome_infos))
            log_task_end(genome_id, file_path)
        else:
            log_data_present(genome_id)

###################
# RNASEQ DATA SETS
###################

def get_baruzzo(dataset, directory):
    zip_name = "{}.tar.bz2".format(dataset["name"])
    url = "http://bp1.s3.amazonaws.com/{}".format(zip_name)
    download_path = directory + "/" + zip_name
    file_utils.download(url, download_path)

    print("Unzipping {}...".format(dataset["name"]), flush=True)
    file_utils.unzip(download_path)

    # Move files to /beers directory
    beers_directory = directory + "/beers/"
    file_utils.create_directory(beers_directory)
    for file_name in os.listdir(directory):
        file_path = directory + "/" + file_name
        if not os.path.isdir(file_path) and not file_path == download_path:
            shutil.move(file_path, beers_directory + file_name)

    # Move FASTA files to root and rename
    dataset_info_path = dataset_jsons_directory + dataset["id"] + ".json"
    with open(dataset_info_path, "r") as dataset_info_file:
        dataset_info = json.load(dataset_info_file)
    for direction, file_specification in dataset_info["data"].items():
        downloaded_file_path = beers_directory + file_specification["name"]
        os.rename(downloaded_file_path, file_specification["path"])

    file_utils.delete(download_path)

# IDs must match IDs in prepared_datasets
rna_seq_data = [
    {
        "id": "simulated_reads_HG19t1r1",
        "getter": get_baruzzo,
        "name": "human_t1r1"
    },
    {
        "id": "simulated_reads_HG19t2r1",
        "getter": get_baruzzo,
        "name": "human_t2r1"
    },
    {
        "id": "simulated_reads_PFALt1r1",
        "getter": get_baruzzo,
        "name": "malaria_t1r1"
    }
]


rna_seq_getters = {
    "baruzzo": get_baruzzo
}

def get_datasets():
    for dataset in rna_seq_data:
        dataset_directory = datasets_directory + dataset["id"]
        dataset_json_origin = dataset_jsons_directory + dataset["id"] + ".json"
        dataset_json_destination = datasets_directory + dataset["id"] + ".json"
        dataset_getter = dataset["getter"]
        if not os.path.isdir(dataset_directory):
            file_utils.create_directory(dataset_directory)
            log_task_start(dataset["name"], dataset_directory)
            dataset_getter(dataset, dataset_directory)
            shutil.copyfile(dataset_json_origin, dataset_json_destination)
            log_task_end(dataset["name"], dataset_directory)
        else:
            log_data_present(dataset["name"])

###################
# SCRIPT EXECUTION
###################

print("", flush=True)
print("Downloading data", flush=True)
print("", flush=True)

file_utils.create_directory(reference_directory)
file_utils.create_directory(datasets_directory)

try:
    get_tools()
    get_genomes()
    get_datasets()
finally:
    for path in started_tasks:
        if not path in finished_tasks:
            print("An error occured, deleting {}".format(path))
            file_utils.delete(path)

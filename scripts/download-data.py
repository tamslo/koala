#!/usr/bin/env python3

import sys, os, shutil, json, yaml
from time import localtime

# Make import from parent directory possible
sys.path.append(
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import modules.file_utils as file_utils

with open("constants.yml", "r") as constants_file:
    constants = yaml.load(constants_file)

# Directory Paths
reference_directory = "data/references/"
datasets_directory = "data/datasets/"

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
fastq_file_ending = ".fastq"
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
    # "hg38": {
    #     "getter": get_human_genome,
    #     "name": "Human (hg38)",
    #     "source": "http://hgdownload.cse.ucsc.edu/downloads.html#human"
    # },
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

def remove_tools():
    for tool_name in tools:
        tool_path = reference_directory + tool_name
        file_utils.delete(tool_path)

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

def write_dataset_json(info):
    dataset_info_path = datasets_directory + info["id"] + ".json"
    info["method"] = constants["dataset"]["FILE"]
    info["layout"] = constants["dataset"]["PAIRED"]
    info["created"] = localtime()
    info["error"] = False
    with open(dataset_info_path, "w") as dataset_info_file:
        json.dump(info, dataset_info_file)

def get_baruzzo(dataset, directory):
    zip_name = "{}.tar.bz2".format(dataset["file_name"])
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

    # Move FASTQ files to root and rename
    def setup_file(direction):
        file_name = "{}.{}.fa".format(dataset["id"], direction)
        file_origin = beers_directory + file_name
        file_destination = "{}/{}{}".format(directory, direction, fastq_file_ending)
        os.rename(file_origin, file_destination)
        return file_name, file_destination

    forward_file_name, forward_file_path = setup_file(constants["dataset"]["FORWARD"])
    reverse_file_name, reverse_file_path = setup_file(constants["dataset"]["REVERSE"])

    # Move CIG file to root and rename
    truth_file_name = "{}.cig".format(dataset["id"])
    truth_file_path = directory + "/truth.cig"
    os.rename(beers_directory + truth_file_name, truth_file_path)

    file_utils.delete(download_path)
    file_utils.delete(beers_directory)

    write_dataset_json({
        "id": dataset["id"],
        "name": dataset["name"],
        "readLength": "100",
        "data": {
            constants["dataset"]["FORWARD"]: {
                "name": forward_file_name,
                "path": forward_file_path,
            },
            constants["dataset"]["REVERSE"]: {
                "name": reverse_file_name,
                "path": reverse_file_path,
            }
        },
        "evaluation": {
            "type": "beers",
            "truth_file": {
                "name": truth_file_name,
                "path": truth_file_path
            }
        }
    })

def get_from_encode(dataset, directory):
    dataset_info = {
        "id": dataset["id"],
        "name": dataset["name"],
        "readLength": "76",
        "data": {
            constants["dataset"]["FORWARD"]: {},
            constants["dataset"]["REVERSE"]: {}
        },
        "evaluation": dataset["evaluation"]
    }

    def get_file(file_id, direction, directory):
        print("Downloading {} file...".format(direction), flush=True)
        zip_name = "{}.fastq.gz".format(file_id)
        url = "https://www.encodeproject.org/files/{}/@@download/{}".format(
            file_id,
            zip_name
        )
        download_path = directory + "/" + zip_name
        file_utils.download(url, download_path)
        print("Unzipping {} file...".format(direction), flush=True)
        file_utils.unzip(download_path)
        file_utils.delete(download_path)

        original_name = "{}.fastq".format(file_id)
        file_origin = "{}/{}".format(directory, original_name)
        file_destination = "{}/{}{}".format(directory, direction, fastq_file_ending)
        os.rename(file_origin, file_destination)
        return original_name, file_destination

    for direction, file_id in dataset["files"].items():
        original_name, file_destination = get_file(file_id, direction, directory)
        dataset_info["data"][direction]["name"] = original_name
        dataset_info["data"][direction]["path"] = file_destination

    write_dataset_json(dataset_info)

# Baruzzo Data Sets
# * id is prefix of unzipped FASTA files
# * file_name is zip name given in download url
rna_seq_data = [
    # {
    #     "id": "GM12878",
    #     "name": "GIAB Pilot Genome",
    #     "getter": get_from_encode,
    #     "files": {
    #         constants["dataset"]["FORWARD"]: "ENCFF000EWJ",
    #         constants["dataset"]["REVERSE"]: "ENCFF000EWX"
    #     },
    #     "evaluation": { "type": "giab" }
    # },
    {
        "id": "simulated_reads_HG19t1r1",
        "getter": get_baruzzo,
        "file_name": "human_t1r1",
        "name": "Simulated Human T1R1"
    },
    # {
    #     "id": "simulated_reads_HG19t1r2",
    #     "getter": get_baruzzo,
    #     "file_name": "human_t1r2",
    #     "name": "Simulated Human T1R2"
    # },
    # {
    #     "id": "simulated_reads_HG19t1r3",
    #     "getter": get_baruzzo,
    #     "file_name": "human_t1r3",
    #     "name": "Simulated Human T1R3"
    # },
    {
        "id": "simulated_reads_HG19t2r1",
        "getter": get_baruzzo,
        "file_name": "human_t2r1",
        "name": "Simulated Human T2R1"
    },
    # {
    #     "id": "simulated_reads_HG19t2r2",
    #     "getter": get_baruzzo,
    #     "file_name": "human_t2r2",
    #     "name": "Simulated Human T2R2"
    # },
    # {
    #     "id": "simulated_reads_HG19t2r3",
    #     "getter": get_baruzzo,
    #     "file_name": "human_t2r3",
    #     "name": "Simulated Human T2R3"
    # },
    {
        "id": "simulated_reads_HG19t3r1",
        "getter": get_baruzzo,
        "file_name": "human_t3r1",
        "name": "Simulated Human T3R1"
    },
    # {
    #     "id": "simulated_reads_HG19t3r2",
    #     "getter": get_baruzzo,
    #     "file_name": "human_t3r2",
    #     "name": "Simulated Human T3R2"
    # },
    # {
    #     "id": "simulated_reads_HG19t3r3",
    #     "getter": get_baruzzo,
    #     "file_name": "human_t3r3",
    #     "name": "Simulated Human T3R3"
    # },
    {
        "id": "simulated_reads_PFALt1r1",
        "getter": get_baruzzo,
        "file_name": "malaria_t1r1",
        "name": "Simulated Malaria T1R1"
    },
    # {
    #     "id": "simulated_reads_PFALt1r2",
    #     "getter": get_baruzzo,
    #     "file_name": "malaria_t1r2",
    #     "name": "Simulated Malaria T1R2"
    # },
    # {
    #     "id": "simulated_reads_PFALt1r3",
    #     "getter": get_baruzzo,
    #     "file_name": "malaria_t1r3",
    #     "name": "Simulated Malaria T1R3"
    # },
    {
        "id": "simulated_reads_PFALt2r1",
        "getter": get_baruzzo,
        "file_name": "malaria_t2r1",
        "name": "Simulated Malaria T2R1"
    },
    # {
    #     "id": "simulated_reads_PFALt2r2",
    #     "getter": get_baruzzo,
    #     "file_name": "malaria_t2r2",
    #     "name": "Simulated Malaria T2R2"
    # },
    # {
    #     "id": "simulated_reads_PFALt2r3",
    #     "getter": get_baruzzo,
    #     "file_name": "malaria_t2r3",
    #     "name": "Simulated Malaria T2R3"
    # },
    {
        "id": "simulated_reads_PFALt3r1",
        "getter": get_baruzzo,
        "file_name": "malaria_t3r1",
        "name": "Simulated Malaria T3R1"
    },
    # {
    #     "id": "simulated_reads_PFALt3r2",
    #     "getter": get_baruzzo,
    #     "file_name": "malaria_t3r2",
    #     "name": "Simulated Malaria T3R2"
    # },
    # {
    #     "id": "simulated_reads_PFALt3r3",
    #     "getter": get_baruzzo,
    #     "file_name": "malaria_t3r3",
    #     "name": "Simulated Malaria T3R3"
    # }
]

def get_datasets():
    for dataset in rna_seq_data:
        dataset_directory = datasets_directory + dataset["id"]
        dataset_getter = dataset["getter"]
        if not os.path.isdir(dataset_directory):
            file_utils.create_directory(dataset_directory)
            log_task_start(dataset["name"], dataset_directory)
            dataset_getter(dataset, dataset_directory)
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
    remove_tools()
finally:
    for path in started_tasks:
        if not path in finished_tasks:
            print("An error occured, deleting {}".format(path))
            file_utils.delete(path)

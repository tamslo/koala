import os
import json
import yaml
import sys
from time import localtime
import shutil

INTERMEDIATE_DIRECTORY = "data/intermediate"
DATASET_DIRECTORY = "data/datasets"

def get_read_length():
    first_data_set = os.listdir(INTERMEDIATE_DIRECTORY)[0]
    first_forward_file = os.path.join(INTERMEDIATE_DIRECTORY, first_data_set, "{}_1.fq".format(first_data_set))
    with open(first_forward_file) as file:
        # Skip first read which seems smaller than succeeding reads
        for skip_line in [file.readline] * 4:
            skip_line()

        read_length = None
        while read_length == None:
            line = file.readline()
            if len(line) < 10 or line.startswith("@"):
                continue
            else:
                read_length = len(line)

    return read_length

READ_LENGTH = get_read_length()

with open("constants.yml", "r") as constants_file:
    constants = yaml.load(constants_file)

def move_files(intermediate_folder, dataset_folder, dataset_id):
    intermediate_forward_name = "{}_1.fq".format(dataset_id)
    intermediate_forward_path = os.path.join(intermediate_folder, intermediate_forward_name)
    forward_file_path = os.path.join(dataset_folder, "forward.fastq")

    intermediate_reverse_name = "{}_2.fq".format(dataset_id)
    intermediate_reverse_path = os.path.join(intermediate_folder, intermediate_reverse_name)
    reverse_file_path = os.path.join(dataset_folder, "reverse.fastq")

    shutil.move(intermediate_forward_path, forward_file_path)
    shutil.move(intermediate_reverse_path, reverse_file_path)

    data_info = {
        constants["dataset"]["FORWARD"]: {
            "name": intermediate_forward_name,
            "path": forward_file_path,
        },
        constants["dataset"]["REVERSE"]: {
            "name": intermediate_reverse_name,
            "path": reverse_file_path,
        }
    }

    return data_info

def create_dataset_info(dataset_id, dataset_folder, data_info):
    dataset_info_path = dataset_folder + ".json"
    dataset_info = {
        "id": "smart_{}".format(dataset_id),
        "name": "SMART {}".format(dataset_id),
        "method": constants["dataset"]["FILE"],
        "layout": constants["dataset"]["PAIRED"],
        "created": localtime(),
        "readLength": str(READ_LENGTH),
        "data": data_info,
        "error": False
    }

    with open(dataset_info_path, "w") as dataset_info_file:
        dataset_info_file.write(json.dumps(dataset_info))

def create_dataset(dataset_id, intermediate_folder, dataset_folder):
    os.mkdir(dataset_folder)
    data_info = move_files(intermediate_folder, dataset_folder, dataset_id)
    create_dataset_info(
        dataset_id,
        dataset_folder,
        data_info
    )

for intermediate_path in os.listdir(INTERMEDIATE_DIRECTORY):
    dataset_id = intermediate_path
    intermediate_folder = os.path.join(INTERMEDIATE_DIRECTORY, dataset_id)
    dataset_folder = os.path.join(DATASET_DIRECTORY, "smart_{}".format(dataset_id))

    if not os.path.isdir(dataset_folder):
        print("Creating data set from {}".format(dataset_id))
        create_dataset(dataset_id, intermediate_folder, dataset_folder)

    shutil.rmtree(intermediate_path)

shutil.rmtree(INTERMEDIATE_DIRECTORY)

#!/usr/bin/python

import csv
import os
import sys

from sam import process_sam_line
from vcf import process_vcf_line

description = "\nCountInsertionsAndDeletions - Counting insertions and " \
    "deletions in position sorted SAM or VCF files. Takes paths to one or " \
    "multiple files to be counted and optionally an output directory as last " \
    "parameter.\n"

def get_parameters():
    input_file_paths = sys.argv[1:]
    if len(input_file_paths) < 1:
        print("[ERROR] The path to one or multiple files to be" \
            " processed is expected in input_file_paths")
        sys.exit(1)

    output_directory = ""
    if os.path.isdir(input_file_paths[-1]):
        output_directory = input_file_paths[-1]
        input_file_paths = input_file_paths[:-1]
        print("[INFO] Treating last parameter '{}' as output directory".format(
            output_directory
        ))
    return input_file_paths, output_directory

def process_file(path, index, output_directory):
    line_processors = {
        "sam": process_sam_line,
        "vcf": process_vcf_line
    }
    file_ending = path.split(".")[-1]

    if not file_ending in line_processors:
        print("[WARNING] Cannot process {}, unknown file format".format(path))
        return None

    print("[INFO] Processing {}...".format(path))
    intermediate_result_path = os.path.join(
        output_directory,
        "intermediate_result_{}.csv".format(str(index))
    )
    process_line = line_processors[file_ending]
    input_file = open(path, "r")
    output_file = open(intermediate_result_path, "w")
    csv_writer = csv.writer(output_file)
    header = ["file", "chr", "pos", "type", "len"]
    csv_writer.writerow(header)
    for line in input_file:
        line_result = process_line(line)
        if line_result != None:
            for row in line_result:
                csv_writer.writerow([path] + row)
    input_file.close()
    output_file.close()
    print("[INFO] Done.")
    return intermediate_result_path

def main():
    print(description)
    input_file_paths, output_directory = get_parameters()
    output_file_path = os.path.join(output_directory, "IndelCounts.csv")

    if os.path.exists(output_file_path):
        print("[WARNING] File {} already exists, it will be overwritten if " \
            "execution is not interrupted before file processing is " \
            "finished".format(output_file_path))

    intermediate_results = []
    for index, path in enumerate(input_file_paths):
        intermediate_result_path = process_file(path, index, output_directory)
        if intermediate_result_path != None:
            intermediate_results.append(intermediate_result_path)

    # TODO: Accumulate results from intermediate_results
    # print("[INFO] Writing CSV file...")
    # utils.overwrite_file(output_file_path)
    # rows = []
    # header = ["file", "chromosome"]
    # for count_name in empty_line_result()["counts"]:
    #     header.append(count_name)
    # rows.append(header)
    # for path, chromosomes in file_results.items():
    #     for chromosome, counts in chromosomes.items():
    #         values = [path, chromosome]
    #         for count_key, count_value in counts.items():
    #             values.append(count_value)
    #         rows.append(values)
    # utils.write_csv(output_file_path, rows)
    # print("[INFO] Done.")

main()

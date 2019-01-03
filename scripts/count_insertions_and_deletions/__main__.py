#!/usr/bin/python

import csv
import os
import sys

import utils
from sam import process_sam_line
from vcf import process_vcf_line

description = "\nCountInsertionsAndDeletions - Counting insertions and " \
    "deletions in position sorted SAM or VCF files. Takes paths to one or " \
    "multiple files to be counted and optionally an output directory as last " \
    "parameter.\n"

intermediate_result_header = ["file", "chr", "pos", "type", "len"]
accumulated_result_header = ["file", "chr", "insertions", "insertion_bases",
    "duplicate_insertions", "deletions", "deletion_bases", "duplicate_deletions"]

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

def write_to_csv(file, row):
    csv_writer = csv.writer(file)
    csv_writer.writerow(row)

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
    write_to_csv(output_file, intermediate_result_header)
    for line in input_file:
        line_result = process_line(line)
        if line_result != None:
            for row in line_result:
                write_to_csv(output_file, [path] + row)
    input_file.close()
    output_file.close()
    print("[INFO] Done.")
    return intermediate_result_path

def has_duplicate(start_position, length, variants):
    if start_position in variants:
        return start_position
    if len(variants) == 0:
        return None

    # Check if current positions are in range of present positions
    sorted_present_positions = sorted(variants)
    end_position = start_position + length
    if end_position < sorted_present_positions[0]:
        return None
    last_present_position = sorted_present_positions[-1]
    last_present_length = max(variants[last_present_position])
    if start_position > last_present_position + last_present_length:
        return None

    # Check if current positions are included in present positions
    duplicate_position = None
    for present_start_position in sorted_present_positions:
        if end_position < present_start_position:
            break # we have found no duplicate
        present_length = variants[present_start_position][0]
        present_end_position = present_start_position + present_length
        if (start_position > present_start_position and
            start_position < present_end_position):
            duplicate_position = present_start_position
            break # duplicate found
    return duplicate_position

def collect_intermediate_results(intermediate_result_path):
    file = None
    intermediate_results = {}
    intermediate_result_file = open(intermediate_result_path)
    csv_reader = csv.reader(intermediate_result_file)

    # Collect results into dictionary with chromosome, type, and position
    # as key to detect duplicates
    for row in csv_reader:
        if row != intermediate_result_header:
            file = row[0]
            chromosome = row[1]
            position = int(row[2])
            type = row[3]
            length = int(row[4])

            # Ensure keys are present
            if not chromosome in intermediate_results:
                intermediate_results[chromosome] = {}
            if not type in intermediate_results[chromosome]:
                intermediate_results[chromosome][type] = {}

            # Initialize lenghts array or append if duplicate or already present
            duplicate_position = None
            if file.endswith("sam"):
                duplicate_position = has_duplicate(
                    position,
                    length,
                    intermediate_results[chromosome][type]
                )
            if duplicate_position != None:
                intermediate_results[chromosome][type][duplicate_position].append(length)
            else:
                intermediate_results[chromosome][type][position] = [length]
    intermediate_result_file.close()
    return file, intermediate_results

def get_counts(variation_results):
    variation_count = 0
    base_count = 0
    duplicate_count = 0
    for position, lengths in variation_results.items():
        variation_count += len(lengths)
        base_count += sum(lengths)
        duplicate_count += len(lengths) - 1
    return variation_count, base_count, duplicate_count

def accumulate_intermediate_results(chromosome_results):
    insertions = 0
    insertion_bases = 0
    insertion_duplicates = 0
    deletions = 0
    deletion_bases = 0
    deletion_duplicates = 0

    if utils.insertion_type in chromosome_results:
        insertions, insertion_bases, insertion_duplicates = get_counts(
            chromosome_results[utils.insertion_type])

    if utils.deletion_type in chromosome_results:
        deletions, deletion_bases, deletion_duplicates = get_counts(
            chromosome_results[utils.deletion_type])

    return [insertions, insertion_bases, insertion_duplicates, deletions,
        deletion_bases, deletion_duplicates]

def acculumate_results(intermediate_results, output_file):
    write_to_csv(output_file, accumulated_result_header)
    for intermediate_result_path in intermediate_results:
        file, intermediate_results = collect_intermediate_results(intermediate_result_path)
        for chromosome, chromosome_results in intermediate_results.items():
            accumulated_result = accumulate_intermediate_results(chromosome_results)
            write_to_csv(
                output_file,
                [file, chromosome] + accumulated_result
            )

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

    print("[INFO] Accumulating output file...")
    output_file = open(output_file_path, "w")
    acculumate_results(intermediate_results, output_file)
    output_file.close()
    print("[INFO] Done.")
main()

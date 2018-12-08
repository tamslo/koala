#!/usr/bin/python

# Script to collect counts of insertions and deletions from SAM or VCF files,
# takes a list of paths to files to be processed and outputs a CSV file  called
# IndelCounts.csv.
# If the last path is a dictionary, the IndelCounts.csv is saved to this
# folder.

import os
import sys
import re
import csv

def raise_error(error_message):
    print(error_message)
    sys.exit(1)

def empty_line_result():
    return {
        "chromosome": None,
        "counts": {
            "insertions": 0,
            "inserted_bases": 0,
            "deletions": 0,
            "deleted_bases": 0
        }
    }

def number_of_indel_bases(cigar_string, index):
    count = 1 # default, if no number is given
    if index > 0:
        count_index = index - 1
        try:
            count = int(cigar_string[count_index])
        except ValueError:
            pass
    return count

def process_sam_line(line):
    line_result = empty_line_result()
    if not line.startswith("@"):
        sam_fields = fields = line.split("\t")
        if len(sam_fields) < 11:
            return -1
        chromosome = sam_fields[2]
        line_result["chromosome"] = chromosome
        cigar_string = sam_fields[5]
        # Process insertions
        insertion_indices = [match.start() for match in re.finditer("I", cigar_string)]
        if len(insertion_indices) > 0:
            line_result["counts"]["insertions"] = len(insertion_indices)
            for index in insertion_indices:
                line_result["counts"]["inserted_bases"] += number_of_indel_bases(cigar_string, index)
        # Process deletions
        deletion_indices = [match.start() for match in re.finditer("D", cigar_string)]
        if len(deletion_indices) > 0:
            line_result["counts"]["deletions"] = len(deletion_indices)
            for index in deletion_indices:
                line_result["counts"]["deleted_bases"] += number_of_indel_bases(cigar_string, index)

    return line_result

def process_vcf_line(line):
    line_result = empty_line_result()
    if not line.startswith("#"):
        vcf_fields = fields = line.split("\t")
        if len(vcf_fields) < 10:
            return -1
        chromosome = vcf_fields[0]
        line_result["chromosome"] = chromosome
        reference_bases = vcf_fields[3]
        alternative_bases = vcf_fields[4]
        # Process insertions
        if len(alternative_bases) > len(reference_bases):
            line_result["counts"]["insertions"] = 1
            line_result["counts"]["inserted_bases"] = len(alternative_bases) - len(reference_bases)
        # Process deletions
        if len(reference_bases) > len(alternative_bases):
            line_result["counts"]["deletions"] = 1
            line_result["counts"]["deleted_bases"] = len(reference_bases) - len(alternative_bases)
    return line_result

def increase_values(values, line_values):
    chromosome = line_values["chromosome"]
    if not chromosome in values:
        values[chromosome] = {}
    for count_key, count_value in line_values["counts"].items():
        if count_key in values[chromosome]:
            values[chromosome][count_key] += count_value
        else:
            values[chromosome][count_key] = count_value
    return values


def main():
    parameters = sys.argv[1:]
    if len(parameters) < 1:
        raise_error("[NO PATH] The path to one or multiple files to be" \
            " processed is expected as parameters")

    output_directory = ""
    if os.path.isdir(parameters[-1]):
        print("[INFO] Treating last path as output directory")
        output_directory = parameters[-1]
        parameters = parameters[:-1]
    output_file_path = output_directory + "IndelCounts.csv"
    if os.path.exists(output_file_path):
        print("[WARNING] Output file {} already exists, it will be " \
            "overwritten".format(output_file_path))

    file_results = {}
    for path in parameters:
        print("[INFO] Processing {}...".format(path))
        results = {}
        file_ending = path.split(".")[-1]
        line_processors = {
            "sam": process_sam_line,
            "vcf": process_vcf_line
        }
        if file_ending not in line_processors:
            raise_error("[INVALID FILE FORMAT] The format '{}' cannot be" \
                " processed".format(file_ending))

        process_line = line_processors[file_ending]
        file = open(path, "r")
        skip_file = False
        for line in file:
            if skip_file:
                continue
            line_result = process_line(line)
            # Error handling
            if line_result == -1:
                print("[MALFORMED FILE] The fields in {} are" \
                    " not tab separated or do not include all required" \
                    " information ".format(path))
                skip_file = True
            if line_result["chromosome"] == None:
                continue
            results = increase_values(results, line_result)
        file.close()
        file_results[path] = results
        print("[INFO] Done.")

    print("[INFO] Writing CSV file...")
    output_file = open(output_file_path, "w")
    csv_writer = csv.writer(output_file)
    header = ["file", "chromosome"]
    for count_name in empty_line_result()["counts"]:
        header.append(count_name)
    csv_writer.writerow(header)
    for path, chromosomes in file_results.items():
        for chromosome, counts in chromosomes.items():
            values = [path, chromosome]
            for count_key, count_value in counts.items():
                values.append(count_value)
            csv_writer.writerow(values)
    output_file.close()
    print("[INFO] Done.")

main()

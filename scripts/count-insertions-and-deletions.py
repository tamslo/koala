#!/usr/bin/python

import sys
import re

def raise_error(error_message):
    print(error_message)
    sys.exit(1)

def empty_counts():
    return {
        "insertions": 0,
        "inserted_bases": 0,
        "deletions": 0,
        "deleted_bases": 0,
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
    counts = empty_counts()
    if not line.startswith("@"):
        sam_fields = fields = line.split("\t")
        if len(sam_fields) < 11:
            return -1
        cigar_string = sam_fields[5]
        # Process insertions
        insertion_indices = [match.start() for match in re.finditer("I", cigar_string)]
        if len(insertion_indices) > 0:
            counts["insertions"] = len(insertion_indices)
            for index in insertion_indices:
                counts["inserted_bases"] += number_of_indel_bases(cigar_string, index)
        # Process deletions
        deletion_indices = [match.start() for match in re.finditer("D", cigar_string)]
        if len(deletion_indices) > 0:
            counts["deletions"] = len(deletion_indices)
            for index in deletion_indices:
                counts["deleted_bases"] += number_of_indel_bases(cigar_string, index)

    return counts

def process_vcf_line(line):
    counts = empty_counts()
    if not line.startswith("#"):
        vcf_fields = fields = line.split("\t")
        if len(vcf_fields) < 10:
            return -1
        reference_bases = vcf_fields[3]
        alternative_bases = vcf_fields[4]
        # Process insertions
        if len(alternative_bases) > len(reference_bases):
            counts["insertions"] = 1
            counts["inserted_bases"] = len(alternative_bases) - len(reference_bases)
        # Process deletions
        if len(reference_bases) > len(alternative_bases):
            counts["deletions"] = 1
            counts["deleted_bases"] = len(reference_bases) - len(alternative_bases)
    return counts

def main():
    parameters = sys.argv[1:]
    if len(parameters) < 1:
        raise_error("[NO PATH] The path to one or multiple files to be" \
            " processed is expected as parameters")

    for path in parameters:
        insertions = 0
        inserted_bases = 0
        deletions = 0
        deleted_bases = 0

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
        for line in file.readlines():
            counts = process_line(line)
            # Error handling
            if counts == -1:
                raise_error("[MALFORMED FILE] The fields in {} are" \
                    " not tab separated".format(path))

            insertions += counts["insertions"]
            inserted_bases += counts["inserted_bases"]
            deletions += counts["deletions"]
            deleted_bases += counts["deleted_bases"]
        file.close()

        print("")
        print("Processed file: {}".format(path))
        print("Number of insertions: {}".format(insertions))
        print("Number of inserted bases: {}".format(inserted_bases))
        print("Number of deletions: {}".format(deletions))
        print("Number of deleted bases: {}".format(deleted_bases))

main()

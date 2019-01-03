import csv
import os
import re

insertion_type = "I"
deletion_type = "D"

def get_fields(line, min_fields):
    fields = line.split("\t")
    if len(fields) < min_fields:
        print("[WARINING] Line is skipped because it is not tab separated")
        return None
    return fields

# def increase_values(values, line_values):
#     # TODO: remove duplicates
#     chromosome = line_values["chromosome"]
#     if not chromosome in values:
#         values[chromosome] = {}
#     for count_key, count_value in line_values["counts"].items():
#         if count_key in values[chromosome]:
#             values[chromosome][count_key] += count_value
#         else:
#             values[chromosome][count_key] = count_value
#     return values

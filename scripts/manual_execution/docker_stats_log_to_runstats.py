# Script to accumulate runstats from logged docker stats, takes path to
# docker stats as parameter

import sys
import os
import re

stats_path = sys.argv[1]

def get_mem_usage(mem_field):
    mem_usage_string = mem_field.split(" / ")[0]
    if mem_usage_string.endswith("KiB"):
        power = 1
    elif mem_usage_string.endswith("MiB"):
        power = 2
    elif mem_usage_string.endswith("GiB"):
        power = 3
    mem_usage_value = float(mem_usage_string[:-3])
    factor = 1024 ** power
    mem_usage_byte = mem_usage_value * factor
    return mem_usage_byte / 10**9

if os.path.exists(stats_path):
    print("Processing input file...")
    stats_file = open(stats_path, "r")
    line_num = 0
    mem_num = 0
    mem_sum = 0
    max_mem = 0
    for line in stats_file:
        line_num +=1
        fields = re.split(r'\s{2,}', line)
        cpu_usage = fields[2]
        if line_num % 2 == 1 or cpu_usage == "0.00%":
            continue
        mem_usage = get_mem_usage(fields[3])
        mem_num += 1
        mem_sum += mem_usage
        if mem_usage > max_mem:
            max_mem = mem_usage

    stats_file.close()
    output_path = stats_path + ".runstats"
    print("Writing output...")
    with open(output_path, "w") as output_file:
        avg_mem = mem_sum / mem_num
        output_file.write(
            "\nAverage memory usage: {} GB\n" \
            "Maximal memory usage: {} GB\n".format(
                round(avg_mem, 2),
                round(max_mem, 2)
            )
        )
    print("Done.")
else:
    print("Given path '{}' does not exist".format(stats_path))

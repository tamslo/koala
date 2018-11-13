import os, urllib.request

giab_version = "3.3.2"
confidence_sets = [
    # {
    #     "reference": "hg38",
    #     "name": "GRCh38",
    #     "bed_file_postfix": "_noCENorHET7"
    # },
    {
        "reference": "hg19",
        "name": "GRCh37"
    }
]

def progress(count, block_size, total_size):
    size = block_size * (count + 1)
    if size >= total_size:
        # Delete last progress if finished
        print("                ", end="\r", flush=True)
    else:
        percent = str(size / total_size * 100)[:4]
        print("Progress: {}%".format(percent), end="\r", flush=True)

def get_confidence_set(confidence_set, directory):
    url = "ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001" \
        "/NISTv{}/{}/".format(giab_version, confidence_set["name"])
    file_prefix = "HG001_{}_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X" \
        "-SOLID_CHROM1-X_v.{}_highconf".format(confidence_set["name"], giab_version)

    vcf_file_url = file_prefix + "_PGandRTGphasetransfer.vcf.gz"
    bed_file_postfix = ""
    if "bed_file_postfix" in confidence_set:
        bed_file_postfix = confidence_set["bed_file_postfix"]
    bed_file_url = file_prefix + "_nosomaticdel{}.bed".format(bed_file_postfix)

    vcf_zip_name = "confidence_calls.vcf.gz"
    download_path = directory + vcf_zip_name
    print("Downloading {}".format(vcf_zip_name), flush=True)
    urllib.request.urlretrieve(url + vcf_file_url, download_path, progress)
    os.system("gunzip {}".format(download_path))
    print("Done", flush=True)
    # os.remove(download_path)

    bed_file_name = "confidence_calls.bed"
    download_path = directory + bed_file_name
    print("Downloading {}".format(bed_file_name), flush=True)
    urllib.request.urlretrieve(url + bed_file_url, download_path, progress)
    print("Done", flush=True)

for confidence_set in confidence_sets:
    directory = "/giab/" + confidence_set["reference"] + "/"
    if not os.path.isdir(directory):
        os.makedirs(directory)
        get_confidence_set(confidence_set, directory)

# Import and Process SMART Data Sets

This manual explains how to import and process (SMART) data sets and documents how it was done initially.

## Import data sets

The data set import is accomplished by two scripts. The import could be done manually on the platform by adding data sets, however the import of big files still fails and this would take a long time. The scripts expect the raw data to have the following format:

```
.
+--+dataset_1
|  +--dataset_1_1.fq[.gz]
|  +--dataset_1_2.fq[.gz]
|  +--...
+--+dataset_2
|  +--dataset_2_1.fq[.gz]
|  +--dataset_2_2.fq[.gz]
|  +--...
+...
```

The folder name (`dataset_1` and `dataset_2`) is the prefix of containing file names, forward FASTQ files have the postfix `_1.fq`, reverse FASTQ files `_2.fq`. The folder name will later be used as the data set ID and name postfix (e.g., `smart_dataset_1` and `SMART dataset_1`). The `smart_` or `SMART` prefix can be changed in `import.py` to comply with other data sets to be imported.

For the import, the following steps need to be executed:

1. Create folder `data/intermediate`
2. Run Rsync (commands shown in `run.sh`) to transfer desired files to `data/intermediate` (maybe not all data sets fit on one machine)
3. Execute `run.sh`, it will extract gzipped files and remove not needed files to reduce file space.

The current setup is:

- Data sets 001-009 on VM2 (192.168.30.19)
- Data sets 010-019 on VM3 (192.168.30.176)
- Data sets 020-033 on VM4 (192.168.30.177)

# Run Alignment

The STAR alignment can be executed on the platform with the imported data sets. Multiple data sets can be selected that emerge in one experiment per data set. However, maybe not all experiments can be executed at once because of the disk space.

# Restrict Alignment by Coverage

To reduce run time for succeeding steps and obtain high-quality variants, the resulting BAM file is limited to regions with a high coverage. This is achieved by running the `services/giab/assets/restrict_bam.sh` script. The steps we conducted are presented as follows:

- Start Docker container with shell `docker run -v ~/code/data:/data -it giab` (paths may vary, the GIAB image is built with deployment; if not, run `docker build -t giab services/giab` from the `code` (repository root) directory)
- In the Docker console, run `./restrict_bam.sh 10 smart` (`10` is the desired coverage, `smart` the prefix of data sets that will be considered). Don't be distracted by errors saying that an integer cannot be parsed from a string in exponential notation; all these lines are included by default since the coverage usually is pretty high.

# Run Whole Pipeline

The whole pipeline can be executed with the Koala platform again, automatically using the Koala platform (you can select multiple data sets again). We ran STAR (it takes the restricted file from cache), GATK Best Practices (runs MarkDuplicates and SplitNCigarReads), GATK HaplotypeCaller (all chromosomes), and GATK VariantFiltration. The resulting VCF files can be downloaded or transferred and used in further experiments.

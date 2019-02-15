# Koala

An extensible platform to run and evaluate RNA-Seq analysis pipelines. The current version is accessible under https://vm-slosarek-01.eaalab.hpi.uni-potsdam.de.

## Setup

Make sure you have the following programs installed, versions given are those used for development (you can run the `scripts/setup.sh` script):

- node (10.2.1)
- npm (5.6.0)
- python (3.5.2)
- pip (8.1.1)
- docker (18.09.0)

Install requirements by running `./scripts/install.sh`, rename the `config.example.yml` to `config.yml`, and adapt the absolute path to the repository (keep the trailing dash).

To start the app in development mode, run `./scripts/dev-start.sh`, the app is accessible under https://localhost:3000.

To start the app in production mode, run `./scripts/prod-start.sh`, the app is accessible under https://localhost:5000.

## Deployment

The deployment server is `vm-slosarek-01.eaalab.hpi.uni-potsdam.de` (192.168.31.121), for a list of all VMs see below.

To deploy, if not done already, add the deploy remote `git remote add deploy deploy@192.168.31.121:git` and push the current version to the deploy remote `git push deploy`.

We are using a git remote on the deployment server that installs and runs the current version using the `post-receive` hook (for instructions see https://gist.github.com/noelboss/3fe13927025b89757f8fb12e9066f2fa). Virtualenv and NVM need to be installed on the deployment server. Make sure the post-receive hook is executable with `chmod a+x ~/git/hooks/post-receive`. To update the hook, run `cp ~/code/scripts/post-receive ~/git/hooks/post-receive`.

Additionally, we forward port 5000 to port 80 with `iptables -t nat -A PREROUTING -i eth0 -p tcp --dport 80 -j REDIRECT --to-port 5000`.

List of VMs:

- VM 1: 192.168.31.121
- VM 2: 192.168.30.19
- VM 3: 192.168.30.176
- VM 4: 192.168.30.177

## Adding a Service

Services are composed of a Dockerfile that defines its dependencies, a configuration file, and a Python class extending the `BaseService`. For one Dockerfile, multiple services can exist, each with an own configuration file and Python class.

Add new Dockerfiles to sub-directories in `services`, where the directory name will be the tag name of the Docker image that will be build from the file with the deploy script.

If only one service depends on this image, create a file called `config.yml`, if multiple services will depend on it, create a file called `<service_name>.config.yml`. The configuration file must include the `name` that will be displayed in the front-end and the `type` of the service, which is one of `aligner`, `alignment_filter`, or `variant_caller`. For some types, additional information is needed, as documented below.

To perform tasks, a class that inherits from `BaseService` needs to be added, look at similar services for examples. The class needs to be imported in `services/__init__.py` and added to the dictionary. For now, all classes are defined in the file `image_name/__init__.py` for shorter imports. If the files get too large, they can of course be split up.

**Important note for aligners**: Aligners are expected to have a maximum quality value of 60, usually options are provided to adapt it if this value is not the default. Additionally, the MD tag needs to be present in order to use Opossum. If no option is provided, SAMtools calmd can be used for post-processing (included in the GATK Docker image).

### Configuration

For services of type `aligner`, it needs to be defined whether the service creates a file by itself or writes to STDOUT by setting the boolean `creates_output` and whether the aligner uses soft clips by setting the boolean `soft_clips_exist`.

For services of type `evaluation`, a list of steps needs to be given that can be evaluated. Possible values are `alignment`, `alignment_filtering`, and `variant_calling` (as defined in `client/src/components/experimentConstants.js`).

## Manual tasks

Some tasks need to or can be executed manually. Also see the [README](scripts/import_smart_datasets/README.md) for SMART data processing.

### NovoAlign and BEDTools genomecov

We suppose Docker for Python cannot cope with the large log output, so these tasks need to be executed manually as explained [here](scripts/manual_execution/README.md).

### Import Data Sets

There are scripts in `scripts/import_smart_datasets` that potentially extract FASTQ files, move them to the right location, and create a JSON file holding the data set information.

Files are expected to be located in `data/intermediate`, with one folder per data set. The folder name is used for the data set's `id` and `name` that are prefixed with `smart_` or `SMART`. The folder is expected to contain paired FASTQ files (potentially gzipped), with the naming scheme `<folder name>_[1|2].fq`.

### Alignment Restriction with Coverage

This can be accomplished using BAMTools that are installed in the GIAB service. It creates a new BAM file, copies it to the default BAM location (Out.bam), and renames the old BAM for further use. Commands are defined in `services/giab/assets/restrict_bam.sh` that can be executed in a GIAB Docker container. It expects the `data` folder to be mounted as a volume under `/data`. It takes the following parameters:

- Minimum coverage (number)
- Prefix for datasets (optional, if not given, all datasets are included)

## Useful Commands for Testing

- Run docker on VM: `docker run -v ~/code/data:/data -it <container_name>`
- Build image manually: `docker build -t <image_name> services/<image_name>`

## Troubleshooting

If the app is not available...

- check if config.yml was created
- check that port forwarding works
- restart

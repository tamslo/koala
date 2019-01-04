# Koala

An extensible platform to run and evaluate RNA-Seq analysis pipelines. The current version is accessible under https://vm-slosarek-01.eaalab.hpi.uni-potsdam.de.

RefFlat annotations and BED files were downloaded from UCSC Genome Table Browser ([hg19](http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=212485719&clade=mammal&org=Human&db=hg19&hgta_group=genes&hgta_track=refGene&hgta_table=refFlat&hgta_regionType=genome&position=chr21%3A33031597-33041570&hgta_outputType=all&hgta_outFileName=refFlat.txt.gz)|[hg38](http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=212485719&clade=mammal&org=Human&db=hg38&hgta_group=genes&hgta_track=refGene&hgta_table=refFlat&hgta_regionType=genome&position=chr21%3A33031597-33041570&hgta_outputType=all&hgta_outFileName=refFlat.txt.gz)).

RNA editing sites for hg38 were translated from those for [hg19](http://lilab.stanford.edu/GokulR/database/Human_AG_all_hg19.bed) with [UCSC hgLiftOver](http://genome.ucsc.edu/cgi-bin/hgLiftOver).

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

The deployment server is `vm-slosarek-01.eaalab.hpi.uni-potsdam.de` (192.168.31.121).

To deploy, if not done already, add the deploy remote `git remote add deploy deploy@192.168.31.121:git` and push the current version to the deploy remote `git push deploy`.

We are using a git remote on the deployment server that installs and runs the current version using the `post-receive` hook (for instructions see https://gist.github.com/noelboss/3fe13927025b89757f8fb12e9066f2fa). Virtualenv and NVM need to be installed on the deployment server. Make sure the post-receive hook is executable with `chmod a+x ~/git/hooks/post-receive`. To update the hook, run `cp ~/code/scripts/post-receive ~/git/hooks/post-receive`.

Additionally, we forward port 5000 to port 80 with `iptables -t nat -A PREROUTING -i eth0 -p tcp --dport 80 -j REDIRECT --to-port 5000`.

## Adding a Service

Services are composed of a Dockerfile that defines its dependencies, a configuration file, and a Python class extending the `BaseService`. For one Dockerfile, multiple services can exist, each with an own configuration file and Python class.

Add new Dockerfiles to sub-directories in `services`, where the directory name will be the tag name of the Docker image that will be build from the file with the deploy script.

If only one service depends on this image, create a file called `config.yml`, if multiple services will depend on it, create a file called `<service_name>.config.yml`. The configuration file must include the `name` that will be displayed in the front-end and the `type` of the service, which is one of `aligner`, `alignment_filter`, or `variant_caller`. For some types, additional information is needed, as documented below.

To perform tasks, a class that inherits from `BaseService` needs to be added, look at similar services for examples. The class needs to be imported in `services/__init__.py` and added to the dictionary. For now, all classes are defined in the file `image_name/__init__.py` for shorter imports. If the files get too large, they can of course be split up.

**Important note for aligners**: Aligners are expected to have a maximum quality value of 60, usually options are provided to adapt it if this value is not the default. Additionally, the MD tag needs to be present in order to use Opossum. If no option is provided, SAMtools calmd can be used for post-processing (included in the GATK Docker image).

### Configuration

For services of type `aligner`, it needs to be defined whether the service creates a file by itself or writes to STDOUT by setting the boolean `creates_output` and whether the aligner uses soft clips by setting the boolean `soft_clips_exist`.

For services of type `evaluation`, a list of steps needs to be given that can be evaluated. Possible values are `alignment`, `alignment_filtering`, and `variant_calling` (as defined in `client/src/components/experimentConstants.js`).

## Useful Commands for Testing

- Run docker on VM: `docker run -v ~/code/data:/data -it <container_name>`
- Build image manually: `docker build -t <image_name> services/<image_name>`

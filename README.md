# Koala

An extensible platform to run and evaluate RNA-Seq analysis pipelines. The current version is accessible under https://vm-slosarek-01.eaalab.hpi.uni-potsdam.de.

## Setup

Make sure you have the following programs installed, versions given are those used for development:

- node (10.2.1)
- npm (5.6.0)
- python (3.5.2)
- pip (8.1.1)
- docker (17.12.1-ce)

Install requirements by running `./scripts/install.sh`, rename the `config.example.yml` to `config.yml`, and adapt the absolute path to the repository (keep the trailing dash).

To start the app in development mode, run `./scripts/dev-start.sh`, the app is accessible under https://localhost:3000.

To start the app in production mode, run `./scripts/prod-start.sh`, the app is accessible under https://localhost:5000.

## Deployment

The deployment server is `vm-slosarek-01.eaalab.hpi.uni-potsdam.de` (192.168.31.121).

To deploy, if not done already, add the deploy remote `git remote add deploy deploy@192.168.31.121:git` and push the current version to the deploy remote `git push deploy`.

We are using a git remote on the deployment server that installs and runs the current version using the `post-receive` hook (for instructions see https://gist.github.com/noelboss/3fe13927025b89757f8fb12e9066f2fa). Virtualenv and NVM need to be installed on the deployment server. Make sure the post-receive hook is executable with `chmod a+x ~/git/hooks/post-receive`. To update the hook, run `cp ~/code/scripts/post-receive ~/git/hooks/post-receive`.

Additionally, we forward port 5000 to port 80 with `iptables -t nat -A PREROUTING -i eth0 -p tcp --dport 80 -j REDIRECT --to-port 5000`.

## Adding a Service

To add a service, add a directory in services, its name will be the `service_id`. To this directory, add a Dockerfile that sets up a service and a YAML with the specification (at least the name that should be displayed in the front-end and the type should be included). To perform tasks, a class that inherits from `BaseService` needs to be added, look at similar services for examples.

## Useful Commands for Testing

- Run docker locally: `docker run -v /c/Users/Tamara/Repos/koala/data:/data -it <container_name>`
- Run docker on VM: `docker run -v ~/code/data:/data -it <container_name>`
- Build test image: `docker build -t test_aligner services/test_aligner`

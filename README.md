# Koala

An extensible platform to run and evaluate RNA-Seq analysis pipelines. The current version is accessible under https://vm-slosarek-01.eaalab.hpi.uni-potsdam.de:5000.

## Setup

Make sure you have the following programs installed, versions given are those used for development:

* node (10.0.0)
* npm (5.6.0)
* python (3.5.2)
* pip (8.1.1)
* docker (1.13.1)
* docker-compose (1.8.0)

Install requirements by running `./scripts/install.sh`.

To start the app in development mode, run `./scripts/dev-start.sh`, the app is accessible under https://localhost:3000.

To start the app in production mode `./scripts/prod-start.sh`, the app is accessible under https://localhost:5000.

## Deployment

**TODO**: Needs to be set up

The deployment server is `vm-slosarek-01.eaalab.hpi.uni-potsdam.de` (192.168.31.121).

To deploy, if not done already, add the deploy remote `git remote add release deploy@192.168.31.121:git` and push the current version to the deploy remote `git push release`.

We are using a git remote on the deployment server `vm-slosarek-01.eaalab.hpi.uni-potsdam.de` (for instructions see https://gist.github.com/noelboss/3fe13927025b89757f8fb12e9066f2fa) that deploys the current version using the `post-receive` hook. Make sure the post-receive hook is executable with `chmod a+x ~/git/hooks/post-receive`.

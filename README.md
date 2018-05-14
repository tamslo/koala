# Koala

An extensible platform to run and evaluate RNA-Seq analysis pipelines. The current version is accessible under https://vm-slosarek-01.eaalab.hpi.uni-potsdam.de:5000.

## Setup

Make sure you have the following programs installed, versions given are those used for development:

* node (10.0.0)
* npm (5.6.0)
* python3 (3.5.2)
* pip3 (8.1.1)
* docker (1.13.1)
* docker-compose (1.8.0)

To install required modules, run `scripts/install.sh`.

To start the app in development mode, run `scripts/start-dev.sh`, the app is accessible under https://localhost:3000.

To start the app in production mode `scripts/start-prod.sh`, , the app is accessible under https://localhost:5000.

## Deployment

**TODO**: Needs to be set up

To deploy, if not done already, add the deploy remote `git remote add deploy deploy@vm-slosarek-01.eaalab.hpi.uni-potsdam.de:git` and push the current version to the deploy remote `git push deploy`.

We are using a git remote on the deployment servervm-slosarek-01.eaalab.hpi.uni-potsdam.de (for instructions see https://gist.github.com/noelboss/3fe13927025b89757f8fb12e9066f2fa) that deploys the current version using the `post-receive` hook. Copy the hook to the server with `cp ~/code/scripts/post-receive ~/git/hooks/post-receive` and make it executable with `chmod a+x ~/git/hooks/post-receive`.

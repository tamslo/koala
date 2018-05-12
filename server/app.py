from flask import Flask, request
from flask_cors import CORS
from collections import OrderedDict
from cache import Cache
from experiments import Experiments
import json, urllib, sys

app = Flask(__name__)
CORS(app)

config = json.load(open("../config.json", "r"), object_pairs_hook=OrderedDict)
data_directory = "../data/"

cache = Cache(data_directory)
experiments = Experiments(data_directory)

@app.route("/context", methods=["GET"])
def get_context():
    return json.dumps({
        "aligners": config["aligners"],
        "experiments": experiments.all()
    })

@app.route("/experiment", methods=["POST"])
@app.route("/experiment/<id>", methods=["DELETE"])
def experiment(id = None):
    if request.method == "POST":
        params = request.get_json()
        experiment = experiments.create(params)
        return json.dumps(experiment)
    else:
        return json.dumps(experiments.delete(id))

@app.route("/data/<id>", methods=["GET"])
def data(id):
    experiment = experiments.select(id)
    experiment_data = cache.get_experiment_data(experiment)
    try:
        if "dataset" in experiment_data:
            dataset_path = experiment_data["dataset"]
            experiment = experiments.add_log_entry(experiment, "load data", one_step = True)
        else:
            experiment = experiments.add_log_entry(experiment, "download data")
            get_data(experiment["dataset"])
            experiment = experiments.log_complete(experiment, "download data")
    except Exception as error:
        experiment = experiments.mark_error(id, error)
    return json.dumps(experiment)

@app.route("/done/<id>", methods=["GET"])
def done(id):
    return json.dumps(experiments.mark_done(id))

def get_data(url):
    file_name, headers = urllib.request.urlretrieve(url, cache.create_dataset(url))
    # TODO throw error if headers not okay
    app.logger.info(headers)
    return file_name

if __name__ == "__main__":
    app.env = "development"
    app.debug = True
    app.run()

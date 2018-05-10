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
def experiment():
    params = request.get_json()
    return experiments.create(params)

@app.route("/data", methods=["GET"])
def data():
    experiment_id = request.args.get("experiment")
    experiment = experiments.select(experiment_id)
    experiment_data = cache.get_experiment_data(experiment)
    try:
        if "dataset" in experiment_data:
            dataset_path = experiment_data["dataset"]
            experiment = experiments.add_status(experiment_id, "Loaded data")
        else:
            dataset_path = get_data(experiment["dataset"])
            experiment = experiments.add_status(experiment_id, "Downloaded data")
    except Exception as error:
        experiment = experiments.mark_error(experiment_id, error)
    return experiment

@app.route("/done", methods=["GET"])
def done():
    experiment_id = request.args.get("experiment")
    experiment = experiments.mark_done(experiment_id)
    return experiment

def get_data(url):
    file_name, headers = urllib.request.urlretrieve(url, cache.create_dataset(url))
    app.logger.info(file_name)
    # TODO throw error if headers not okay
    app.logger.info(headers)
    return file_name

if __name__ == "__main__":
    app.env = "development"
    app.debug = True
    app.run()

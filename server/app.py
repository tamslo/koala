from flask import Flask, request
from flask_cors import CORS
from collections import OrderedDict
from cache import Cache
from experiments import Experiments
import json, urllib, sys

app = Flask(__name__)
CORS(app)

config = json.load(open("../config.json", "r"), object_pairs_hook=OrderedDict)
data_directory = "data/"

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

@app.route("/run", methods=["GET"])
def run():
    experiment_id = request.args.get("experiment")
    print("ID")
    print(experiment_id)
    sys.stdout.flush()
    experiment = experiments.select(experiment_id)
    experiment_data = cache.get_experiment_data(experiment)
    dataset_path = experiment_data["dataset"] or get_data(experiment["dataset"])
    # TODO run experiment, evaluate, return report
    report = {}
    return json.dumps(report)

def get_data(url):
    file_name, headers = urllib.request.urlretrieve(url, cache.create_dataset(url))
    app.logger.info(file_name)
    app.logger.info(headers)
    return file_name

if __name__ == "__main__":
    app.env = "development"
    app.debug = True
    app.run()

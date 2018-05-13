from flask import Flask, request
from flask_cors import CORS
from collections import OrderedDict
from cache import Cache
from experiments import Experiments
import json, urllib, atexit

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

@app.route("/experiment", methods=["POST", "PUT"])
@app.route("/experiment/<id>", methods=["DELETE"])
def experiment(id = None):
    if request.method == "POST":
        params = request.get_json()
        experiment = experiments.create(params)
        return json.dumps(experiment)
    elif request.method == "PUT":
        params = request.get_json()
        experiment = experiments.update(params)
        return json.dumps(experiment)
    else:
        return json.dumps(experiments.delete(id))

@app.route("/data/<id>", methods=["GET"])
def data(id):
    experiment = experiments.select(id)
    experiment_data = cache.get_experiment_data(experiment)
    try:
        if experiment_data["dataset"]:
            dataset_path = experiment_data["dataset"]
            experiment = experiments.add_log_entry(experiment, "cached", one_step = True)
        else:
            experiment = experiments.add_log_entry(experiment, "dataset")
            get_data(experiment["dataset"])
            experiment = experiments.log_complete(experiment, "dataset")
    except Exception as error:
        cache.clean_up("dataset", experiment["dataset"])
        experiment = experiments.mark_error(id, error)
    return json.dumps(experiment)

@app.route("/done/<id>", methods=["GET"])
def done(id):
    return json.dumps(experiments.mark_done(id))

def get_data(url):
    file_name, headers = urllib.request.urlretrieve(url, cache.create_dataset(url))
    return file_name

def clean_up():
    # Not working properly if debug=True, see https://stackoverflow.com/questions/37064595/handling-atexit-for-multiple-app-objects-with-flask-dev-server-reloader
    for experiment_id, experiment in experiments.all().items():
        last_log_entry = experiment["log"][-1]
        if (
            experiment["done"] or
            experiment["error"] or
            last_log_entry["completed"]
        ):
            continue
        else:
            experiments.mark_interrupted(experiment_id)
            action = last_log_entry["action"]
            if action in experiment:
                id = experiment[action]
                try:
                    cache.clean_up(action, id)
                except Exception as error:
                    app.logger.error("Manual cache-cleanup needed for {} {} ({})".format(
                        action,
                        id,
                        error
                    ))

atexit.register(clean_up)

if __name__ == "__main__":
    app.run(host = "0.0.0.0", debug = True)

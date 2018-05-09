from flask import Flask, request
from flask_cors import CORS
from collections import OrderedDict
from cache import Cache
import json, urllib, sys

app = Flask(__name__)
CORS(app)

config = json.load(open("../config.json", "r"), object_pairs_hook=OrderedDict)
cache = Cache()

@app.route("/context", methods=["GET"])
def get_context():
    return json.dumps({
        "aligners": config["aligners"],
        "experiments": [] # TODO store and load
    })

@app.route("/run", methods=["POST"])
def run():
    params = request.get_json()
    experiment = cache.get_experiment(params)
    dataset = experiment["dataset"] or get_data(params["dataset"])
    # TODO run experiment, evaluate, return report
    return json.dumps(experiment)

def get_data(url):
    file_name, headers = urllib.request.urlretrieve(url, cache.create_dataset(url))
    app.logger.info(file_name)
    app.logger.info(headers)
    return file_name

if __name__ == "__main__":
    app.env = "development"
    app.debug = True
    app.run()

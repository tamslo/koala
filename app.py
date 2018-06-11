import json, atexit, zipfile, time, os, optparse, shutil
from flask import Flask, request, send_file, send_from_directory
from flask_cors import CORS
from collections import OrderedDict
from modules.datasets import Datasets
from modules.experiments import Experiments
from modules.runner import Runner
from modules.downloader import Downloader

app = Flask(__name__)
CORS(app)

services = json.load(open("services.json", "r"), object_pairs_hook=OrderedDict)
data_directory = "data/"
download_directory = "data/tmp/"

datasets = Datasets(data_directory)
experiments = Experiments(data_directory)
runner = Runner(datasets, experiments, data_directory)
downloader = Downloader(download_directory)

@app.route("/ping")
def ping():
    return time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())

@app.route("/context", methods=["GET"])
def get_context():
    return json.dumps({
        "services": services,
        "experiments": experiments.all(),
        "datasets": datasets.get_datasets()
    })

@app.route("/dataset", methods=["POST"])
def dataset():
    params = request.get_json()
    dataset = datasets.create(params)
    return json.dumps(dataset)

@app.route("/experiment", methods=["POST", "PUT", "DELETE"])
def experiment():
    if request.method == "POST":
        params = request.get_json()
        experiment = experiments.create(params)
        return json.dumps(experiment)
    elif request.method == "PUT":
        params = request.get_json()
        experiment = experiments.update(params)
        return json.dumps(experiment)
    else:
        experiment_id = request.args.get("id")
        return json.dumps(experiments.delete(experiment_id))

@app.route("/execute", methods=["GET"])
def data():
    action = request.args.get("action")
    experiment_id = request.args.get("experiment")
    experiment = runner.execute(action, experiment_id)
    return json.dumps(experiment)

@app.route("/done", methods=["GET"])
def done():
    experiment_id = request.args.get("experiment")
    return json.dumps(experiments.mark_done(experiment_id))

@app.route("/download", methods=["GET"])
def download():
    experiment_id = request.args.get("experiment")
    path = request.args.get("path")
    if path != None:
        download_path = path
        download_name = downloader.file_name(path)
    elif experiment_id != None:
        experiment = experiments.select(experiment_id)
        download_path, download_name = downloader.zip(experiment)
    return send_file(
        download_path,
        as_attachment=True,
        attachment_filename=download_name
    )


# Routes for client built with `npm run build`

@app.route("/")
def serve():
    return send_from_directory("client/build/", "index.html")

@app.route("/koala.ico")
def servefav():
    return send_from_directory("client/build/", "koala.ico")

@app.route("/static/js/<path:path>")
def servejs(path):
    return send_from_directory("client/build/static/js/", path)

@app.route("/static/media/<path:path>")
def servemedia(path):
    return send_from_directory("client/build/static/media/", path)


def clean_up():
    shutil.rmtree(download_directory)
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
                try:
                    datasets.clean_up(action, experiment)
                except Exception as error:
                    app.logger.error("Manual cleanup of datasets needed for {} {} ({})".format(
                        action,
                        id,
                        error
                    ))

atexit.register(clean_up)

if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option("-d", "--debug", action="store_true", dest="debug", help=optparse.SUPPRESS_HELP)
    options, _ = parser.parse_args()
    app.run(host="0.0.0.0", port="5000", debug=bool(options.debug))

import json, urllib, atexit, zipfile, time, os, optparse
from flask import Flask, request, send_file, send_from_directory
from flask_cors import CORS
from collections import OrderedDict
from cache import Cache
from experiments import Experiments

app = Flask(__name__)
CORS(app)

config = json.load(open("config.json", "r"), object_pairs_hook=OrderedDict)
data_directory = "data/"
archive_directory = data_directory + "zips/"
if not os.path.isdir(archive_directory):
    os.makedirs(archive_directory)

cache = Cache(data_directory)
experiments = Experiments(data_directory)

@app.route("/ping")
def ping():
    return time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())

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
            dataset_path, headers = urllib.request.urlretrieve(
                experiment["dataset"],
                cache.create_dataset(experiment["dataset"])
            )
            experiment = experiments.log_complete(experiment, "dataset")
        experiment = experiments.add_download(experiment, "dataset", dataset_path)
    except Exception as error:
        cache.clean_up("dataset", experiment["dataset"])
        experiment = experiments.mark_error(id, error)
    return json.dumps(experiment)

@app.route("/done/<id>", methods=["GET"])
def done(id):
    return json.dumps(experiments.mark_done(id))

@app.route("/download/all", methods=["GET"])
def download_all():
    experiment_id = request.args.get("experiment")
    experiment = experiments.select(experiment_id)
    archive_name = experiment["name"] + ".zip"
    archive_path = archive_directory + archive_name
    archive = zipfile.ZipFile(archive_path, "w")
    for key, url in experiment["downloads"].items():
        path = url.split("path=")[1]
        file_name = path.split("/")[-1]
        archive.write(path, file_name)
    archive.close()
    return send_file(archive_path, as_attachment=True, attachment_filename=archive_name)

@app.route("/download", methods=["GET"])
def download():
    path = request.args.get("path")
    file_name = path.split("/")[-1]
    return send_file(path, as_attachment=True, attachment_filename=file_name)

# Routes for client built with `npm run build`

@app.route("/")
def serve():
    return send_from_directory("client/build/", "index.html")

@app.route("/static/js/<path:path>")
def servejs(path):
    return send_from_directory("client/build/static/js/", path)

@app.route("/static/media/<path:path>")
def servemedia(path):
    return send_from_directory("client/build/static/media/", path)


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
    parser = optparse.OptionParser()
    parser.add_option("-d", "--debug", action="store_true", dest="debug", help=optparse.SUPPRESS_HELP)
    options, _ = parser.parse_args()
    app.run(host = "0.0.0.0", port = "5000", debug = bool(options.debug))

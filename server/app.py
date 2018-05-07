import os, json, urllib, uuid
from flask import Flask, request
from collections import OrderedDict

app = Flask(__name__)

config = json.load(open("../config.json", "r"), object_pairs_hook=OrderedDict)

@app.route("/data", methods=["GET"])
def data():
    data_directory = "../data/"
    cache_index_path = data_directory + "index.json"

    if not os.path.isdir(data_directory):
        os.mkdir(data_directory)
        cache_index = open(cache_index_path, "w")
        cache_index.write(json.dumps({}));
        cache_index.close()

    cached_data = json.load(open(cache_index_path, "r"))
    url = request.args.get("url")

    if url in cached_data:
        file_path = cached_data[url]
    else:
        file_path = data_directory + str(uuid.uuid4())
        cached_data[url] = file_path;
        urllib.request.urlretrieve(url, file_path)
        cache_index = open(data_directory + cache_index, "w")
        cache_index.write(json.dumps(cached_data))
        cache_index.close()

    file_content = open(file_path, "r").read()
    return file_content

if __name__ == "__main__":
    app.env = "development"
    app.debug = True
    app.run()

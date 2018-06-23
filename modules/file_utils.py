import os, json, urllib

def create_directory(path):
    if not os.path.isdir(path):
        os.makedirs(path)

def write(content, path):
    with open(path, "w") as file:
        file.write(content)

def delete(path):
    try:
        os.remove(path)
    except OSError:
        pass

def download(url, destination):
    urllib.request.urlretrieve(
        url,
        destination
    )

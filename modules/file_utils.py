import os, json, urllib.request

def create_directory(path):
    if not os.path.isdir(path):
        os.makedirs(path)

def create_file(path):
    with open(path, "w") as file:
        return None

def write(content, path):
    with open(path, "w") as file:
        file.write(content)

def delete(path):
    try:
        os.remove(path)
    except OSError:
        pass

def download(url, destination):
    def progress(count, block_size, total_size):
        size = block_size * (count + 1)
        if size >= total_size:
            # Delete last progress if finished
            print("                ", end="\r", flush=True)
        else:
            percent = str(size / total_size * 100)[:4]
            print("Progress: {}%".format(percent), end="\r", flush=True)

    urllib.request.urlretrieve(
        url,
        destination,
        progress
    )

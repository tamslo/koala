import os, shutil, json, urllib.request

def validate_file_content(path):
    file_exists = os.path.exists(path)
    file_is_empty = file_exists and os.stat(path).st_size == 0
    if file_is_empty or not file_exists:
        raise Exception("Output {} not written".format(path))

def create_directory(path):
    if not os.path.isdir(path):
        os.makedirs(path)

def create_file(path):
    file = open(path, "w")
    file.close()

def write(content, path):
    with open(path, "w") as file:
        file.write(content)

def delete(path):
    if os.path.isdir(path):
        shutil.rmtree(path, ignore_errors=True)
    else:
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

def unzip(path):
    directory = path.rsplit("/", 1)[0]
    if path.endswith(".tar.bz2"):
        os.system("tar xjC {} -f {}".format(directory, path))
        unzipped_path = path.replace(".tar.bz2", "")
    elif path.endswith(".tar.gz"):
        os.system("tar xzC {} -f {}".format(directory, path))
        unzipped_path = path.replace(".tar.gz", "")
    elif path.endswith(".gz"):
        os.system("gunzip {}".format(path))
        unzipped_path = path.replace(".gz", "")
    else:
        # Unkown compression
        print("[file_utils] Unknown compression type for {}".format(path))
        return path

    return unzipped_path

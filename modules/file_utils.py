import os, shutil, json, urllib.request

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
    zip_directory = path.rsplit("/", 1)[0]
    path_parts = path.split(".", 1)
    unzipped_path = path_parts[0]
    zip_type = path_parts[1]

    if zip_type == "tar.bz2":
        os.system("tar xjC {} -f {}".format(zip_directory, path))
    else:
        # Unkown compression
        print("[file_utils] Unknown compression type: {}".format(zip_type))
        return path

    return unzipped_path

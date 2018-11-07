import docker
import yaml
import time
import datetime
import threading
import modules.file_utils as file_utils

class Docker:
    def __init__(self, data_directory):
        self.docker_client = docker.from_env()
        with open("config.yml", "r") as config_file:
            config = yaml.load(config_file)
            self.absolute_data_path = config["absolute_repo_path"] + data_directory

    def run(self, container_name, command, auto_remove=True, **kwargs):
        if command == None or command == "":
            print("[WARNING] Executing an empty command in {}".format(container_name), flush=True)

        return self.docker_client.containers.run(
            container_name,
            command,
            volumes={
                self.absolute_data_path: {
                    "bind": "/data",
                    "mode": "rw"
                }
            },
            auto_remove=auto_remove,
            **kwargs
        )

    def run_and_write_logs(self, container_name, command, stdout_file_path=None,
        stderr_file_path=None, stats_file_path=None):

        def write_stats(text):
            if stats_file_path != None:
                with open(stats_file_path, "a") as stats_file:
                    stats_file.write("{}\n".format(text))

        def log_output(docker_container):
            if stdout_file_path != None:
                out_file = open(stdout_file_path, "wb")
                for line in docker_container.logs(stdout=True, stderr=False, stream=True):
                    out_file.write(line)
                out_file.close()
            if stderr_file_path != None:
                log_file = open(stderr_file_path, "wb")
                for line in docker_container.logs(stdout=False, stderr=True, stream=True):
                    log_file.write(line)
                log_file.close()
            docker_container.reload()
            if docker_container.status != "exited":
                docker_container.stop()
            # If log file is empty, delete it. If only one file is written it is
            # expected to be the log file, if both files are written, it is expected
            # to be stderr_file_path
            log_file_path = stderr_file_path or stdout_file_path
            if not file_utils.file_has_content(log_file_path):
                file_utils.delete(log_file_path)

        def log_stats(docker_container):
            def gigabytes_string(bytes):
                return "{} GB".format(round(bytes / 10**9, 1))

            max_memory_usage = 0
            sum_memory_usage = 0
            num_memory_usages = 0
            while docker_container.status != "exited":
                docker_container.reload()
                stats = docker_container.stats(stream=False)
                if len(stats["memory_stats"].keys()) > 0:
                    max_memory_usage = stats["memory_stats"]["max_usage"]
                    sum_memory_usage += stats["memory_stats"]["usage"]
                    num_memory_usages += 1

            if num_memory_usages > 0:
                avg_memory_usage = sum_memory_usage / num_memory_usages
                write_stats("Average memory: {}".format(
                    gigabytes_string(avg_memory_usage)
                ))
                write_stats("Maximum memory: {}".format(
                    gigabytes_string(max_memory_usage)
                ))

        write_stats("{}\n".format(command))
        start_time = time.time()
        write_stats("Start time: {}".format(
            time.strftime("%d %b %Y %H:%M:%S", time.localtime(start_time))
        ))

        docker_container = self.run(
            container_name,
            command,
            stderr=stderr_file_path != None,
            detach=True,
            auto_remove=False
        )
        output = threading.Thread(target=log_output, args=(docker_container,))
        stats = threading.Thread(target=log_stats, args=(docker_container,))
        output.start()
        stats.start()
        output.join()
        stats.join()
        docker_container.remove()

        end_time = time.time()
        write_stats("End time: {}".format(
            time.strftime("%d %b %Y %H:%M:%S", time.localtime(end_time))
        ))
        write_stats("Total time: {}\n".format(
            str(datetime.timedelta(seconds=end_time - start_time))
        ))

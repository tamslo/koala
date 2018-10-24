import os
import yaml
from .base_service import BaseService

# Register services with functionality
from .test.aligner import TestAlignerWritesFile
from .test.aligner import TestAlignerWritesLog
from .star.aligner import Star
from .novoalign.aligner import NovoAlign
ServiceClasses = {
    "test_aligner_writes_file": TestAlignerWritesFile,
    "test_aligner_writes_log": TestAlignerWritesLog,
    "star": Star,
    "novoalign": NovoAlign
}

def getServices():
    services = []
    for directory in os.listdir("services"):
        service_directory = "services/{}/".format(directory)
        if os.path.isdir(service_directory):
            for file in os.listdir(service_directory):
                if file.endswith("config.yml"):
                    config_path = service_directory + file
                    if file == "config.yml":
                        service_id = directory
                    else:
                        service_id = directory + "_" + file.split(".config.yml")[0]
                    if service_id in ServiceClasses:
                        ServiceClass = ServiceClasses[service_id]
                    else:
                        ServiceClass = BaseService
                    services.append(ServiceClass(config_path, service_id))
    return services

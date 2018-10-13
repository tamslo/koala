import os
import yaml
from .base_service import BaseService

# Register services with functionality
from .test_aligner.aligner import TestAligner
from .star.aligner import Star
from .novoalign.aligner import NovoAlign
ServiceClasses = {
    "test_aligner": TestAligner,
    "star": Star,
    "novoalign": NovoAlign
}

services = []
for service in os.listdir("services"):
    service_directory = "services/{}/".format(service)
    config_path = service_directory + "config.yml"
    if os.path.isdir(service_directory) and os.path.exists(config_path):
        if service in ServiceClasses:
            ServiceClass = ServiceClasses[service]
        else:
            ServiceClass = BaseService
        services.append(ServiceClass(config_path, service))

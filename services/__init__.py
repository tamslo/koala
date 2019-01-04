import os
import yaml
from .base_service import BaseService

# Register services

from .test import TestAlignerWritesFile
from .test import TestAlignerWritesLog
from .star import Star
from .novoalign import NovoAlign
from .novoalign import NovoAlignIndelSensitive
from .gatk import GatkFilters
from .opossum import Opossum
from .beers import BeersEvaluator
from .gatk import HaplotypeCaller
from .gatk import HaplotypeCallerFiltered
from .giab import GiabEvaluator

with open("config.yml") as config_file:
    config = yaml.load(config_file)
    environment = config["environment"]

ServiceClasses = {
    "gatk_filters": GatkFilters,
    "opossum": Opossum,
    "beers": BeersEvaluator,
    "gatk_haplotypecaller": HaplotypeCaller,
    "gatk_haplotypecaller_chr1": HaplotypeCallerFiltered,
    "gatk_haplotypecaller_chr3": HaplotypeCallerFiltered,
    "gatk_haplotypecaller_chr4": HaplotypeCallerFiltered,
    "gatk_haplotypecaller_chr5": HaplotypeCallerFiltered,
    "gatk_haplotypecaller_chr17": HaplotypeCallerFiltered,
    "gatk_haplotypecaller_chr21": HaplotypeCallerFiltered,
    "giab": GiabEvaluator
}

if environment == "test":
    ServiceClasses["test_aligner_writes_file"] = TestAlignerWritesFile
    ServiceClasses["test_aligner_writes_log"] = TestAlignerWritesLog
    ServiceClasses["novoalign"] = NovoAlign # for testing

if environment == "production":
    ServiceClasses["star"] = Star
    ServiceClasses["novoalign"] = NovoAlign
    ServiceClasses["novoalign_tweaked"] = NovoAlignIndelSensitive

def get_services():
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
                        services.append(ServiceClass(config_path, service_id, image_name=directory))
    return services

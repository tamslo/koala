#!/bin/bash

function send_request() {
  experiment=$1
  curl --request POST \
    --header "Content-Type: application/json" \
    --data "$experiment" \
    "http://localhost:5000/experiment"
}

send_request '{
  "name": "STAR Human",
  "reference": "hg19",
  "datasets": ["simulated_reads_HG19t1r1", "simulated_reads_HG19t2r1", "simulated_reads_HG19t3r1"],
  "pipeline": {"alignment": "star"}
}'

send_request '{
  "name": "NovoAlign Human",
  "reference": "hg19",
  "datasets": ["simulated_reads_HG19t1r1", "simulated_reads_HG19t2r1", "simulated_reads_HG19t3r1"],
  "pipeline": {"alignment": "novoalign"}
}'

send_request '{
  "name": "STAR Malaria",
  "reference": "pfal",
  "datasets": ["simulated_reads_PFALt1r1", "simulated_reads_PFALt2r1", "simulated_reads_PFALt3r1"],
  "pipeline": {"alignment": "star"}
}'

send_request '{
  "name": "NovoAlign Malaria",
  "reference": "pfal",
  "datasets": ["simulated_reads_PFALt1r1", "simulated_reads_PFALt2r1", "simulated_reads_PFALt3r1"],
  "pipeline": {"alignment": "novoalign"}
}'

exit 0

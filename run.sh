#!/bin/bash

bash scripts/install.sh
cp config.example.yml config.yml
python app.py

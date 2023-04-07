#! /usr/bin/env python
import yaml

config_file = "config.yml"
with open(config_file, "r") as ymlfile:
    cfg = yaml.safe_load(ymlfile)

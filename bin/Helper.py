#!/usr/bin/env python
import yaml
import os


def import_yaml(yaml_path: os.path):
    """
    Opens yaml file containing hyper parameters.

    :param yaml_path: File path to yaml
    :return: dictionary with parameters
    """
    try:
        with open(yaml_path, 'r') as stream:
            return yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)


def save_yaml(d, filepath):
    with open(filepath, 'w') as file:
        documents = yaml.dump(d, file)


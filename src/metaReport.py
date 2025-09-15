#!/usr/bin/env python

from json2html import *
import argparse
#TODO fix import issue from Helper import import_yaml
import os
import yaml

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

def generate_report(input, output):

    data = import_yaml(input)
    # Example JSON data
    data_example = {
        "name": "John Doe",
        "age": 30,
        "address": {
            "street": "123 Main St",
            "city": "New York",
            "state": "NY"
        },
        "hobbies": ["Reading", "Traveling", "Swimming"]
    }

    # Convert JSON to HTML
    html_content = json2html.convert(json=data)

    html_template = f"""
    <html>
    <head>
        <title>JSON Data Report</title>
        <style>
            body {{ font-family: Arial, sans-serif; }}
            table {{ width: 100%; border-collapse: collapse; }}
            th, td {{ padding: 8px 12px; border: 1px solid #ddd; }}
            th {{ background-color: #f2f2f2; }}
        </style>
    </head>
    <body>
        <h2>Report</h2>
        {html_content}
    </body>
    </html>
    """

    # Save the styled HTML content to a file
    with open(output, "w") as html_file:
        html_file.write(html_template)

    print("Styled HTML report generated successfully.")


    print("HTML report generated successfully.")


def parse_arguments():
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument('--params', required=False, help='Topo file', default='../demo_plotly/config/params.yml')
    parser.add_argument('--sims', required=False, help='Topo file', default='../demo_plotly/config/params.yml')


    # Output
    parser.add_argument('--output', required=False, default='test/metaReport.html', help='')

    return parser.parse_args()

# Example of running the function
if __name__ == '__main__':
    args = parse_arguments()

    generate_report(args.params, args.output)

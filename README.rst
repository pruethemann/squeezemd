=========
squeezeMD - A Comprehensive Molecular Dynamics Analysis Workflow
=========


.. image:: https://img.shields.io/pypi/v/squeezemd.svg
        :target: https://pypi.python.org/pypi/squeezemd

.. image:: https://img.shields.io/travis/pruethemann/squeezemd.svg
        :target: https://travis-ci.com/pruethemann/squeezemd

.. image:: https://readthedocs.org/projects/squeezemd/badge/?version=latest
        :target: https://squeezemd.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status


* Free software: GNU General Public License v3
* Documentation: https://squeezemd.readthedocs.io.


Summary
--------
This workflow provides an integrated solution for conducting comprehensive molecular dynamics (MD) analysis. It encompasses a range of functionalities from mutation analysis, MD simulations, explorative trajectory analysis, to interaction fingerprinting and global interaction analysis. Designed to streamline the process of analyzing complex molecular dynamics simulations, this workflow offers tools for detailed examination of molecular interactions, stability, and conformational changes.

Detailed Summary
Snakefile: Serves as the backbone of the workflow, orchestrating the execution of various analysis scripts. It ensures that the analysis pipeline is executed efficiently and in the correct order, managing dependencies between different analysis steps.

Mutagensis
~~~~~~~~~~
This script performs mutation analysis, identifying and characterizing the effects of mutations on the structure and function of molecules. It is essential for understanding the impact of specific amino acid changes on protein stability and activity.

MD
~~~~~~~~~~
Conducts molecular dynamics simulations, providing insights into the physical movements and conformational changes of molecules over time. This script is crucial for exploring the dynamic nature of molecular systems.

Explorative Trajectory Analysis
~~~~~~~~~~
Offers tools for in-depth trajectory analysis, enabling the identification of key molecular events and interactions throughout the simulation. It aids in uncovering the mechanisms driving molecular behavior.

Centering Trajectorys
~~~~~~~~~~

Centers and aligns molecular dynamics trajectories, facilitating the comparison and analysis of simulation results. This preprocessing step is vital for accurate analysis of molecular motions.

Interaction Analysis
~~~~~~~~~~
A specialized tool for analyzing molecular interactions, focusing on the contributions of individual atoms and residues to the overall interaction network within the molecule.
Expands on the analysis of molecular interactions by examining global interaction patterns, providing a comprehensive view of the interaction landscape within the molecular system.

Interaction Fingerprints
~~~~~~~~~~
Generates interaction fingerprints, summarizing the interaction patterns between molecules in a compact, easily interpretable format. This tool is useful for comparing molecular interactions across different simulations or conditions.
Offers advanced analysis of interaction fingerprints, facilitating the comparison of global interaction patterns across multiple simulations. It aids in identifying consistent or aberrant interaction motifs.

Interaction Surface
~~~~~~~~~~
Focuses on the analysis of interaction surfaces, characterizing the areas of molecules involved in binding or interaction with other molecules. This script is essential for understanding the molecular basis of recognition and binding processes.

Demo workflow
----

1. The workflow can be tested by performing the following commands:
```bash
cd demo
# Perform a dry run
squeeze --resources gpu=1 -j4 -n

# Perform the demo production run
squeeze --resources gpu=1 -j4
```
2. If this works run the pipeline
```
squeeze --resources gpu=1 -j4
```

1.

Infos
----

- Python Package and terminal: https://python-packaging.readthedocs.io/en/latest/command-line-scripts.html
- Github workflow pypi: https://github.com/pypa/packaging.python.org/blob/main/source/guides/github-actions-ci-cd-sample/publish-to-test-pypi.yml

Execute
----

```
python3 setup.py sdist && pip3 install --upgrade .
twine upload --verbose dist/squeezemd-0.1.5.tar.gz
username: __token__
pw: pyPi token
```

Sync results / data
```
rsync -avz --include '*/' --include '/pdb/***' --include '/config/***' --exclude '*' peter@gpu:/home/peter/MD/6_BD001-K10-Q9
rsync -avz --include '/pdb/***' --include '/config/***' peter@gpu:/home/peter/MD/6_BD001-K10-Q9 .
rsync -avz  peter@gpu:/home/peter/MD/6_BD001-K10-Q9/config .
```


Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage

TODO
--------

- Adjust MD according to ChatGPT
- Add mutations to Surface by introducting single parameters and derive location of last_frame. Do proper error handling
- Use jinja2 for template preparation
- Use pepfiles for simulation information

CHEATSHEET
------
Code version for execution
> snakemake -R `snakemake --list-code-changes`

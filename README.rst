=========
squeezeMD
=========


.. image:: https://img.shields.io/pypi/v/squeezemd.svg
        :target: https://pypi.python.org/pypi/squeezemd

.. image:: https://img.shields.io/travis/pruethemann/squeezemd.svg
        :target: https://travis-ci.com/pruethemann/squeezemd

.. image:: https://readthedocs.org/projects/squeezemd/badge/?version=latest
        :target: https://squeezemd.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status




Molecular Dynamics pipeline


* Free software: GNU General Public License v3
* Documentation: https://squeezemd.readthedocs.io.


Features
--------

* TODO

Installation
----

1. Download this github repository:
 > git clone ..
2. Install the squeezemd by executing
> `cd sequeezeMD \
conda env create -f squeezeMD.yml`
or use mamba
`conda install -c conda-forge mamba`
`mamba env create -f squeezeMD.yml`

1. Download this GitHub repository:
```git clone ..```

2. Install the squeezeMD by executing
```
cd squeezeMD \
conda env create -f squeezeMD.yml
```


Demo workflow
----

1. The workflow can be tested by performing the following commands:
```bash
squeeze --resources gpu=1 -j4 -n
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
```


Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage


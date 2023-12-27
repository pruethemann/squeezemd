#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [ ]

test_requirements = ['pytest>=3', ]

setup(
    author="Rüthemann Peter",
    author_email='peter.ruethemann@gmail.com',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Molecular Dynamics pipeline",
    entry_points={
        'console_scripts': [
            'squeezemd=squeezemd.cli:main',
        ],
    },
    install_requires=requirements,
    license="GNU General Public License v3",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    package_data={'': ['bin/7_interaction_csv.x', 'bin/7_interaction-analyzer.x']},
    data_files=[('', ['Snakefile'])],
    keywords='squeezemd',
    name='squeezemd',
    packages=find_packages(include=['squeezemd', 'squeezemd.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/pruethemann/squeezemd',
    version='0.1.4',
    zip_safe=False,
    scripts=['bin/squeeze','bin/1_mutation.py', 'bin/2_createleap.py', 'bin/2_SplitChains.py', 'bin/3_MD.py',
             'bin/4_center.py','bin/5_interactionFingerprint.py', 'bin/6_ExplorativeTrajectoryAnalysis.py', 'bin/run.py'],
)

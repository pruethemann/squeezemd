#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = []

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
    long_description_content_type="text/markdown",
    long_description=readme + '\n\n' + history,
    entry_points={
        'console_scripts': [
            'squeezemd=squeezemd.cli:main',
        ],
    },
    install_requires=requirements,
    license="GNU General Public License v3",

    include_package_data=True,
    data_files=[('', ['Snakefile', 'config/pymol_template.pml', 'config/RMSD.rst', 'config/RMSF.rst', 'config/interaction.rst', 'config/posco.rst','config/pymol_template.pml'])],
    keywords='squeezemd',
    name='squeezemd',
    packages=find_packages(include=['squeezemd', 'squeezemd.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/pruethemann/squeezemd',
    version='0.4.4',
    zip_safe=False,
    scripts=['bin/squeeze','bin/1_mutation.py', 'bin/2_MD.py', 'bin/3_ExplorativeTrajectoryAnalysis.py', 'bin/4_centerMDTraj.py', 'bin/5_PoseScoring.py',
             'bin/6_Ana_PoseScoring.py', 'bin/7_interactionFingerprint.py', 'bin/8_GlobalFinterprintAnalysis.py', 'bin/9_FreeEnergyStats.py',
             'bin/10_InteractionSurface.py' ,'bin/Helper.py', 'install/install_bins_mac.sh', 'install/install_bins_linux.sh','install/upgrade_package.sh'],
)

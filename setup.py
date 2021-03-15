#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2020 dejunlin <dejun.lin@gmail.com>
# Created at 2020-02-12 23:01 on dejunlin@threonine.gs.washington.edu
# Usage: setup.py
# Description: Compute HiCRep reproducibility stratum-corrected correlation score (SCCS).
# Reference: Genome Res. 2017 Nov;27(11):1939-1949. doi: 10.1101/gr.220640.117
# The algorithm first normalizes the input contact matrices by the total
# number of contacts and then for each chromosome: 1) mean-filter the input
# matrices with an input window size; 2) exclude common zero entries in
# the input matrices; 3) compute the SCC score. It doesn't have the
# procedure to bootstrap the window-size parameter
#
# Distributed under terms of the GNU General Public License v3.0.
import setuptools

setuptools.setup(
    name="hicrep",
    python_requires=">=3.7.6",
    version="0.2.4",
    description="Python implementation of HiCRep stratum-adjusted correlation coefficient of Hi-C data with sparse contact matrix support",
    long_description="see https://github.com/dejunlin/hicrep",
    long_description_content_type="text/markdown",
    url="https://github.com/dejunlin/hicrep.git",
    author="Dejun Lin",
    author_email="dejun.lin@gmail.com",
    license="GPLv3",
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
    ],
    packages=["hicrep"],
    include_package_data=True,
    install_requires=[
        "Deprecated",
        "numpy>=1.17.0",
        "scipy",
        "cooler",
        "pandas",
        "h5py",
    ],
    entry_points={"console_scripts": ["hicrep=hicrep:main"]},
    data_files = [("", ["LICENSE.txt"])]
)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Sat Feb 21 07:24:07 2026
"""
Pre-processor script to fetch files for tests
"""

import pooch

# data subdir in the notebook folder
DATASETS_DIR = u'data/'

alex1c = pooch.create(path=DATASETS_DIR, base_url='doi:10.5281/zenodo.20038738')
alex1c.load_registry_from_doi()
_ = alex1c.fetch("12d_New_30p_320mW_steer_3.hdf5")
_ = alex1c.fetch("0023uLRpitc_NTP_20dT_0.5GndCl.hdf5")
# _ = phdf5.fetch("dsdna_d7_d17_50_50_1.hdf5")

mphmm = pooch.create(path=DATASETS_DIR, base_url='doi:10.5281/zenodo.5902313')
mphmm.load_registry_from_doi()
hp3_files = ['HP3_TE150_SPC630.hdf5', 'HP3_TE200_SPC630.hdf5']
for fn in hp3_files:
    mphmm.fetch(fn)


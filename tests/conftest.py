#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 18 12:02:27 2026

@author: paul
"""
import pooch
import pytest
import os

from fretbursts import loader, bg


# data subdir in the notebook folder
DATASETS_DIR = u'data/'

# alex1c = pooch.create(path=DATASETS_DIR, base_url='doi:10.6084/m9.figshare.1019906.v26')
# alex1c.load_registry_from_doi()

# phdf5 = pooch.create(path=DATASETS_DIR, base_url='doi:10.6084/m9.figshare.1456362.v15')
# phdf5.load_registry_from_doi()

# mphmm = pooch.create(path=DATASETS_DIR, base_url='doi:10.5281/zenodo.5902313')
# mphmm.load_registry_from_doi()

fetchers = dict()

def fetch_if(name:str, doi:str)->str:
    fn = DATASETS_DIR + name
    if not os.path.exists(fn):
        if doi not in fetchers:
            fetch = pooch.create(path=DATASETS_DIR, base_url=doi)
            fetch.load_registry_from_doi()
            fetchers[doi] = fetch
        fetchers[doi].fetch(name)
    return fn


def _alex_process(d):
    loader.alex_apply_period(d)
    d.calc_bg(bg.exp_fit, time_s=30, tail_min_us=300)
    d.burst_search(L=10, m=10, F=7)


def load_dataset_1ch(process=True):
    fname = fetch_if("0023uLRpitc_NTP_20dT_0.5GndCl.hdf5", 'doi:10.5281/zenodo.20038738')
    d = loader.photon_hdf5(fname)
    if process:
        _alex_process(d)
    return d


def load_dataset_1ch_nsalex(process=True):
    fname = fetch_if("HP3_TE150_SPC630.hdf5", 'doi:10.5281/zenodo.5902313')
    d = loader.photon_hdf5(fname)
    if process:
        _alex_process(d)
    return d

@pytest.fixture
def dataset_1ch_file():
    return fetch_if("0023uLRpitc_NTP_20dT_0.5GndCl.hdf5", 'doi:10.5281/zenodo.20038738')

def load_dataset_8ch():
    fname = fetch_if("12d_New_30p_320mW_steer_3.hdf5", 'doi:10.5281/zenodo.20038738')
    d = loader.photon_hdf5(fname)
    d.calc_bg(bg.exp_fit, time_s=30, tail_min_us=300)
    d.burst_search(L=10, m=10, F=7)
    return d

@pytest.fixture
def fake_pax_file():
    return fetch_if("0023uLRpitc_NTP_20dT_0.5GndCl.hdf5", 'doi:10.5281/zenodo.20038738')

def load_fake_pax():
    fname = fetch_if("0023uLRpitc_NTP_20dT_0.5GndCl.hdf5", 'doi:10.5281/zenodo.20038738')
    d = loader.photon_hdf5(fname)
    d.add(ALEX=False, meas_type='PAX')
    loader.alex_apply_period(d)
    d.calc_bg(bg.exp_fit, time_s=30, tail_min_us='auto')
    d.burst_search(L=10, m=10, F=6, pax=True)
    return d

    
def load_dataset_grouped(process=True):
    fn = ['HP3_TE150_SPC630.hdf5', 'HP3_TE200_SPC630.hdf5']
    fn = [fetch_if(f, 'doi:10.5281/zenodo.5902313') for f in fn]
    d = loader.photon_hdf5(fn)
    if process:
        _alex_process(d)
    return d


@pytest.fixture(scope="module")
def data_8ch(request):
    d = load_dataset_8ch()
    return d

@pytest.fixture(scope="module")
def data_1ch(request):
    d = load_dataset_1ch()
    return d

@pytest.fixture(scope="module")
def data_1ch_nsalex(request):
    d = load_dataset_1ch_nsalex()
    return d

@pytest.fixture(scope="module")
def data_grouped(request):
    d = load_dataset_grouped()
    return d

@pytest.fixture(scope="module", params=[
                                    load_dataset_1ch,
                                    load_dataset_1ch_nsalex,
                                    load_dataset_8ch,
                                    load_dataset_grouped
                                    ])
def data(request):
    load_func = request.param
    d = load_func()
    return d

@pytest.fixture(scope="module", params=[
                                    load_dataset_8ch,
                                    load_dataset_grouped
                                    ])
def data_mch(request):
    load_func = request.param
    d = load_func()
    return d

@pytest.fixture(scope='module', params=[load_dataset_1ch,
                                        load_dataset_1ch_nsalex,
                                        load_dataset_grouped])
def data_alex(request):
    load_func = request.param
    d = load_func()
    return d

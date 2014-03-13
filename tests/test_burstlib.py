#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
Module containing automated unit tests for FRETBursts.

Running the tests require `py.test.
"""

import pytest

import numpy as np
from loaders import load_multispot8, load_usalex, usalex_apply_period
from fretbursts_path_def import data_dir
import background as bg
import burstlib as bl


@pytest.fixture(scope="module")
def data():
    #fn = "12d_New_30p_320mW_steer_3.dat"
    #dir_ = "2013-05-15/"
    #fname = data_dir+dir_+fn
    #BT = 0.038
    #gamma = 0.43
    #d = load_multispot8(fname=fname, BT=BT, gamma=gamma)
    
    fn = "0023uLRpitc_NTP_20dT_0.5GndCl.sm"
    dir_ = u'/Users/anto/Dropbox/notebooks/frebursts_notebooks/data/'
    fname = dir_ + fn
    d = load_usalex(fname=fname, BT=0.11, gamma=1.)
    d.add(det_donor_accept=(0, 1), alex_period=4000, 
          D_ON=(2850, 580), A_ON=(900, 2580))
    usalex_apply_period(d)

    d.calc_bg(bg.exp_fit, time_s=30, tail_min_us=300)
    d.burst_search_t(L=10, m=10, P=None, F=7, ph_sel='DA')
    return d


##
# Test functions
#
def test_b_functions(data):
    itstart, iwidth, inum_ph, iistart, iiend, itend = 0, 1, 2, 3, 4, 5    
    d = data
    for mb in d.mburst:
        assert (bl.b_start(mb) == mb[:, itstart]).all()
        assert (bl.b_end(mb) == mb[:, itend]).all()
        assert (bl.b_width(mb) == mb[:, iwidth]).all()
        assert (bl.b_istart(mb) == mb[:, iistart]).all()
        assert (bl.b_iend(mb) == mb[:, iiend]).all()
        assert (bl.b_size(mb) == mb[:, inum_ph]).all()

        rate = 1.*mb[:, inum_ph]/mb[:, iwidth]        
        assert (bl.b_rate(mb) == rate).all()

        separation = mb[1:, itstart] - mb[:-1, itend]     
        assert (bl.b_separation(mb) == separation).all()
        
        assert (bl.b_end(mb) > bl.b_start(mb)).all()

        
def test_b_end_b_iend(data):
    """Test coherence between b_end() and b_iend()"""
    d = data    
    for ph, mb in zip(d.ph_times_m, d.mburst):
        assert (ph[bl.b_iend(mb)] == bl.b_end(mb)).all()

def test_monotonic_burst_start(data):
    """Test for monotonic burst_start."""
    d = data
    for i in xrange(d.nch):
        assert (np.diff(bl.b_start(d.mburst[i])) > 0).all()

def test_monotonic_burst_end(data):
    """Test for monotonic burst_end."""
    d = data
    for mb in d.mburst:
        assert (np.diff(bl.b_end(mb)) > 0).all()

def test_burst_start_end_size(data):
    """Test consistency between burst istart, iend and size"""
    d = data
    for mb in d.mburst:
        size = mb[:, bl.iiend] - mb[:, bl.iistart] + 1
        assert (size == mb[:, bl.inum_ph]).all()
        size2 = bl.b_iend(mb) - bl.b_istart(mb) + 1
        assert (size2 == bl.b_size(mb)).all()

def test_burst_fuse_0ms(data):
    """Test that after fusing with ms=0 the sum of bursts sizes is that same
    as the number of ph in bursts (via burst selection).
    """
    d = data
    if not hasattr(d, 'fuse'):
        df = d.fuse_bursts(ms=0)
        for ph, mb in zip(df.ph_times_m, df.mburst):
            m = bl.ph_select(ph, mb)
            assert m.sum() == bl.b_size(mb).sum()
#            
if __name__ == '__main__':
    pytest.main("-x -v tests/test_burstlib.py")
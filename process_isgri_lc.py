from __future__ import print_function
from copy import deepcopy
import dataanalysis.callback

import ddosa
from astropy.io import fits
try:
    from pscolors import render
except:
    from bcolors import render
import subprocess
import os

import dataanalysis as da
from dataanalysis import graphtools

from numpy import *
import numpy as np
from collections import defaultdict

from astropy.io.fits import Column

try:
    import crab
except:
    pass

# Importing the library
import psutil
  
import logging

logging.basicConfig(level='INFO')

logger = logging.getLogger()



import linecache
import os
import tracemalloc

def display_top(snapshot, key_type='lineno', limit=3):
    snapshot = snapshot.filter_traces((
        tracemalloc.Filter(False, "<frozen importlib._bootstrap>"),
        tracemalloc.Filter(False, "<unknown>"),
    ))
    top_stats = snapshot.statistics(key_type)

    def cprint(*args):
        print("\033[34m", *args, "\033[0m")

    cprint("Top %s lines" % limit)
    for index, stat in enumerate(top_stats[:limit], 1):
        frame = stat.traceback[0]
        # replace "/path/to/module/file.py" with "module/file.py"
        filename = os.sep.join(frame.filename.split(os.sep)[-2:])
        cprint("#%s: %s:%s: %.1f KiB"
              % (index, filename, frame.lineno, stat.size / 1024))
        line = linecache.getline(frame.filename, frame.lineno).strip()
        if line:
            cprint('    %s' % line)

    other = top_stats[limit:]
    if other:
        size = sum(stat.size for stat in other)
        cprint("%s other: %.1f KiB" % (len(other), size / 1024))
    total = sum(stat.size for stat in top_stats)
    cprint("Total allocated size: \033[31m%.1f MiB\033[0m" % (total / 1024 / 1024))



def get_open_fds():
    '''
    return the number of open file descriptors for current process

    .. warning: will only work on UNIX-like os-es.
    '''

    pid = os.getpid()
    procs = subprocess.check_output(
        # [ "lsof", '-w', '-Ff', "-p", str( pid ) ] )
        ["lsof", '-w', '-Fn', "-p", str(pid)])

    # files=filter(
    #        lambda s: s and s[ 0 ] == 'f' and s[1: ].isdigit(),
    files = procs.decode().split('\n')

    # print(files)

    nprocs = len(
        files
    )
    return nprocs


class FactorLCfromScW(graphtools.Factorize):
    root = 'ii_lc_extract'
    leaves = ["ScWData", "Revolution"]


class ScWLCList(ddosa.DataAnalysis):
    input_scwlist = ddosa.RevScWList
    copy_cached_input = False
    input_lcsummary = FactorLCfromScW

    allow_alias = True

    version = "allthem"

    def main(self):
        self.lcs = [[ddosa.ii_lc_extract(assume=scw)]
                    for scw in self.input_scwlist.scwlistdata]

        if len(self.lcs) == 0:
            raise ddosa.EmptyScWList()


class ISGRILCSum(ddosa.DataAnalysis):
    input_lclist = ScWLCList

    copy_cached_input = False

    cached = True

    version = "v1.1.7"

    sources = ["Crab"]
    extract_all = True
    # save_spec=True

    def get_version(self):
        v = self.get_signature()+"."+self.version
        if self.extract_all:
            v += ".extractall"
        else:
            v+".extract_"+("_".join(self.sources))
        return v

    def patch_isgri_lc_xax_e(self, lc):
        timedel = lc.header['TIMEDEL']
        time = lc.data['TIME']

        n_lines = len(lc.data)

        new_time = zeros(n_lines)
        new_time[0] = timedel/2
    
        for i in range(1, n_lines):
            x = time[i] - time[i-1] - new_time[i-1]
            if x <= 0:
                logger.warning(
                    "lc_pick cleaning ISGRI, inconsistent binning at index %d", i)
                continue

            if x < timedel/2:
                new_time[i] = x
            else:
                new_time[i] = timedel/2


        nc = lc.columns.add_col(fits.Column('XAX_E', format='1E', array=new_time))
        nlc = fits.BinTableHDU.from_columns(nc)

        for k,v in lc.header.items():
            if k not in nlc.header:
                nlc.header[k] = v

        return nlc
		
    def main(self):
        tracemalloc.start()

        lcs = {}

        choice = self.input_lclist.lcs

        allsource_summary = []

        def sig(x, y): return (
            ((x/y)[~np.isnan(y) & (y != 0) & ~np.isinf(y) & ~np.isinf(x) & ~np.isnan(x)])**2).sum()**0.5

        import time
        t0 = time.time()
        i_lc = 1

        used_fns = []

        for lc, in choice:
            if hasattr(lc, 'empty_results'):
                print("skipping, for clearly empty:", lc)
                continue

            if not hasattr(lc, 'lightcurve'):
                print("skipping, for whatever reason no data:", lc)
                print(dir(lc))
                continue

            snapshot_pre_loop = tracemalloc.take_snapshot()

            fn = lc.lightcurve.get_path()

            if fn in used_fns:
                raise RuntimeError(f'found duplicate LC file {fn}, used so far {used_fns}')
            else:
                used_fns.append(fn)

            print("%i/%i" % (i_lc, len(choice)))
            tc = time.time()
            print("seconds per lc:", (tc-t0)/i_lc, "will be ready in %.5lg seconds" %
                  ((len(choice)-i_lc)*(tc-t0)/i_lc))
            i_lc += 1
            print("lc from", fn)
                        
            f = fits.open(fn, memmap=False)
            print("proceeding to parse", f.filename())

            t1, t2 = f[1].header['TSTART'], f[1].header['TSTOP']
            print(t1, t2)

            for _e in f:                    
                e = deepcopy(_e.copy())
                del _e

                if e.header.get('EXTNAME', 'unnamed') != "ISGR-SRC.-LCR":
                    continue
                
                name = e.header.get('NAME', "Unnamed")

                allsource_summary.append(
                    [name, t1, t2, e.data['RATE'], e.data['ERROR']])

                if (name in self.sources) or (self.extract_all):
                    rate = e.data['RATE']
                    err = e.data['ERROR']
                    if name not in lcs:
                        print("new lcs[name]", name)
                        lcs[name] = e
                    else:                        
                        print("lcs[name].data of", len(lcs[name].data), "e.data of", len(e.data))
                        lcs[name].data = concatenate((lcs[name].data, e.data))

                    print(render("{BLUE}%.20s{/}" % name), "%.4lg sigma" % (sig(rate, err)),
                        "total %.4lg" % (sig(lcs[name].data['RATE'], lcs[name].data['ERROR'])))

                    print("\033[31msize of lcs[name]", lcs[name].size/1024/1024, "Mb" , lcs[name].data.size * lcs[name].data.itemsize/1024/1024, "Mb\033[0m")

            del f


            print('RAM memory % used:', psutil.virtual_memory()[2])
            print('RAM memory:', psutil.virtual_memory())

            try:
                print("get_open_fds", get_open_fds())
            except Exception as e:
                print("unable to check open fds", e)

            snapshot = tracemalloc.take_snapshot()
            display_top(snapshot)


            diff = snapshot_pre_loop.compare_to(snapshot, 'lineno')
            
            print("[ Top 10 differences ]")
            for stat in diff[:10]:
                print("\033[32m", stat, "\033[0m")
            # assert diff[0].size_diff*u.B < 1*u.MB

            # pick the biggest memory block
            stat = diff[0]
            print("%s memory blocks: %.1f KiB" % (stat.count, stat.size / 1024))
            for line in stat.traceback.format():
                print(line)



        # self.lcs=lcs

        # source_results = []
        self.extracted_sources = []

        for name, lc in lcs.items():
            source_short_name = name.strip().replace(" ", "_")

            prepatch_fn = "isgri_sum_lc_prepatch_%s.fits" % source_short_name
            fn = "isgri_sum_lc_%s.fits" % source_short_name

            lc.writeto(prepatch_fn, clobber=True)

            reopened_lc = fits.open(prepatch_fn)
            patched_lc = fits.HDUList([fits.PrimaryHDU(), self.patch_isgri_lc_xax_e(reopened_lc[1])])
            patched_lc.writeto(fn, overwrite=True)

            reopened_lc.close()

            # test
            test_lc = fits.open(fn)
            print("patched lc columns:", test_lc[1].data.columns)
            assert 'XAX_E' in [ c.name for c in test_lc[1].data.columns ]
            test_lc.close()

            attr = fn.replace(".fits", "")
            self.extracted_sources.append([name, attr])

            setattr(self, attr, da.DataFile(fn))

            snapshot = tracemalloc.take_snapshot()
            display_top(snapshot)


dataanalysis.callback.default_callback_filter.set_callback_accepted_classes(
    [ddosa.mosaic_ii_skyimage, ddosa.ii_skyimage, ddosa.BinEventsImage, ddosa.ibis_gti, ddosa.ibis_dead, ddosa.ISGRIEvents, ddosa.ii_spectra_extract, ddosa.BinEventsSpectra, ddosa.ii_lc_extract, ddosa.BinEventsLC, ISGRILCSum])

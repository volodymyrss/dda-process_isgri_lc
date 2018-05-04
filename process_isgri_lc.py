from __future__ import print_function

import ddosa 
from astropy.io import fits 
from bcolors import render
import subprocess
import os

import dataanalysis as da
from dataanalysis import graphtools

from numpy import *
from collections import defaultdict

try:
    import crab
except:
    pass

def get_open_fds():
    '''
    return the number of open file descriptors for current process

    .. warning: will only work on UNIX-like os-es.
    '''

    pid = os.getpid()
    procs = subprocess.check_output( 
        #[ "lsof", '-w', '-Ff', "-p", str( pid ) ] )
        [ "lsof", '-w', '-Fn', "-p", str( pid ) ] )

    #files=filter( 
    #        lambda s: s and s[ 0 ] == 'f' and s[1: ].isdigit(),
    files=        procs.split( '\n' ) 

    #print(files)

    nprocs = len( 
            files
        )
    return nprocs


class FactorLCfromScW(graphtools.Factorize):
    root='ii_lc_extract'
    leaves=["ScWData","Revolution"]

class ScWLCList(ddosa.DataAnalysis):
    input_scwlist=ddosa.RevScWList
    copy_cached_input=False
    input_lcsummary=FactorLCfromScW

    allow_alias=True

    version="allthem"
    
    def main(self):
        self.lcs=[[ddosa.ii_lc_extract(assume=scw)] for scw in self.input_scwlist.scwlistdata]

        if len(self.lcs)==0:
            raise ddosa.EmptyScWList()
        

class ISGRILCSum(ddosa.DataAnalysis):
    input_lclist=ScWLCList

    copy_cached_input=False


    cached=True

    version="v1"

    sources=["Crab"]
    extract_all=True
    #save_spec=True

    def get_version(self):
        v=self.get_signature()+"."+self.version
        if self.extract_all:
            v+=".extractall"
        else:
            v+".extract_"+("_".join(self.sources))
        return v

    def main(self):
        lcs={}

        choice=self.input_lclist.lcs

        allsource_summary=[]


        sig=lambda x,y:(((x/y)[~isnan(y) & (y!=0) & ~isinf(y) & ~isinf(x) & ~isnan(x)])**2).sum()**0.5

        import time
        t0=time.time()
        i_lc=1

        for lc, in choice:
            if hasattr(lc,'empty_results'):
                print("skipping, for clearly empty:",lc)
                continue

            if not hasattr(lc,'lightcurve'):
                print("skipping, for whatever reason no data:",lc)
                print(dir(lc))
                continue


            fn=lc.lightcurve.get_path()
            print("%i/%i"%(i_lc,len(choice)))
            tc=time.time()
            print("seconds per lc:",(tc-t0)/i_lc,"will be ready in %.5lg seconds"%((len(choice)-i_lc)*(tc-t0)/i_lc))
            i_lc+=1
            print("lc from",fn)

            f=fits.open(fn)

            t1,t2=f[1].header['TSTART'],f[1].header['TSTOP']
            print(t1,t2)

            for e in f:
                try:
                    if e.header['EXTNAME']!="ISGR-SRC.-LCR": continue
                except:
                    continue

                try:
                    name=e.header['NAME']
                except:
                    name="Unnamed"

                allsource_summary.append([name,t1,t2,copy(e.data['RATE']),copy(e.data['ERROR'])])
                if (name in self.sources) or (self.extract_all):
                    rate=e.data['RATE']
                    err=e.data['ERROR']
                    exposure=e.header['EXPOSURE']
                    if name not in lcs:
                        lcs[name]=e
                        preserve_file=True
                    else:
                        lcs[name].data=concatenate((lcs[name].data,e.data))

                    print(render("{BLUE}%.20s{/}"%name),"%.4lg sigma"%(sig(rate,err)),"total %.4lg"%(sig(lcs[name].data['RATE'],lcs[name].data['ERROR'])))

            f.close()

            try:
                print(get_open_fds())
            except Exception as e:
                print("unable to check open fds")

        #self.lcs=lcs

        source_results=[]
        self.extracted_sources=[]

        for name,lc in lcs.items():
            source_short_name=name.strip().replace(" ","_")


            fn="isgri_sum_lc_%s.fits"%source_short_name
            lc.writeto(fn,clobber=True)

            setattr(self,fn.replace(".fits",""),da.DataFile(fn))





import dataanalysis.callback

dataanalysis.callback.default_callback_filter.set_callback_accepted_classes([ddosa.mosaic_ii_skyimage, ddosa.ii_skyimage, ddosa.BinEventsImage, ddosa.ibis_gti, ddosa.ibis_dead, ddosa.ISGRIEvents, ddosa.ii_spectra_extract, ddosa.BinEventsSpectra, ddosa.ii_lc_extract, ddosa.BinEventsLC, ISGRILCSum])



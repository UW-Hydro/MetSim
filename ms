#!/usr/bin/env python3

"""
Command line tool for MetSim
"""

import sys
import pprint
import multiprocessing

def main():
    import metsim
    from metsim import io
    from metsim import options
    from metsim import methods
    from metsim.metsim import MetSim 
    options.init(options.parse(sys.argv[1:])) 
  
    ms = MetSim(methods.mtclim, metsim.forcing_files, metsim.proc_count)
    ms.run()

if __name__ == '__main__':
    main()
    

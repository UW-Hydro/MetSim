#!/usr/bin/env python3

"""
Command line tool for MetSim
"""

import sys
import pprint
import multiprocessing

def main():
    import metsim
    from metsim import options
    from metsim import methods
    from metsim import io
    options.init(options.parse(sys.argv[1:])) 
  
    metsim.manager = multiprocessing.Manager()
    metsim.write_queue = metsim.manager.Queue()
    process_handles = [
            multiprocessing.Process(target=metsim.method, args=(flist, metsim.write_queue))
            for flist in metsim.forcing_chunks
            ]
    for p in process_handles:
        p.start()
    for p in process_handles:
        p.join()


if __name__ == '__main__':
    main()
    

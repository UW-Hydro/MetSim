#!/usr/bin/env python3

"""
Command line tool for MetSim
"""

import sys
import multiprocessing

def main():
    import metsim
    from metsim import options
    from metsim import methods
    from metsim import io
    options.init(options.parse(sys.argv[1:])) 
    io.read_config(metsim.config)
    metsim.manager = multiprocessing.Manager()
    metsim.write_queue = manager.Queue()
    process_handles = [] #TODO
    for p in process_handles:
        p.start()
    for p in process_handles:
        p.join()


if __name__ == '__main__':
    main()
    

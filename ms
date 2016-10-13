#!/usr/bin/env python3

"""
Command line tool for MetSim
"""

import sys

def main():
    import metsim
    from metsim import io
    from metsim import options
    from metsim.metsim import MetSim 
    options.init(options.parse(sys.argv[1:])) 
 
    ms = MetSim(metsim.method, metsim.forcing_files, metsim.proc_count)
    ms.launch_processes()

if __name__ == '__main__':
    main()
    

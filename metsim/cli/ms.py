#!/usr/bin/env python3
"""
Command line tool for MetSim
"""
# Meteorology Simulator
# Copyright (C) 2017  The Computational Hydrology Group, Department of Civil
# and Environmental Engineering, University of Washington.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import argparse
import json
import logging
import os
import sys
from collections import OrderedDict
from configparser import SafeConfigParser


def _is_valid_file(parser, arg):
    if not os.path.isfile(arg):
        message = "Configuration file {} does not exist.\n".format(arg)
        parser.print_help()
        parser.exit(status=1, message=message)
    return arg


def parse(args):
    """Parse the command line arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument('config', type=lambda x: _is_valid_file(parser, x),
                        help='Input configuration file')
    parser.add_argument('-n', '--n-processes', default=1, type=int,
                        help='Parallel mode: number of processes to use')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Increase the verbosity of MetSim')
    parser.add_argument('-tg', '--time_grouper', nargs='?',
                        const='1AS', default=None, type=str,
                        help='Pandas TimeGrouper string (e.g. `1AS`)')
    return parser.parse_args()


def init(opts):
    """Initialize some information based on the options & config"""
    config = SafeConfigParser()
    config.optionxform = str
    config.read(opts.config)
    conf = OrderedDict(config['MetSim'])
    conf['forcing_vars'] = OrderedDict(config['forcing_vars'])
    conf['domain_vars'] = OrderedDict(config['domain_vars'])
    conf['state_vars'] = OrderedDict(config['state_vars'])
    out_dir = conf['out_dir']
    out_state = conf.get('out_state', None)
    if out_state is None:
        out_state = os.path.join(out_dir, 'state.nc')

    method = conf['method']

    # If the forcing variable is a directory, scan it for files
    if os.path.isdir(conf['forcing']):
        forcing_files = [os.path.join(conf['forcing'], fn) for fn in
                         next(os.walk(conf['forcing']))[2]]
    else:
        forcing_files = conf['forcing']

    # We assume there is only one domain file and one state file
    domain_file = conf['domain']
    state_file = conf['state']

    def to_list(s):
        return json.loads(s.replace("'", '"'))

    conf.update({"calendar": conf.get('calendar', 'standard'),
                 "nprocs": opts.n_processes,
                 "method": method,
                 "out_dir": out_dir,
                 "out_state": out_state,
                 "state": state_file,
                 "domain": domain_file,
                 "forcing": forcing_files,
                 "verbose": opts.verbose * logging.INFO})
    conf['out_vars'] = to_list(conf.get('out_vars', '[]'))
    conf['iter_dims'] = to_list(conf.get('iter_dims', '["lat", "lon"]'))
    if opts.time_grouper is not None:
        conf['time_grouper'] = opts.time_grouper
    conf = {k: v for k, v in conf.items() if v != []}
    return conf


def main():
    """Runs MetSim"""
    from metsim.metsim import MetSim
    setup = init(parse(sys.argv[1:]))
    ms = MetSim(setup)
    if ms.params['nprocs'] > 1:
        ms.launch()
    else:
        ms.run()


if __name__ == '__main__':
    main()

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
from configparser import ConfigParser


def _is_valid_file(parser, arg):
    if not os.path.isfile(arg):
        message = "Configuration file {} does not exist.\n".format(arg)
        parser.print_help()
        parser.exit(status=1, message=message)
    return arg


def parse(args):
    """Parse the command line arguments"""
    from metsim import __name__, __version__
    parser = argparse.ArgumentParser()
    parser.add_argument('config', type=lambda x: _is_valid_file(parser, x),
                        help='Input configuration file')
    parser.add_argument('-n', '--num_workers', default=1, type=int,
                        help='Parallel mode: number of processes to use')
    parser.add_argument('-s', '--scheduler', default='distributed', type=str,
                        help='Dask scheduler to use')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Increase the verbosity of MetSim')
    parser.add_argument('--version', action='version',
                        version='{} {}'.format(__name__, __version__),
                        help='Name and version number')
    return parser.parse_args()


def init(opts):
    """Initialize some information based on the options & config"""
    config = ConfigParser()
    config.optionxform = str
    config.read(opts.config)
    conf = OrderedDict(config['MetSim'])

    def invert_dict(d):
        return OrderedDict([reversed(item) for item in d.items()])

    def to_list(s):
        return json.loads(s.replace("'", '"').split('#')[0])

    conf['forcing_vars'] = OrderedDict(config['forcing_vars'])
    if conf['forcing_fmt'] != 'binary':
        conf['forcing_vars'] = invert_dict(conf['forcing_vars'])
    conf['domain_vars'] = invert_dict(OrderedDict(config['domain_vars']))
    conf['state_vars'] = invert_dict(OrderedDict(config['state_vars']))
    conf['chunks'] = OrderedDict(config['chunks'])
    if 'constant_vars' in config:
        conf['constant_vars'] = OrderedDict(config['constant_vars'])

    # If the forcing variable is a directory, scan it for files
    if os.path.isdir(conf['forcing']):
        forcing_files = [os.path.join(conf['forcing'], fn) for fn in
                         next(os.walk(conf['forcing']))[2]]
    else:
        forcing_files = conf['forcing']

    # Update the full configuration
    conf.update({"calendar": conf.get('calendar', 'standard'),
                 "scheduler": opts.scheduler,
                 "num_workers": opts.num_workers,
                 "verbose": logging.DEBUG if opts.verbose else logging.INFO,
                 "forcing": forcing_files,
                 "out_dir": os.path.abspath(conf['out_dir']),
                 "prec_type": conf.get('prec_type', 'uniform')})
    conf['out_vars'] = to_list(conf.get('out_vars', '[]'))

    conf = {k: v for k, v in conf.items() if v != []}
    return conf


def main():
    """Runs MetSim"""
    from metsim.metsim import MetSim
    setup = init(parse(sys.argv[1:]))
    ms = MetSim(setup)
    ms.run()


if __name__ == '__main__':
    main()

#!/usr/bin/env python

import xarray as xr
from metsim import MetSim
import os.path
import sys, getopt
import pandas as pd
from collections import OrderedDict
import datetime
from configparser import SafeConfigParser
import json

def str_to_bool(s):
    if s == 'True' or s == '1':
         return True
    elif s == 'False' or s == '0':
         return False
    else:
         raise ValueError('spinup must be set to True, 1, False, or 0')

def init(config_file):
    """Initialize some information based on the options & config"""
    if not os.path.isfile(config_file):
        exit("Configuration file {} does not exist.".format(config_file)
             + "Use `python wrap_metsim_over_time.py -h` "
               "for more information.")
    config = SafeConfigParser()
    config.optionxform = str
    config.read(config_file)
    conf = OrderedDict(config['MetSim'])
    conf['forcing_vars'] = OrderedDict(config['forcing_vars'])
    conf['domain_vars'] = OrderedDict(config['domain_vars'])
    conf['state_vars'] = OrderedDict(config['state_vars'])

    method = conf['method']
    startdate = conf.get('start')
    stopdate = conf.get('stop')
    forcing_fmt = conf.get('forcing_fmt',None)
    domain_fmt = conf.get('domain_fmt',None)
    state_fmt = conf.get('state_fmt',None)
    in_fmt = conf.get('in_fmt',None)
    out_fmt = conf.get('out_fmt',None)
    out_dir = conf['out_dir']
    time_step = conf.get('time_step',None)
    out_prefix = conf.get('out_prefix', None)
    out_state_tmp = out_dir + '/' + out_prefix + '.state_tmp.nc'
    out_state = conf['out_state']
    spinup = str_to_bool(conf['spinup'])
    in_dir = conf['in_dir']
    in_prefix = conf['in_prefix']

    prec_type = conf.get('prec_type', None)
    if prec_type is None:
        prec_type = 'uniform'

    domain = conf['domain']
    state_initial = conf['state']

    def to_list(s):
        return json.loads(s.replace("'", '"'))

    conf.update({"calendar": conf.get('calendar', 'standard'),
                 "method": method,
                 "start": startdate,
                 "stop": stopdate,
                 "out_dir": out_dir,
                 "out_state": out_state_tmp,
                 "out_prefix": out_prefix + '.forcing_tmp',
                 "state": state_initial,
                 "domain": domain,
                 "prec_type": prec_type,
                 "time_step": time_step,
                 "forcing_fmt": forcing_fmt,
                 "domain_fmt": domain_fmt,
                 "state_fmt": state_fmt,
                 "out_fmt": out_fmt,
                 "annual": False})

    conf['out_vars'] = to_list(conf.get('out_vars', '[]'))
    conf['iter_dims'] = to_list(conf.get('iter_dims', '["lat", "lon"]'))
    conf = {k: v for k, v in conf.items() if v != []}
    return conf, startdate, stopdate, out_dir, out_prefix, spinup, in_dir, \
           in_prefix, state_initial, out_state, out_state_tmp

def main():
    computation_start = datetime.datetime.now()
    config_file = ''

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hc:",
                                                 ["config_file="])
    except getopt.GetoptError:
        print(sys.argv[0], '-c <config_file>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], '-c <config_file>')
            sys.exit()
        elif opt in ("-c", "--config_file"):
            config_file = arg

    params, startdate, stopdate, out_dir, out_prefix, spinup, in_dir, \
    in_prefix, state_initial, out_state, out_state_tmp = init(config_file)

    try:
        startdate, junk = startdate.split(':')
    except:
        pass
    try:
        stopdate, junk = stopdate.split(':')
    except:
        pass
    dates = pd.date_range(start=startdate, end=stopdate, freq='MS')
    first_year, first_month, first_day = startdate.split('/')
    first_month = int(first_month)
    first_day = int(first_day)
    first_year = int(first_year)
    last_year, last_month, last_day = stopdate.split('/')
    last_month = int(last_month)
    last_day = int(last_day)
    last_year = int(last_year)

    forcing_file = in_dir + '/' + \
                   in_prefix + ".%04d%02d" % (first_year, first_month) + \
                   '.nc'

    params['forcing'] = forcing_file
    year = first_year
    date_i = 0

    forcing_outfile_tmp = out_dir + '/' + out_prefix + '.forcing_tmp_total.nc'

    if spinup:
        print('Doing spinup year')
    while year <= last_year:
        for month_i in range(12):
            startdate = dates[date_i]
            month_current = startdate.month
            year_current = startdate.year
            if month_current < 12:
                stopdate = pd.to_datetime(("%04d%02d%02d" % (year_current,
                                                            month_current +
                                                            1, 1)), format='%Y%m%d') - pd.Timedelta('1 days')
            else:
                stopdate = pd.to_datetime(("%04d%02d%02d" % (year_current +
                                                            1, 1, 1)), format='%Y%m%d') - pd.Timedelta('1 days')

            if year == first_year and month_i == 0:

                # delete any previous temporary state files created with
                # ms.run()
                try:
                    os.remove(out_state_tmp)
                except:
                    pass
                params["forcing"] = forcing_file
                params["start"] = startdate
                params["stop"] = stopdate
                print('Running forcings from start and stop dates:', startdate,
                      stopdate)
                ms = MetSim(params)
                ms.run()

                # delete any previous temporary state files created in wrapper
                # script and save the current
                try:
                    os.remove(out_state)
                except:
                    pass
                ms.state.to_netcdf(out_state)
            else:
                forcing_filename = (in_prefix +
                                    ".%04d%02d" % (year_current,
                                                   month_current) + '.nc')
                infile = (in_dir + '/' + forcing_filename)

                # delete any previous temporary state files created with
                # ms.run()
                try:
                    os.remove(out_state_tmp)
                except:
                    pass

                ms = []
                params['state'] = out_state
                params['start'] = startdate
                params['stop'] = stopdate
                params['forcing'] = infile
                ms = MetSim(params)
                print('Running forcings from start and stop dates:', startdate,
                      stopdate)
                ms.run()

                # delete any previous temporary state files created in wrapper
                # script and save the current
                try:
                    os.remove(out_state)
                except:
                    pass
                ms.state.to_netcdf(out_state)

            if not spinup:

                print('Concatenating monthly dataset to the yearly dataset')

                # concat the month's dataset to the yearly dataset
                forcing_current = xr.open_dataset(forcing_outfile_tmp)
                if month_i == 0:
                    forcing_concat = forcing_current
                else:
                    forcing_concat = xr.concat([forcing_concat,
                                                forcing_current], dim='time')
                #forcing_current.close()
            # delete the current month's forcing files
            try:
                os.remove(forcing_outfile_tmp)
            except:
                pass

            date_i += 1
            print('Completed forcings from start and stop dates:',
                  startdate, stopdate)
            if date_i > len(dates) - 1:
                break

        if spinup:
            # save the state file at the end of spinup
            print('Saving state file at end of spinup periodnc')
            ds = xr.open_dataset(state_initial)
            ms.state.coords['time'] = ds.coords['time']
            ds.close()
            state_spinup_end_file = out_dir + '/' + out_prefix + \
                                    '.state_spinup_end.nc'
            try:
                os.remove(state_spinup_end_file)
            except:
                pass
            ms.state.to_netcdf(state_spinup_end_file)
            spinup = 0
            date_i = 0
            params['state'] = state_spinup_end_file

        else:
            # output concatenated dataset
            outfile = (out_dir + '/' + out_prefix + ".%04d" % (year) + '.nc')
            print('Saving yearly datasets to outfile:', outfile)
            try:
                os.remove(outfile)
            except:
                pass
            forcing_concat.to_netcdf(outfile, engine='scipy')
            year = year + 1

    # delete any previous temporary state files created with ms.run()
    try:
        os.remove(out_dir + '/' + out_prefix + '.state_tmp.nc')
    except:
        pass

    # delete the last month's forcing files
    try:
        os.remove(forcing_outfile_tmp)
    except:
        pass
    computation_end = datetime.datetime.now()
    computation_total = computation_end - computation_start
    print('Total computation time:', computation_total)
    forcing_current.close()

if __name__ == "__main__":
    main()

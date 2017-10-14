#!/usr/bin/env python

import xarray as xr
from metsim import MetSim
import os.path
import sys, getopt
import pandas as pd
from collections import OrderedDict
import datetime


def main():
    computation_start = datetime.datetime.now()
    infile_directory = ''
    spinup = ''
    domainfile = ''
    statefile_initial = ''
    infile_prefix = ''
    outfile_prefix = ''
    outfile_directory = ''
    startdate = ''
    enddate = ''
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:p:d:s:f:o:u:a:b:",
                                                 ["infile_directory=", "spinup="
                                                  "domainfile=",
                                                  "statefile_initial=",
                                                  "infile_prefix=",
                                                  "outfile_prefix=",
                                                  "outfile_directory=",
                                                  "startdate",
                                                  "enddate="])
    except getopt.GetoptError:
        print(sys.argv[0], '-i <infile_directory> -p <spinup>'
                           '-d <domainfile> -f <infile_prefix> '
                           '-s <statefile_initial> -o <outfile_prefix> -u '
                           '<outfile_directory> -a <startdate> -b <enddate>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], ' -i <infile_directory> -p <spinup> '
                               '-d <domainfile> -s <statefile_initial> '
                               '-f <infile_prefix> -o <outfile_prefix> -u '
                               '<outfile_directory> -a <startdate>'
                               ' -b <enddate>')
            sys.exit()
        elif opt in ("-i", "--infile_directory"):
            infile_directory = arg
        elif opt in ("-p", "--spinup"):
            spinup = int(arg)
        elif opt in ("-d", "--domainfile"):
            domainfile = arg
        elif opt in ("-s", "--statefile_initial"):
            statefile_initial = arg
        elif opt in ("-f", "--infile_prefix"):
            infile_prefix = arg
        elif opt in ("-o", "--outfile_prefix"):
            outfile_prefix = arg
        elif opt in ("-u", "--outfile_directory"):
            outfile_directory = arg
        elif opt in ("-a", "--startdate"):
            startdate = arg
        elif opt in ("-b", "--enddate"):
            enddate = arg

    # infiles = os.listdir(infile_directory)
    outfile_directory = outfile_directory + '/'
    state = xr.open_dataset(statefile_initial)
    domain = xr.open_dataset(domainfile)

    dates = pd.date_range(start=startdate, end=enddate, freq='MS')
    first_year, first_month, first_day = startdate.split('-')
    first_month = int(first_month)
    first_day = int(first_day)
    first_year = int(first_year)
    last_year, last_month, last_day = enddate.split('-')
    last_month = int(last_month)
    last_day = int(last_day)
    last_year = int(last_year)

    forcing_file = infile_directory + '/' + \
                   infile_prefix + "%04d%02d" % (first_year, first_month) + \
                   '.nc'
    forcing_intial = xr.open_dataset(forcing_file)
    outfile_prefix_month_1 = "%s%04d%02d" % (outfile_prefix, first_year, first_month)

    year = first_year
    date_i = 0

    params = {
        "method": 'mtclim',
        "domain": domainfile,
        "state": statefile_initial,
        "forcing_fmt": 'netcdf',
        "domain_fmt": 'netcdf',
        "state_fmt": 'netcdf',
        "out_dir": outfile_directory,
        "out_state": outfile_directory + 'state_tmp.nc',
        "out_prefix": 'forcing_tmp',
        "time_step": '60',
        "out_fmt": 'netcdf',
        "annual": False,
        'forcing': forcing_file,
        'state_vars': OrderedDict({'prec': 'prec',
                                   't_max': 't_max',
                                   't_min': 't_min',
                                   'swe': 'swe'}),
        'forcing_vars': OrderedDict({'Prec': 'prec',
                                     'Tmax': 't_max',
                                     'Tmin': 't_min',
                                     'wind': 'wind'}),
        'domain_vars': OrderedDict({'lat': 'lat',
                                    'lon': 'lon',
                                    'mask': 'mask',
                                    'elev': 'elev',
                                    't_pk': 't_pk',
                                    'dur': 'dur'})
    }

    forcing_outfile_tmp = outfile_directory + 'forcing_tmp' + '_total.nc'
    state_outfile = outfile_directory + 'state.nc'

    if spinup:
        print('Doing spinup year')
    while year <= last_year:
        for month_i in range(12):
            startdate = dates[date_i]
            month_current = startdate.month
            year_current = startdate.year
            if month_current < 12:
                enddate = pd.to_datetime(("%04d%02d%02d" % (year_current,
                                                            month_current +
                                                            1, 1)), format='%Y%m%d') - pd.Timedelta('1 days')
            else:
                enddate = pd.to_datetime(("%04d%02d%02d" % (year_current +
                                                            1, 1, 1)), format='%Y%m%d') - pd.Timedelta('1 days')

            if year == first_year and month_i == 0:

                # delete any previous temporary state files created with
                # ms.run()
                try:
                    os.remove(outfile_directory + 'state_tmp.nc')
                except:
                    pass
                params["forcing"] = forcing_file
                params["start"] = startdate
                params["stop"] = enddate
                print('Running forcings from start and stop dates:', startdate,
                      enddate)
                ms = MetSim(params)
                ms.run()

                # delete any previous temporary state files created in wrapper
                # script and save the current
                try:
                    os.remove(state_outfile)
                except:
                    pass
                ms.state.to_netcdf(state_outfile)
            else:
                forcing_filename = (infile_prefix +
                                    "%04d%02d" % (year_current,
                                                   month_current) + '.nc')
                infile = (infile_directory + '/' + forcing_filename)

                # delete any previous temporary state files created with
                # ms.run()
                try:
                    os.remove(outfile_directory + 'state_tmp.nc')
                except:
                    pass

                ms = []
                params['state'] = state_outfile
                params['start'] = startdate
                params['stop'] = enddate
                params['forcing'] = infile

                ms = MetSim(params)
                print('Running forcings from start and stop dates:', startdate,
                      enddate)
                ms.run()

                # delete any previous temporary state files created in wrapper
                # script and save the current
                try:
                    os.remove(state_outfile)
                except:
                    pass
                ms.state.to_netcdf(state_outfile)

            if not spinup:

                print('Concatenating monthly dataset to the yearly dataset')

                # concat the month's dataset to the yearly dataset
                forcing_current = xr.open_dataset(forcing_outfile_tmp)
                if month_i == 0:
                    forcing_concat = forcing_current
                else:
                    forcing_concat = xr.concat([forcing_concat,
                                                forcing_current], dim='time')

            # delete the current month's forcing files
            try:
                os.remove(forcing_outfile_tmp)
            except:
                pass

            date_i += 1
            print('Completed forcings from start and stop dates:',
                  startdate, enddate)
        if spinup:
            # save the state file at the end of spinup
            print('Saving spinup month 12 state file')
            ds = xr.open_dataset(statefile_initial)
            ms.state.coords['time'] = ds.coords['time']

            try:
                os.remove(outfile_directory + 'state_spinup_end.nc')
            except:
                pass
            state_spinup_end_file = outfile_directory + 'state_spinup_end.nc'
            ms.state.to_netcdf(state_spinup_end_file)

            spinup = 0
            date_i = 0
            params['state'] = state_spinup_end_file
        else:
            # output concatenated dataset
            outfile = (outfile_directory + '/' + outfile_prefix +
                       "%04d_" % (year) + 'total.nc')
            print('Saving yearly datasets to outfile:', outfile)
            try:
                os.remove(outfile)
            except:
                pass
            forcing_concat.to_netcdf(outfile, engine='scipy')

            year = year + 1

    # delete any previous temporary state files created with ms.run()
    try:
        os.remove(outfile_directory + 'state_tmp.nc')
    except:
        pass
    state.close()
    forcing_intial.close()
    domain.close()
    computation_end = datetime.datetime.now()
    computation_total = computation_end - computation_start
    print('Total computation time:', computation_total)


if __name__ == "__main__":
    main()

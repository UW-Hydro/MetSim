"""
MTCLIM
"""

import numpy as np

from metsim import params

def run(fname):
    """
    TODO
    """
    print("Trying to do mtclim")
    # Read the data in here 
    # data = read_data(fname)
    data = None
    calc_t_air(data)
    calc_precip(data)
    calc_snowpack(data)
    calc_srad_hum_it(data)
    calc_longwave(data)
    # Write out the data in here
    # write_data(data, new_fname)


def calc_t_air(df):
    """
    TODO
    """
    delta_z = df['site_elev'] - df['base_elev']
    lapse_rates = [df['tmin_lr'], df['tmax_lr']]
    df['s_t_min'] = df['tmin'] + delta_z * lapse_rates[0]
    df['s_t_max'] = df['tmax'] + delta_z * lapse_rates[1]
    df['s_t_min'].where(df['s_t_min'] > df['s_t_max'], 
                        other=df['s_t_max'] - 0.5, 
                        inplace=True)
    t_mean = mean(s_t_min + s_t_max)
    df['s_t_day'] = ((df['s_t_max'] - t_mean) * df['t_air_coeff']) + t_mean


def calc_precip(df):
    """
    TODO
    """
    try:
        factor = df['site_isoh'] / df['base_isoh']
    except:
        factor = 1
    df['s_prcp'] = df['prcp'] * factor


def calc_snowpack(df):
    """
    TODO
    """
    # Some initialization
    df['s_swe'] = 0.
    _simple_snowpack(df)
    n_days = len(df['prcp'])
    
    swe_sum = 0.0
    count = 0
    for i in range(n_days):
        if False: #FIXME
            count += 1
            swe_sum += df['s_swe'][i]

    if count:
        sp = swe_sum/count
        _simple_snowpack(df, snowpack=sp)
     

def _simple_snowpack(df, snowpack=0.0):
    """
    TODO
    """
    n_days = len(df['prcp']) 
    for i in range(n_days):
        if df['s_t_min'][i] <= params['SNOW_TCRIT']:
            snowpack += df['s_prcp'][i]
        else:
            snowpack -= (params['SNOW_TRATE'] * 
                         (df['s_t_min'][i] - params['SNOW_TCRIT']))
        snowpack = np.maximum(snowpack, 0.0)
        df['s_swe'][i] = snowpack


def calc_srad_hum_it():
    pass


def calc_longwave():
    air_temp = 1.0 #FIXME
    emissivity_calc = {
            'TVA'        : lambda x : 0.74 + 0.0049 * x,
            'ANDERSON'   : lambda x : 0.68 + 0.036 * np.power(x, 0.5),
            'BRUTSAERT'  : lambda x : 1.24 * np.power(x/air_temp, 0.14285714),
            'SATTERLUND' : lambda x : 1.08 * (1 - np.exp(-1 * np.power(x, (air_temp/2016)))),
            'IDSO'       : lambda x : 0.7 + 5.95e-5 * x * np.exp(1500/air_temp),
            'PRATA'      : lambda x : (1 - (1 + (46.5 * x/air_temp)) *
                np.exp(-1*np.power((1.2 + 3 * (46.5*x/air_temp)), 0.5))),
            }



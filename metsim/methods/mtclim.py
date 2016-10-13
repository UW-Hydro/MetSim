"""
MTCLIM
"""

import numpy as np

from metsim import PARAMS  as params
from metsim import CONSTS  as consts
from metsim import OPTIONS as options

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


def calc_srad_hum_it(df, tol=0.01, win_type='boxcar'):
    """
    TODO
    """
    n_days = len(df['prcp'])
    daylength = np.zeros(366)
    window = np.zeros(n_days + 90)
    ttmax0 = np.zeros(366)
    flat_potrad = np.zeros(366) 
    slope_potrad = np.zeros(366) 
    t_fmax = np.zeros(n_days)
    df['s_tfmax'] = 0.0

    df['tmax'] = np.maximum(df['tmax'], df['tmin'])
    dtr = df['tmax'] - df['tmin']
    sm_dtr = pd.rolling_window(dtr, window=30, freq='D', 
                               win_type=win_type).fillna(method='bfill')
    
    if n_days <= 30:
        warn('Timeseries is shorter than rolling mean window, filling '
             'missing values with unsmoothed data')
        sm_dtr.fillna(dtr, inplace=True)

    sum_prcp = df['s_prcp'].values.sum()
    ann_prcp = (sum_prcp / n_days) * consts['DAYS_PER_YEAR']
    if ann_prcp == 0.0:
        ann_prcp = 1.0

    if n_days <= 90:
        sum_prcp = df['s_prcp'].values.sum()
        eff_ann_prcp = (sum_prcp / n_days) * consts['DAYS_PER_YEAR']
        eff_ann_prcp = np.maximum(eff_ann_prcp, 8.0)
        parray = eff_ann_prcp
    else:
        start_yday = df.index.dayofyear[0]
        end_yday = df.index.dayofyear[n_days - 1]
        if start_yday != 1:
            if end_yday == start_yday-1:
                isloop = True
        else:
            if end_yday == 365 or end_yday == 366:
                isloop = True
        
        for i in range(90):
            if isloop:
                window[i] = df['s_prcp'][n_days - 90 + i]
            else:
                window[i] = df['s_prcp'][i]
        window[90:] = df['s_prcp']

        for i in range(n_days):
            sum_prcp = 0.0
            for j in range(90):
                sum_prcp += window[i+j]
                sum_prcp = (sum_prcp / 90.) * consts['DAYS_PER_YEAR']


    #------------------------------------------------------------------------------------
    # From here on to end note no data is needed from input - can we put this elsewhere?
    #------------------------------------------------------------------------------------


    trans1 = _calc_trans()
    lat = np.clip(params['site_lat']*consts['RADPERDEG'], -np.pi/2., np.pi/2.0)
    sinlat = np.sin(lat)
    cosslp = np.cos(params['site_slope']  * consts['RADPERDEG'])
    sinslp = np.sin(params['site_slope']  * consts['RADPERDEG'])
    cosasp = np.cos(params['site_aspect'] * consts['RADPERDEG'])
    sinasp = np.sin(params['site_aspect'] * consts['RADPERDEG'])
    coszeh = np.cos(np.pi / 2.-(params['site_east_horiz'] * consts['RADPERDEG']))
    coszwh = np.cos(np.pi / 2.-(params['site_west_horiz'] * consts['RADPERDEG']))
    dt = consts['SRADDT']  
    dh = dt / consts['SECPERRAD']  
    tiny_step_per_day = 86400 / consts['SRADDT']
    tiny_rad_fract = np.zeros(shape=(366, tiny_step_per_day), dtype=np.float64)

   
    for i in range(365):
        decl = consts['MINDECL'] * np.cos((i + consts['DAYSOFF']) *
                                             consts['RADPERDAY'])
        cosdecl = np.cos(decl)
        sindecl = np.sin(decl)

        bsg1 = -sinslp * sinasp * cosdecl
        bsg2 = (-cosasp * sinslp * sinlat + cosslp * coslat) * cosdecl
        bsg3 = (cosasp * sinslp * coslat + cosslp * sinlat) * sindecl
        cosegeom = coslat * cosdecl
        sinegeom = sinlat * sindecl
        coshss = np.clip(-sinegeom / cosegeom, -1, 1)
        hss = np.cos(coshss)  
        daylength[i] = np.maximum(2.0 * hss * consts['SECPERRAD'], 86400)
        sc = 1368.0 + 45.5 * np.sin((2.0 * np.pi * i / 365.25) + 1.7)
        dir_beam_topa = sc * dt
        sum_trans = 0.
        sum_flat_potrad = 0.
        sum_slope_potrad = 0.

        for h in np.arange(-hss, hss, dh):
            cosh = np.cos(h)
            sinh = np.sin(h)
            cza = cosegeom * cosh + sinegeom
            cbsa = sinh * bsg1 + cosh * bsg2 + bsg3

            if (cza > 0.):
                dir_flat_topa = dir_beam_topa * cza

                am = 1.0 / (cza + 0.0000001)
                if (am > 2.9):
                    ami = int((np.cos(cza) / consts['RADPERDEG'])) - 69
                    if (ami < 0):
                        ami = 0
                    if (ami > 20):
                        ami = 20
                    am = consts['OPTAM'][ami]

                trans2 = np.power(trans1, am)
                sum_trans += trans2 * dir_flat_topa
                sum_flat_potrad += dir_flat_topa

                if ((h < 0. and cza > coszeh and cbsa > 0.) or
                        (h >= 0. and cza > coszwh and cbsa > 0.)):
                    sum_slope_potrad += dir_beam_topa * cbsa
            else:
                dir_flat_topa = -1

            tinystep = np.clip(((12 * 3600 + h * consts['SECPERRAD']) /
                                consts['SRADDT']),
                               0, tiny_step_per_day - 1)

            if dir_flat_topa > 0:
                tiny_rad_fract[i, tinystep] = dir_flat_topa
            else:
                tiny_rad_fract[i, tinystep] = 0

        if daylength[i] and sum_flat_potrad > 0:
            tiny_rad_fract[i] /= sum_flat_potrad

        if daylength[i]:
            ttmax0[i] = sum_trans / sum_flat_potrad
            flat_potrad[i] = sum_flat_potrad / daylength[i]
            slope_potrad[i] = sum_slope_potrad / daylength[i]
        else:
            ttmax0[i] = 0.
            flat_potrad[i] = 0.
            slope_potrad[i] = 0.

    ttmax0[365] = ttmax0[364]
    flat_potrad[365] = flat_potrad[364]
    slope_potrad[365] = slope_potrad[364]
    daylength[365] = daylength[364]

    tiny_rad_fract[365] = tiny_rad_fract[364]

    avg_horizon = (params['site_east_horiz'] + params['site_west_horiz']) / 2.0
    horizon_scalar = 1.0 - np.sin(avg_horizon * consts['RADPERDEG'])
    
    if (params['site_slope'] > avg_horizon):
        slope_excess = params['site_slope'] - avg_horizon
    else:
        slope_excess = 0.

    if (2.0 * avg_horizon > 180.):
        slope_scalar = 0.
    else:
        slope_scalar = np.clip(1.-(slope_excess / (180.0-2.0 * avg_horizon)), 0, None)

    sky_prop = horizon_scalar * slope_scalar
    b = params['B0'] + params['B1'] * np.exp(-params['B2'] * sm_dtr)
    t_fmax = 1.0 - 0.9 * np.exp(-b * np.power(dtr, params['C']))


    #------------------------------------------------------------------------------------
    # End portion of the code that doesn't have any reference to the input dataframe
    #------------------------------------------------------------------------------------


    inds = np.nonzero(df['prcp'] > options['SW_PREC_THRESH'])[0]
    t_fmax[inds] *= params['RAIN_SCALAR']
    df['s_tfmax'] = t_fmax

    if 'tdew' in df:
        tdew = df['tdew']
    else:
        tdew = df['s_tmin']
    
    if 's_hum' in df:
        pva = df['s_hum']
    else:
        pva = svp(tdew)

    pa = atm_pres(params['site_elev'])
    yday = df.index.dayofyear - 1
    df['s_dayl'] = daylength[yday]
    tdew_save = tdew
    pva_save = pva

    # FIXME: This function has lots of inputs and outputs 
    tdew, pva, pet = _compute_srad_humidity_onetime(
        tdew, pva, ttmax0, flat_potrad, slope_potrad, sky_prop, daylength,
        parray, pa, dtr)

    sum_pet = pet.values.sum()
    ann_pet = (sum_pet / self.ndays) * consts['DAYS_PER_YEAR'] 

    if (('tdew' in df) or ('s_hum' in df) or
            (options['VP_ITER'].upper() == 'VP_ITER_ANNUAL' and
             ann_pet / ann_prcp >= 2.5)):
        tdew = tdew_save[:]
        pva = pva_save[:]

    if (options['VP_ITER'].upper() == 'VP_ITER_ALWAYS' or
        (options['VP_ITER'].upper() == 'VP_ITER_ANNUAL' and
         ann_pet / ann_prcp >= 2.5) or
            options['VP_ITER'].upper() == 'VP_ITER_CONVERGE'):
        if (options['VP_ITER'].upper() == 'VP_ITER_CONVERGE'):
            max_iter = 100
        else:
            max_iter = 2
    else:
        max_iter = 1

    rmse_tdew = tol + 1

    #FIXME: Strange internal function here just returns another function call
    def f(tdew, *args):
        rmse_tdew = rmse(self._compute_srad_humidity_onetime(tdew, *args),
                         tdew)
        return rmse_tdew

    res = minimize_scalar(f, tdew, args=(pva, ttmax0, flat_potrad,
                                         slope_potrad, sky_prop, daylength,
                                         parray, pa, dtr),
                          tol=rmse_tdew, options={'maxiter': max_iter})
    tdew = res.x
    pva = svp(tdew)
    
    if 's_hum' not in df:
        df['s_hum'] = pva
    
    pvs = svp(df['s_tday'])
    vpd = pvs - pva
    df['s_vpd'] = np.maximum(vpd, 0.)


def _calc_trans():
    """
    TODO
    """
    pratio = np.power((1.0-(consts['LR_STD'] * params['site_elev'])/consts['T_STD']),
                      (consts['G_STD'] / (consts['LR_STD'] * (consts['R'] / consts['MA']))))
    return np.power(params['TBASE'], pratio)


#FIXME: This function has lots of inputs and outputs (see above for call)
def _compute_srad_humidity_onetime(self, tdew, pva, ttmax0, flat_potrad,
                                   slope_potrad, sky_prop, daylength,
                                   parray, pa, dtr):
    """
    TODO
    """
    return (0,0,0)


def calc_longwave(df):
    emissivity_calc = {
            'TVA'        : lambda x : 0.74 + 0.0049 * x,
            'ANDERSON'   : lambda x : 0.68 + 0.036 * np.power(x, 0.5),
            'BRUTSAERT'  : lambda x : 1.24 * np.power(x/air_temp, 0.14285714),
            'SATTERLUND' : lambda x : 1.08 * (1 - np.exp(-1 * np.power(x, (air_temp/2016)))),
            'IDSO'       : lambda x : 0.7 + 5.95e-5 * x * np.exp(1500/air_temp),
            'PRATA'      : lambda x : (1 - (1 + (46.5 * x/air_temp)) *
                np.exp(-1*np.power((1.2 + 3 * (46.5*x/air_temp)), 0.5))),
            }

    air_temp = df['s_t_day'] + consts['KELVIN'] 
    emissivity_clear = emissivity_calc[options['LW_TYPE'].upper()](air_temp)
    tskc = df['s_tsck']

    if options['LW_CLOUD'].upper() == "LW_CLOUD_DEARDORFF":
        emissivity = tsck * 1.0 + (1-tsck) * emissivity_clear
    else:
        emissivity = (1.0 + (0.17 * tsck**2)) * emissivity_clear

    df['s_lwrad'] = emissivity * consts['STEFAN_B'] * np.power(air_temp, 4)



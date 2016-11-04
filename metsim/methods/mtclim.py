"""
MTCLIM
"""

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from statsmodels.tools.eval_measures import rmse

import metsim.disaggregate as disaggregate
from metsim.defaults import PARAMS  as params
from metsim.defaults import CONSTS  as consts
from metsim.defaults import OPTIONS as options

from metsim.physics import svp, calc_pet, atm_pres

n_days = 0

def run(f):
    """
    TODO
    """
    # f is input forcings
    print("Trying to do mtclim")
    # FIXME: Put this somewhere else
    global n_days 
    n_days = len(f['day_of_year'])
    f = f.join(calc_t_air(f['t_min'], f['t_max']))
    f = f.join(calc_precip(f['precip'], f.get('site_isoh', 1), f.get('base_isoh', 1)))
    f = f.join(calc_snowpack(f['day_of_year'], f['s_precip'], f['s_t_min']))
    calc_srad_hum_it(f)
    f = f.join(calc_longwave(f['s_t_day'], f['s_tskc']))


def calc_t_air(t_min, t_max):
    df = pd.DataFrame()
    dZ = params['site_elev'] - params['base_elev']
    lapse_rates = [params['t_min_lr'], params['t_max_lr']]
    df['s_t_max'] = t_max + dZ * lapse_rates[1]
    df['s_t_min'] = np.maximum(t_min + dZ * lapse_rates[0], df['s_t_max']-0.5)
    t_mean = np.mean(df['s_t_min'] + df['s_t_max'])
    df['s_t_day'] = ((df['s_t_max'] - t_mean) * params['TDAYCOEF']) + t_mean
    return df


def calc_precip(precip, site_isoh, base_isoh):
    df = pd.DataFrame()
    df['s_precip'] = precip * (site_isoh / base_isoh) 
    return df


def calc_snowpack(yday, s_precip, s_t_min):

    def _simple_snowpack(s_precip, s_t_min, snowpack=0.0):
        s_swe = np.array(np.ones(n_days) * snowpack)
        accum = s_t_min <= params['SNOW_TCRIT']
        melt  = s_t_min >  params['SNOW_TCRIT']
        s_swe[accum] += s_precip.where(accum)
        s_swe[melt]  -= params['SNOW_TRATE'] * (s_t_min.where(melt) - params['SNOW_TCRIT'])
        s_swe = np.maximum(np.cumsum(s_swe), 0.0) 
        return s_swe 
   
    df = pd.DataFrame()
    df['s_swe'] = _simple_snowpack(s_precip, s_t_min)
    start = (yday == yday[0])
    end = (yday == (start-2)%365 + 1) 
    loop_days = np.logical_or(start, end)
    loop_swe  = sum(df['s_swe'].where(loop_days))/sum(loop_days)

    if np.any(loop_swe):
        df['s_swe'] = _simple_snowpack(s_precip, s_t_min, snowpack=loop_swe)
    return df 


def calc_solar_geom():
    tt_max0 = np.zeros(366)
    daylength = np.zeros(366)
    flat_potrad = np.zeros(366) 
    slope_potrad = np.zeros(366) 
    trans = np.power(params['TBASE'], np.power((1.0-(consts['LR_STD'] * params['site_elev'])/consts['T_STD']),
                 (consts['G_STD'] / (consts['LR_STD'] * (consts['R'] / consts['MA'])))))
    lat    = np.clip(params['site_lat']*consts['RADPERDEG'], -np.pi/2., np.pi/2.0)
    coslat = np.cos(lat)
    sinlat = np.sin(lat)
    cosslp = np.cos(params['site_slope']  * consts['RADPERDEG'])
    sinslp = np.sin(params['site_slope']  * consts['RADPERDEG'])
    cosasp = np.cos(params['site_aspect'] * consts['RADPERDEG'])
    sinasp = np.sin(params['site_aspect'] * consts['RADPERDEG'])
    coszeh = np.cos(np.pi / 2.-(params['site_east_horiz'] * consts['RADPERDEG']))
    coszwh = np.cos(np.pi / 2.-(params['site_west_horiz'] * consts['RADPERDEG']))
    dt     = consts['SRADDT']  
    dh     = dt / consts['SECPERRAD']  
    tiny_step_per_day = 86400 / consts['SRADDT']
    tiny_rad_fract    = np.zeros(shape=(366, tiny_step_per_day), dtype=np.float64)
    for i in range(365):
        decl = consts['MINDECL'] * np.cos((i + consts['DAYSOFF']) * consts['RADPERDAY'])
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
        dir_beam_topa = 1368.0 + 45.5 * np.sin((2.0 * np.pi * i / 365.25) + 1.7) * dt
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
                    ami = min(max(int((np.cos(cza) / consts['RADPERDEG'])) - 69,0),20)
                    am = consts['OPTAM'][ami]
                sum_trans += np.power(trans, am) * dir_flat_topa
                sum_flat_potrad += dir_flat_topa
                # FIXME: This is a long conditional
                if ((h < 0. and cza > coszeh and cbsa > 0.) or
                        (h >= 0. and cza > coszwh and cbsa > 0.)):
                    sum_slope_potrad += dir_beam_topa * cbsa
            else:
                dir_flat_topa = -1
            tinystep = np.clip(((12 * 3600 + h * consts['SECPERRAD'])/dt), 0, tiny_step_per_day - 1)
            tiny_rad_fract[i, tinystep] = max(dir_flat_topa, 0)
        if daylength[i] and sum_flat_potrad > 0:
            tiny_rad_fract[i] /= sum_flat_potrad
        if daylength[i]:
            tt_max0[i] = sum_trans / sum_flat_potrad
            flat_potrad[i] = sum_flat_potrad / daylength[i]
            slope_potrad[i] = sum_slope_potrad / daylength[i]
        else:
            tt_max0[i] = 0.
            flat_potrad[i] = 0.
            slope_potrad[i] = 0.
    tt_max0[365] = tt_max0[364]
    flat_potrad[365] = flat_potrad[364]
    slope_potrad[365] = slope_potrad[364]
    daylength[365] = daylength[364]
    tiny_rad_fract[365] = tiny_rad_fract[364]
    return tt_max0, flat_potrad, slope_potrad, daylength, tiny_rad_fract


def calc_shortwave(df):
    pass


def calc_srad_hum_it(df, tol=0.01, win_type='boxcar'):
    """
    TODO
    """
    window = np.zeros(n_days + 90)
    t_fmax = np.zeros(n_days)
    df['s_tfmax'] = 0.0

    df['t_max'] = np.maximum(df['t_max'], df['t_min'])
    dtr = df['t_max'] - df['t_min']
    sm_dtr = pd.rolling_window(dtr, window=30, freq='D', 
                               win_type=win_type).fillna(method='bfill')
    if n_days <= 30:
        print('Timeseries is shorter than rolling mean window, filling ')
        print('missing values with unsmoothed data')
        sm_dtr.fillna(dtr, inplace=True)

    sum_precip = df['s_precip'].values.sum()
    ann_precip = (sum_precip / n_days) * consts['DAYS_PER_YEAR']
    if ann_precip == 0.0:
        ann_precip = 1.0

    if n_days <= 90:
        sum_precip = df['s_precip'].values.sum()
        eff_ann_precip = (sum_precip / n_days) * consts['DAYS_PER_YEAR']
        eff_ann_precip = np.maximum(eff_ann_precip, 8.0)
        parray = eff_ann_precip
    else:
        parray = np.zeros(n_days)
        start_yday = df['day_of_year'][0]
        end_yday = df['day_of_year'][-1]
        if start_yday != 1:
            if end_yday == start_yday-1:
                isloop = True
        else:
            if end_yday == 365 or end_yday == 366:
                isloop = True
        
        if isloop:
            for i in range(90):
                window[i] = df['s_precip'][n_days - 90 + i]
        else:
            for i in range(90):
                window[i] = df['s_precip'][i]
        window[90:] = df['s_precip']

        for i in range(n_days):
            sum_precip = 0.0
            for j in range(90):
                sum_precip += window[i+j]
                sum_precip = (sum_precip / 90.) * consts['DAYS_PER_YEAR']
            sum_precip = np.maximum(sum_precip, 8.0)
            parray[i] = sum_precip

    # FIXME: This is still bad form
    tt_max0, flat_potrad, slope_potrad, daylength, tiny_rad_fract = calc_solar_geom()
    # NOTE: Be careful with this one!
    disaggregate.tiny_rad_fract = tiny_rad_fract

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
    inds = np.nonzero(df['precip'] > options['SW_PREC_THRESH'])[0]
    t_fmax[inds] *= params['RAIN_SCALAR']
    df['s_tfmax'] = t_fmax

    tdew = df.get('tdew', df['s_t_min'])
    pva = df['s_hum'] if 's_hum' in df else svp(tdew)

    pa = atm_pres(params['site_elev'])
    yday = df['day_of_year'] - 1 
    df['s_dayl'] = daylength[yday]
    tdew_save = tdew
    pva_save = pva

    # FIXME: This function has lots of inputs and outputs 
    tdew, pva, pet = _compute_srad_humidity_onetime(
               tdew, pva, tt_max0, flat_potrad, 
               slope_potrad, sky_prop, daylength,
        parray, pa, dtr, df)

    sum_pet = pet.values.sum()
    ann_pet = (sum_pet / n_days) * consts['DAYS_PER_YEAR'] 

    # FIXME: Another really long conditional
    if (('tdew' in df) or ('s_hum' in df) or
            (options['VP_ITER'].upper() == 'VP_ITER_ANNUAL' and
             ann_pet / ann_precip >= 2.5)):
        tdew = tdew_save[:]
        pva = pva_save[:]

    # FIXME: Another really long conditional
    #if (options['VP_ITER'].upper() == 'VP_ITER_ALWAYS' or
    #    (options['VP_ITER'].upper() == 'VP_ITER_ANNUAL' and
    #     ann_pet / ann_precip >= 2.5) or
    #        options['VP_ITER'].upper() == 'VP_ITER_CONVERGE'):
    #    if (options['VP_ITER'].upper() == 'VP_ITER_CONVERGE'):
    #        max_iter = 100
    #    else:
    #        max_iter = 2
    #else:
    #    max_iter = 1
    
    #FIXME Still want to reduce the number of args here
    #FIXME This also takes up the majority of the mtclim runtime
    rmse_tdew = tol + 1
    #f = lambda x : rmse(_compute_srad_humidity_onetime(x, pva, tt_max0, flat_potrad,
    #                                     slope_potrad, sky_prop, daylength,
    #                                     parray, pa, dtr, df)[0], tdew)
    def f(x):
        tdew_calc = _compute_srad_humidity_onetime(x, pva, tt_max0, flat_potrad,
                                         slope_potrad, sky_prop, daylength,
                                         parray, pa, dtr, df)[0]
        print(tdew_calc - tdew)
        err = rmse(tdew_calc, tdew)
        print(err)
        return err

    res = minimize(f, tdew, tol=rmse_tdew)
    tdew = res.x
    pva = svp(tdew)
    if 's_hum' not in df:
        df['s_hum'] = pva
    
    pvs = svp(df['s_t_day'])
    vpd = pvs - pva
    df['s_vpd'] = np.maximum(vpd, 0.)


#FIXME: This function has lots of inputs and outputs (see above for call)
#FIXME: This needs to not use module level variables inside
def _compute_srad_humidity_onetime(tdew, pva, tt_max0, flat_potrad,
                                   slope_potrad, sky_prop, daylength,
                                   parray, pa, dtr, df):
    """
    TODO
    """
    yday = df['day_of_year'] - 1
    t_tmax = np.minimum(tt_max0[yday] + (params['ABASE'] * pva), 0.0001)
    df['s_ttmax'] = t_tmax
    t_final = t_tmax * df['s_tfmax']
    df['s_fdir'] = 1.0 - np.clip(-1.25 * t_final * 1.25, 0., 1.) 

    srad1 = slope_potrad[yday] * t_final * df['s_fdir']
    srad2 = (flat_potrad[yday] * t_final * 1 - df['s_fdir']) * \
        (sky_prop + params['DIF_ALB'] * (1.0 - sky_prop))

    sc = np.zeros_like(df['s_swe'])
    if (options['MTCLIM_SWE_CORR']):
        inds = np.nonzero(df['s_swe'] > 0. & daylength[yday] > 0.)
        sc[inds] = (1.32 + 0.096 * df['s_swe'][inds]) *\
            1.0e6 / daylength[yday][inds]
        sc = np.maximum(sc, 100.)  # JJH - this is fishy 
                                   # ARB - I agree

    if 's_swrad' in df:
        potrad = (srad1 + srad2 + sc) * daylength[yday] / t_final / 86400
        df['s_tfmax'] = np.ones(len(sc)) 
        inds = np.nonzero((potrad > 0.) & 
               (df['s_swrad'] > 0.) & 
               (daylength[yday] > 0))[0]
        df['s_tfmax'][inds] = np.maximum((df['s_swrad'][inds] / (potrad[inds] * t_tmax[inds])), 1.)
    else:
        df['s_swrad'] = srad1 + srad2 + sc

    if (options['LW_CLOUD'].upper() == 'CLOUD_DEARDORFF'):
        df['s_tskc'] = (1. - df['s_tfmax'])
    else:
        df['s_tskc'] = np.sqrt((1. - df['s_tfmax']) / 0.65)

    # Compute PET using SW radiation estimate, and update Tdew, pva **
    tmink = df['s_t_min'] + consts['KELVIN']
    pet = calc_pet(df['s_swrad'], df['s_t_day'], pa,
                   df['s_dayl'])

    # calculate ratio (PET/effann_prcp) and correct the dewpoint
    ratio = pet / parray
    df['s_ppratio'] = ratio * 365.25
    tdewk = tmink * (-0.127 + 1.121 *
                     (1.003 - 1.444 * ratio + 12.312 *
                      np.power(ratio, 2) - 32.766 * np.power(ratio, 3)) +
                     0.0006 * dtr)
    tdew_tmp = tdewk - consts['KELVIN']
    return tdew_tmp, pva, pet


def calc_longwave(s_t_day, s_tskc):
    emissivity_calc = {
            'DEFAULT'    : lambda x : x,
            'TVA'        : lambda x : 0.74 + 0.0049 * x,
            'ANDERSON'   : lambda x : 0.68 + 0.036 * np.power(x, 0.5),
            'BRUTSAERT'  : lambda x : 1.24 * np.power(x/air_temp, 0.14285714),
            'SATTERLUND' : lambda x : 1.08 * (1 - np.exp(-1 * np.power(x, (air_temp/2016)))),
            'IDSO'       : lambda x : 0.7 + 5.95e-5 * x * np.exp(1500/air_temp),
            'PRATA'      : lambda x : (1 - (1 + (46.5 * x/air_temp)) *
                np.exp(-1*np.power((1.2 + 3 * (46.5*x/air_temp)), 0.5))),
            }
    cloud_calc = {
            'DEFAULT' : lambda x : (1.0 + (0.17 * s_tskc**2)) * x,
            'CLOUD_DEARDORFF' : lambda x : s_tskc * 1.0 + (1-s_tskc) * x
            }
    df = pd.DataFrame()
    air_temp = s_t_day + consts['KELVIN'] 
    emissivity_clear = emissivity_calc[options['LW_TYPE'].upper()](air_temp)
    emissivity = cloud_calc[options['LW_CLOUD'].upper()](emissivity_clear) 
    df['s_lwrad'] = emissivity * consts['STEFAN_B'] * np.power(air_temp, 4)
    return df



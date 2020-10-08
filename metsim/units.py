import metsim.constants as cnst

converters = {
        'pet': {
            'mm timestep-1': lambda x, ts: x,
            'mm h-1': lambda x, ts: x / (ts / cnst.MIN_PER_HOUR),
            'mm s-1': lambda x, ts: x / (ts * cnst.SEC_PER_MIN),
            },
        'prec': {
            'mm timestep-1': lambda x, ts: x,
            'mm h-1': lambda x, ts: x / (ts / cnst.MIN_PER_HOUR),
            'mm s-1': lambda x, ts: x / (ts * cnst.SEC_PER_MIN),
            },
        'shortwave': {
            'W m-2': lambda x, ts: x,
            },
        'longwave': {
            'W m-2': lambda x, ts: x,
            },
        't_max': {
            'C': lambda x, ts: x,
            'K': lambda x, ts: x + 273.15,
            },
        't_min': {
            'C': lambda x, ts: x,
            'K': lambda x, ts: x + 273.15,
            },
        'temp': {
            'C': lambda x, ts: x,
            'K': lambda x, ts: x + 273.15,
            },
        'vapor_pressure': {
            'Pa': lambda x, ts: x,
            'hPa': lambda x, ts: x / 100.,
            'kPa': lambda x, ts: x / 1000.,
            },
        'air_pressure': {
            'kPa': lambda x, ts: x,
            'hPa': lambda x, ts: x * 100.,
            'Pa': lambda x, ts: x * 1000.,
            },
        'tskc': {
            'fraction': lambda x, ts: x,
            '%': lambda x, ts: x * 100.,
            },
        'rel_humid': {
            '%': lambda x, ts: x,
            'fraction': lambda x, ts: x / 100.,
            },
        'spec_humid': {
            'g g-1': lambda x, ts: x,
            },
        'wind': {
            'm s-1': lambda x, ts: x,
            },
        'daylength': {
            's': lambda x, ts: x,
            }
        }

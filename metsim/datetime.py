import re
import numpy as np
import pandas as pd
import xarray as xr
from netCDF4 import num2date, date2num
DEFAULT_ORIGIN = '0001-01-01'


def date_range(start=None, end=None, periods=None, freq='D', tz=None,
               normalize=False, name=None, closed=None, calendar='standard',
               **kwargs):
    ''' Return a fixed frequency datetime index, with day (calendar) as the
    default frequency

    Parameters
    ----------
    start : string or datetime-like, default None
        Left bound for generating dates
    end : string or datetime-like, default None
        Right bound for generating dates
    periods : integer or None, default None
        If None, must specify start and end
    freq : string or DateOffset, default 'D' (calendar daily)
        Frequency strings can have multiples, e.g. '5H'
    tz : string or None
        Time zone name for returning localized DatetimeIndex, for example
        Asia/Hong_Kong
    normalize : bool, default False
        Normalize start/end dates to midnight before generating date range
    name : str, default None
        Name of the resulting index
    closed : string or None, default None
        Make the interval closed with respect to the given frequency to
        the 'left', 'right', or both sides (None)
    calendar : string
        Describes the calendar used in the time calculations. Default is a the
        standard calendar (with leap years)

    Notes
    -----
    2 of start, end, or periods must be specified
    To learn more about the frequency strings, please see `this link
    <http://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases>`__.

    Returns
    -------
    rng : DatetimeIndex
    '''
    if calendar in ['standard', 'gregorian', 'propoleptic_gregorian']:
        return pd.date_range(start=start, end=end, periods=periods,
                             freq=freq, tz=tz, normalize=normalize, name=name,
                             closed=closed, **kwargs)
    else:
        # start and end are give
        if (start is not None) and (end is not None) and (periods is None):

            steps, units = decode_freq(freq)
            start_num, end_num = date2num(
                pd.to_datetime([start, end]).to_pydatetime(),
                units, calendar=calendar)
            periods = int((end_num - start_num) / steps)  # Todo divide by freq

            times = num2date(
                np.linspace(start_num, end_num, periods,
                            endpoint=False,
                            dtype=np.float128), units, calendar)

            index = pd.DatetimeIndex(xr.conventions.nctime_to_nptime(times))

            return index

        else:
            raise NotImplementedError(
                'Specified arguments are not valid for this calendar')


def decode_freq(freq):
    if len(freq) > 1:
        r = re.compile('([0-9]+)([a-zA-Z]+)')
        step, unit = r.match(freq).groups()
    else:
        step = 1
        unit = freq
    return (int(step), units_from_freq(unit))


def units_from_freq(freq, origin=DEFAULT_ORIGIN):
    if 'H' in freq:
        return 'hours since %s' % origin
    elif 'D' in freq:
        return 'days since %s' % origin
    else:
        raise NotImplementedError()

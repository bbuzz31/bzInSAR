""" Set of functions to deal with mintpy objects """
from VLM.bzFRInGE import *
import h5py
import pandas as pd


# mostly copied from bbMP get_ifg_info
# close, but not same as the coherenceSpatialAvg.txt from mintpy
def stack2df(path_stack, drop=True, ret_dat=False):
    """ Convert ifgramStack to a dataframe """
    with h5py.File(path_stack, 'r') as h5:
        dt    = h5['date'][:]
        refs  = [j.decode('utf-8') for j in dt[:, 0]]
        secs  = [j.decode('utf-8') for j in dt[:, 1]]
        dropI = h5['dropIfgram'][:]
        bperp = h5['bperp'][:]
        try:
            coh   = h5['coherence'][:]
            coh   = np.nanmean(np.where(np.isclose(coh, 0), np.nan, coh), axis=(1, 2))
        except:
            pass

        if ret_dat:
            data = h5[ret_dat][:]
        else:
            data = None

    if drop:
        refs  = np.array(refs)[dropI].tolist()
        secs  = np.array(secs)[dropI].tolist()
        bperp = np.array(bperp)[dropI].tolist()
        try:
            coh   = coh[dropI]
        except:
            pass

        if ret_dat:
            data = data[dropI]


    df = pd.DataFrame({'ref': refs, 'sec': secs, 'bperp': bperp})

    try:
        df['coh'] = coh
    except:
        pass

    df = df.sort_values('ref sec'.split())
    return df, data


def get_ifg_dates(path_stack):
    """ Get ifg dates from an ifgramStack """
    with h5py.File(path_stack, 'r') as h5:
        dates = h5['date'][:]
    ifgs = [f"{dt[0].decode('utf-8')}_{dt[1].decode('utf-8')}" for dt in dates]
    return ifgs

def get_ts_dates(Exp):
    """ Get the dates of the timeseries for GACOS  """
    with h5py.File(Exp.path_ts, 'r') as h5:
        dates = h5['date'][:]
    dates = [dt.decode() for dt in dates] # gets rid of the leading b
    return dates

if __name__ == '__main__':
    Exp = ExpBase(NYC_Base)
    dates = get_ts_dates(Exp)
    with open ('./tmp.txt', 'w') as fh:
        [print (dt, file=fh) for dt in dates]
    print (Exp.SNWE)

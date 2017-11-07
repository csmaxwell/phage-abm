
from pandas import HDFStore



def read_data(model, table):
    with HDFStore(model) as hdf:
        tst = hdf.select(table)
    try:
        #Step is a variable used by the mesa package
        tst.Step = tst.Step.astype(int)
    except AttributeError:
        pass
    tst.re_degrade_foreign = tst.re_degrade_foreign.astype(float)
    return tst

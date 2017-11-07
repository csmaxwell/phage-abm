import os
from pandas import HDFStore
from uuid import uuid4

def hdf_exist_p(store, data_name):
    """Tests if hdf5 storage exists and contains the data_name"""
    if not data_name.startswith("/"):
        data_name = "/" + data_name
    if os.path.isfile(store):
        with HDFStore(store) as hdf:
            return data_name in hdf.keys()
    else:
        return False
    
    
def write_to_hdf(store_path, df, data_name, hash_id=None):
    """Adds a hash_ID to the results, writes to hdf5. Initializes file if
    needed.

    """
    if not hash_id:
        df["hash"] = uuid4().hex
    else:
        df["hash"] = hash_id
    cols = df.columns
    if not hdf_exist_p(store_path, data_name):
        with HDFStore(store_path) as hdf:
            hdf.put(data_name, df, columns = cols, 
                    format='table', data_columns=True)
    else:
        with HDFStore(store_path) as hdf:
            hdf.append(data_name, df, format='table', data_columns=True)

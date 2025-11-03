import os
# import glob
import h5py
from datetime import datetime
from numpy import ndarray
from pint import UnitRegistry

u = UnitRegistry()
Q_ = u.Quantity


timestamp_format = '%Y-%m-%d-%H-%M-%S'  # used by `datetime.strftime` in `new_timestamp()`

def new_timestamp():
    return datetime.strftime(datetime.now(),timestamp_format)

"""
HDF5 utilities for unitful quantities and arrays
"""
def dump_hdf5(ds,fpath,open_mode='a'):
    with h5py.File(fpath, open_mode) as f:
        for k,v in ds.items():
            print("dumping item: " + k)
            try:
                units = v.units
                v = v.m
            except:
                units = False
            if type(v) is ndarray:
                h5_ds = f.create_dataset(k,
                                        v.shape,
                                        dtype=v.dtype,
                                        data=v,
                                        compression='gzip',
                                        )
                if units:
                    h5_ds.attrs['units'] = str(units)
            else:  #
                if type(v) is list:
                    for ind,item in enumerate(v):
                        try:
                            item_units = item.units
                            item_val = item.m
                        except:
                            item_units = False
                        if item_units:
                            v[ind] = (item_val,str(item_units))
                if units:
                    print("non-array, non-list item with units:")
                    print(units)
                    # f.attrs[k] = (v,str(units))
                    f.attrs[k] = v
                    f.attrs["_".join([k,"units"])] = str(units)
                else:
                    print("non-array, non-list item without units")
                    f.attrs[k] = v
        f.flush()

def load_hdf5(fpath=None,dir=None,filter=None,sim_index=None):
    if fpath is None:
        # sim_id = ''.join(str(ind)+'-' for ind in self.param_index_combinations[sim_index])[:-1]
        file_list =  glob(os.path.normpath(dir)+os.path.normpath('/'+ sim_id + '*'))
        fpath = max(file_list,key=os.path.getctime)
    print('loading file: ' + os.path.basename(fpath))
    # with open(fpath, "rb") as f:
    #     ds = pickle.load(f)
    # for k,v in ds.items():
    #     try:
    #         ds[k] = u.Quantity.from_tuple(v)
    #     except:
    #         ds[k] = v
    ds = {}
    with h5py.File(fpath, "r") as f:
        for h5_ds_name in f:
            print('importing ' + h5_ds_name + '...')
            ds[h5_ds_name] = _load_hdf5_item(f[h5_ds_name])
        for key,val in f.attrs.items():
            print('importing attr ' + key + '...')
            # try:
            #     ds[key] = u.Quantity.from_tuple(val)
            # except:
            #     ds[key] = val
            if "_".join([key,"units"]) in f.attrs.keys():
                units_str = str(f.attrs["_".join([key,"units"])])
                print(units_str)
                ds[key] = u.Quantity(val,units_str)
            elif key.split("_")[-1]=="units":
                pass
            else:
                ds[key] = val
    return ds

def _load_hdf5_item(item):
    if issubclass(type(item),h5py.Group):
        ds = {subname: _load_hdf5_item(item[subname]) for subname in item}
        for key,val in item.attrs.items():
            try:
                ds[key] = u.Quantity.from_tuple(val)
            except:
                ds[key] = val
        return ds
    else:
        if 'units' in item.attrs:
            # ds = Q_(item[()],item.attrs['units'])

            try:
                ds = Q_(item[()],item.attrs['units'])
            except:
                ds = item[()]
            return ds
        else:
            return item[()]

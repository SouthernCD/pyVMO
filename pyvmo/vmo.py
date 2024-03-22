from pyvmo.format import vcf_to_vmo
from pyvmo.operation import extract_sub_vmo, get_valid_variant_index
from toolbiox.lib.common.os import mkdir
from toolbiox.lib.common.util import pickle_load
import numpy as np
import pandas as pd


class VMO(object):
    """
    VariantMemmapObject
    """

    def __init__(self, vmo_dir):
        self.vmo_dir = vmo_dir
        mkdir(self.vmo_dir, True)

        self.info_store = "%s/info_store.h5" % self.vmo_dir
        self.variants_info_store = "%s/variants_info_store.h5" % self.vmo_dir
        self.samples_info_store = "%s/samples_info_store.h5" % self.vmo_dir
        self.memmap_file_path = "%s/memmap.dat" % self.vmo_dir
        self.memmap_size_file_path = "%s/memmap_size.pkl" % self.vmo_dir

    def store_vcf_to_vmo(self, vcf_file, chunk_size=10000, force_update=False):
        vcf_to_vmo(vcf_file, self, chunk_size, force_update)

    # sample info
    def get_sample_info(self):
        return pd.read_hdf(self.samples_info_store, key='samples').samples.tolist()

    def get_samples_index(self, target_sample_ids):
        sample_ids = self.get_sample_info()
        return [sample_ids.index(id) for id in target_sample_ids]

    # variant info
    def get_variant_info(self):
        return pd.read_hdf(self.variants_info_store, key='varants').apply(lambda x: x.astype(str), axis=0)

    # genotype matrix
    def get_matrix(self):
        shape, dtype = pickle_load(self.memmap_size_file_path)
        return np.memmap(self.memmap_file_path, dtype=dtype, mode='r', shape=shape)

    # build submatrix
    def extract_sub_vmo(self, sub_vmo_dir, var_idx_list=None, spl_idx_list=None, chunk_size=10000, force_update=False):
        sub_vmo_dir = extract_sub_vmo(self, sub_vmo_dir, var_idx_list, spl_idx_list, chunk_size, force_update)
        return VMO(sub_vmo_dir)

    # site filter
    def site_filter(self, maf_thr=0.05, mis_thr=0.5, chunk_size=1000, n_jobs=8):
        m = self.get_matrix()
        valid_variant_bool_index = get_valid_variant_index(
            m, maf_thr=maf_thr, mis_thr=mis_thr, chunk_size=chunk_size, n_jobs=n_jobs)
        return np.where(valid_variant_bool_index)[0]


if __name__ == '__main__':
    vcf_file = "my_vcf_file"
    vmo_dir = "my_vmo_dir"
    sample_id_list = ["sample1", "sample2", "sample3"]

    # Store raw vcf to vmo, vcf to hdf5 about 5h (64618398, 1757, 2) and matrix to memmap about 45min
    vmo = VMO(vmo_dir)
    vmo.store_vcf_to_vmo(vcf_file)

    # Extract submatrix for used samples
    spl_idx_list = vmo.get_samples_index(sample_id_list)
    sub_vmo_dir = "sub_vmo_dir"
    # about 20min (331 samples)
    sub_vmo = vmo.extract_sub_vmo(sub_vmo_dir, spl_idx_list=spl_idx_list, chunk_size=10000)

    # get_valid_variant_index about 2hour
    valid_variant_idx_list = sub_vmo.site_filter(maf_thr=0.05, mis_thr=0.5)
    len(valid_variant_idx_list)

    # Extract submatrix for valid variants
    passed_vmo_dir = "passed_vmo_dir"
    # about 5min
    passed_vmo = sub_vmo.extract_sub_vmo(passed_vmo_dir, var_idx_list=valid_variant_idx_list, chunk_size=10000)

    # get matrix from vmo
    m = passed_vmo.get_matrix()
    # force load matrix into memory
    mm = np.array(m)

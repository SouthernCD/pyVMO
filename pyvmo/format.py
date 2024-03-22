from cyvcf2 import VCF
from toolbiox.lib.common.os import mkdir, have_file
from toolbiox.lib.common.util import pickle_dump
from pyvmo.operation import get_mis_ref_alt_num_parallel
import allel
import gc
import h5py
import numpy as np
import pandas as pd
import time


def vcf_to_vmo(input_vcf_file, output_vmo, chunk_size=10000, force_update=False):
    mkdir(output_vmo.vmo_dir, False if force_update else True)
    log_file = "%s/store.log" % output_vmo.vmo_dir

    # Store vcf to hdf5
    if not have_file(output_vmo.info_store) or force_update:
        with open(log_file, 'w') as logger_handle:
            allel.vcf_to_hdf5(input_vcf_file, output_vmo.info_store,
                              fields='*', log=logger_handle)

    # Store vcf to memmap
    if not have_file(output_vmo.memmap_file_path) or force_update:
        with open(log_file, 'a') as logger_handle:
            block_parsed_num = 0
            last_printed_ratio = 0
            with h5py.File(output_vmo.info_store, 'r') as f:
                # 获取数据集
                dataset = f['calldata/GT']

                # 如果未指定行和列，则使用所有行和列
                rows = range(len(dataset))
                columns = range(dataset.shape[1])

                # 创建一个内存映射文件
                memmap = np.memmap(output_vmo.memmap_file_path, dtype=dataset.dtype,
                                   mode='w+', shape=(len(rows), len(columns), 2))
                del memmap
                gc.collect()

                # 分块读写数据
                total_chunks = len(
                    range(0, len(rows), chunk_size)) * len(range(0, len(columns), chunk_size))

                for i in range(0, len(rows), chunk_size):
                    for j in range(0, len(columns), chunk_size):
                        # 创建一个内存映射文件
                        memmap = np.memmap(
                            output_vmo.memmap_file_path, dtype=dataset.dtype, mode='r+', shape=(len(rows), len(columns), 2))
                        row_slice = slice(
                            i, min(i + chunk_size, len(rows)))
                        col_slice = slice(
                            j, min(j + chunk_size, len(columns)))
                        memmap[row_slice,
                               col_slice] = dataset[row_slice, col_slice]
                        memmap.flush()
                        del memmap
                        gc.collect()
                        # dataset.read_direct(memmap, np.s_[row_slice, col_slice], np.s_[row_slice, col_slice])
                        block_parsed_num += 1
                        parsed_ratio = block_parsed_num/total_chunks*100
                        if parsed_ratio - last_printed_ratio >= 1:
                            logger_handle.write("%s: parsed %d blocks (%.2f%%)\n" % (time.strftime(
                                "%Y-%m-%d %H:%M:%S", time.localtime()), block_parsed_num, block_parsed_num*chunk_size*chunk_size/len(rows)/len(columns)*100))
                            logger_handle.flush()
                            last_printed_ratio = parsed_ratio

    # Store memmap size
    if not have_file(output_vmo.memmap_size_file_path) or force_update:
        with h5py.File(output_vmo.info_store, 'r') as f:
            dataset = f['calldata/GT']
            shape = dataset.shape
            dtype = dataset.dtype

        pickle_dump((shape, dtype), output_vmo.memmap_size_file_path)

    # Store variant info into df hdf5

    if not have_file(output_vmo.variants_info_store) or force_update:
        with open(log_file, 'a') as logger_handle:
            vcf_reader = VCF(input_vcf_file)
            total_num = vcf_reader.num_records
            num = 0
            variants_list = []
            for variant in vcf_reader():
                # variant
                variants_list.append({
                    'CHROM': variant.CHROM,
                    'POS': variant.POS,
                    'REF': variant.REF,
                    'ALT': ','.join([str(a) for a in variant.ALT]),
                })

                num += 1
                if num % 100000 == 0:
                    logger_handle.write("%s: read %d/%d variants (%.2f%%)\n" % (time.strftime(
                        "%Y-%m-%d %H:%M:%S", time.localtime()), num, total_num, num/total_num*100))
                    logger_handle.flush()

            df = pd.DataFrame(variants_list)
            df.index.name = 'index'
            df.to_hdf(output_vmo.variants_info_store,
                      key='varants', mode='w', format='table')

    # Store sample info into df hdf5
    if not have_file(output_vmo.samples_info_store) or force_update:
        with h5py.File(output_vmo.info_store, 'r') as f:
            sample_ids = f['samples'][:]

        df = pd.DataFrame([i.decode()
                           for i in sample_ids], columns=['samples'])
        df.to_hdf(output_vmo.samples_info_store, key='samples', mode='w')

    return output_vmo


def vmo_to_vcf(vmo):
    pass


def vmo_to_bimbam(vmo, bimbam_file, chunk_size=1000, n_jobs=8):

    m = vmo.get_matrix()

    var_df = vmo.get_variant_info()

    mis_num, ref_num, alt_num = get_mis_ref_alt_num_parallel(m, chunk_size=chunk_size, n_jobs=n_jobs)

    with open(bimbam_file, 'w') as bimbam_handle:

        for i in range(m.shape[0]):
            chr_id = var_df.iloc[i]['CHROM']
            pos = int(var_df.iloc[i]['POS'])
            ref = var_df.iloc[i]['REF']
            alt = var_df.iloc[i]['ALT']

            minor_allele = alt if ref_num[i] > alt_num[i] else ref
            major_allele = ref if ref_num[i] > alt_num[i] else alt
            major_allele_is_ref = ref_num[i] > alt_num[i]

            var_id = "%s_%d_%s/%s" % (chr_id, pos, minor_allele, major_allele)

            bimbam_col_list = [var_id, minor_allele, major_allele] +  list(np.sum(m[i], axis=1) if major_allele_is_ref else  np.sum(np.abs(m[i] - 1),axis=1))

            bimbam_handle.write(",".join([str(i) for i in bimbam_col_list]) + "\n")


if __name__ == "__main__":
    from pyvmo import VMO
    from pyvmo.operation import get_mis_ref_alt_num_parallel
    import numpy as np
    
    vmo = VMO(
        "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/1.raw_data/beagle/reseq_landrace_passed_vmo")

    bimbam_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/1.raw_data/beagle/reseq_landrace_passed.bimbam"
    vmo_to_bimbam(vmo, bimbam_file, chunk_size=1000, n_jobs=20)
from joblib import Parallel, delayed
from toolbiox.lib.common.os import have_file
from toolbiox.lib.common.util import pickle_dump
import gc
import h5py
import numpy as np
import pandas as pd
import time

# build submatrix


def extract_sub_vmo(vmo, sub_vmo_dir, var_idx_list=None, spl_idx_list=None, chunk_size=10000, force_update=False):

    variants_info_store = "%s/variants_info_store.h5" % sub_vmo_dir
    samples_info_store = "%s/samples_info_store.h5" % sub_vmo_dir
    memmap_file_path = "%s/memmap.dat" % sub_vmo_dir
    memmap_size_file_path = "%s/memmap_size.pkl" % sub_vmo_dir

    log_file = "%s/store.log" % sub_vmo_dir

    # build submatrix
    if not have_file(memmap_file_path) or force_update:
        # 获取数据集的大小和类型
        vmo_m = vmo.get_matrix()
        vmo_v, vmo_s, vmo_h = vmo_m.shape
        vmo_dtype = vmo_m.dtype
        del vmo_m
        gc.collect()

        # 如果未指定行和列，则使用所有行和列
        if var_idx_list is None:
            var_idx_list = range(vmo_v)
        if spl_idx_list is None:
            spl_idx_list = range(vmo_s)

        # 创建一个内存映射文件
        memmap = np.memmap(memmap_file_path, dtype=vmo_dtype,
                           mode='w+', shape=(len(var_idx_list), len(spl_idx_list), vmo_h))

        # 分块读写数据
        block_parsed_num = 0
        last_printed_ratio = 0

        total_chunks = len(range(0, len(var_idx_list), chunk_size)) * \
            len(range(0, len(spl_idx_list), chunk_size))

        with open(log_file, 'a') as logger_handle:
            for i in range(0, len(var_idx_list), chunk_size):
                for j in range(0, len(spl_idx_list), chunk_size):
                    row_slice = var_idx_list[i:min(
                        i + chunk_size, len(var_idx_list))]
                    col_slice = spl_idx_list[j:min(
                        j + chunk_size, len(spl_idx_list))]
                    memmap = np.memmap(memmap_file_path, dtype=vmo_dtype,
                                       mode='r+', shape=(len(var_idx_list), len(spl_idx_list), vmo_h))
                    vmo_m = vmo.get_matrix()
                    memmap[i:min(i + chunk_size, len(var_idx_list)), j:min(
                        j + chunk_size, len(spl_idx_list)), :] = vmo_m[row_slice, :][:, col_slice, :]
                    memmap.flush()
                    del memmap
                    del vmo_m
                    gc.collect()
                    block_parsed_num += 1
                    parsed_ratio = block_parsed_num/total_chunks*100
                    if parsed_ratio - last_printed_ratio >= 1:
                        logger_handle.write("%s: parsed %d blocks (%.2f%%)\n" % (time.strftime(
                            "%Y-%m-%d %H:%M:%S", time.localtime()), block_parsed_num, parsed_ratio))
                        logger_handle.flush()
                        last_printed_ratio = parsed_ratio

        pickle_dump(((len(var_idx_list), len(spl_idx_list), vmo_h),
                    vmo_dtype), memmap_size_file_path)

    # build variant info
    if not have_file(variants_info_store) or force_update:
        with h5py.File(vmo.variants_info_store, 'r') as f:
            df = pd.read_hdf(vmo.variants_info_store, key='varants')
            df.iloc[var_idx_list].to_hdf(variants_info_store,
                                         key='varants', mode='w', format='table')

    # build sample info
    if not have_file(samples_info_store) or force_update:
        with h5py.File(vmo.samples_info_store, 'r') as f:
            df = pd.read_hdf(vmo.samples_info_store, key='samples')
            df.iloc[spl_idx_list].to_hdf(samples_info_store,
                                         key='samples', mode='w', format='table')

    return sub_vmo_dir


def get_mis_ref_alt_num(genotype_matrix):
    # Count the number of each element in the matrix
    counts = np.apply_along_axis(lambda x: np.histogram(
        x, bins=[-1.5, -0.5, 0.5, 1.5])[0], 1, genotype_matrix)
    counts = np.apply_along_axis(lambda x: np.sum(x), 2, counts)

    # Calculate the number of each type of element
    mis_num = counts[:, 0]
    ref_num = counts[:, 1]
    alt_num = counts[:, 2]

    return mis_num, ref_num, alt_num


def get_mis_ref_alt_num_chunk(k, chunk_size, genotype_matrix):
    k_end = min(k + chunk_size, genotype_matrix.shape[0])
    genotype_matrix_chunk = genotype_matrix[k:k_end]
    mis_num_chunk, ref_num_chunk, alt_num_chunk = get_mis_ref_alt_num(
        genotype_matrix_chunk)
    return mis_num_chunk, ref_num_chunk, alt_num_chunk


def get_mis_ref_alt_num_parallel(genotype_matrix, chunk_size=1000, n_jobs=8):
    # Split the genotype_matrix into chunks
    mis_num, ref_num, alt_num = [], [], []

    total_chunks = genotype_matrix.shape[0] // chunk_size + \
        (1 if genotype_matrix.shape[0] % chunk_size != 0 else 0)
    batch_size = n_jobs * 100
    processed_chunks = 0

    # Count the number of each element in the chunks in parallel
    for start_chunk in range(0, total_chunks, batch_size):
        end_chunk = min(start_chunk + batch_size, total_chunks)

        results = Parallel(n_jobs=n_jobs)(
            delayed(get_mis_ref_alt_num_chunk)(
                k * chunk_size, chunk_size, genotype_matrix)
            for k in range(start_chunk, end_chunk)
        )

        processed_chunks += len(results)

        # Unpack the results
        mis_num_chunk, ref_num_chunk, alt_num_chunk = zip(*results)

        # Concatenate the results
        mis_num.extend(np.concatenate(mis_num_chunk))
        ref_num.extend(np.concatenate(ref_num_chunk))
        alt_num.extend(np.concatenate(alt_num_chunk))

        print("Time: %s, processed %d/%d chunks, %.2f%%" % (time.strftime(
            "%Y-%m-%d %H:%M:%S", time.localtime()), processed_chunks, total_chunks, processed_chunks/total_chunks*100))

    return np.array(mis_num), np.array(ref_num), np.array(alt_num)


def get_valid_variant_index(genotype_matrix, maf_thr=0.05, mis_thr=0.5, chunk_size=200, n_jobs=8):
    # # Count the number of each element in the matrix
    # counts = np.apply_along_axis(lambda x: np.histogram(
    #     x, bins=[-1.5, -0.5, 0.5, 1.5])[0], 1, genotype_matrix)
    # counts = np.apply_along_axis(lambda x: np.sum(x), 2, counts)

    # # Calculate the number of each type of element
    # mis_num = counts[:, 0]
    # ref_num = counts[:, 1]
    # alt_num = counts[:, 2]

    mis_num, ref_num, alt_num = get_mis_ref_alt_num_parallel(
        genotype_matrix, chunk_size=chunk_size, n_jobs=n_jobs)

    # Calculate the frequencies
    all_allele_num = genotype_matrix.shape[1]*genotype_matrix.shape[2]
    minor_freq = np.minimum(alt_num, ref_num) / (all_allele_num - mis_num)
    missing_freq = mis_num / all_allele_num

    # Determine which variants are valid
    valid_variant_bool_index = (
        minor_freq >= maf_thr) & (missing_freq <= mis_thr)

    return valid_variant_bool_index

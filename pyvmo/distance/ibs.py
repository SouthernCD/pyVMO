from itertools import combinations
from joblib import Parallel, delayed
import numpy as np
import time


def identity_by_state(s1, s2):
    return (1 - np.abs(s1 - s2) / 2).mean()


def get_IBS_matrix_loop(genotype_matrix, logger=None):
    if logger:
        logger.info("Start to calculate IBS matrix")
    genotype_matrix = np.sum(genotype_matrix, axis=2)
    if logger:
        logger.info("Finished to sum genotype matrix")
    results = np.zeros((genotype_matrix.shape[1], genotype_matrix.shape[1]))
    num = 0
    total_num = len(list(combinations(range(genotype_matrix.shape[1]), 2)))
    for i in range(genotype_matrix.shape[1]):
        gi = genotype_matrix[:, i]
        for j in range(genotype_matrix.shape[1]):
            if i > j:
                gj = genotype_matrix[:, j]
                results[i, j] = identity_by_state(gi, gj)
                results[j, i] = results[i, j]
                num += 1
                # print(num)
                if num % 100 == 0:
                    if logger:
                        logger.info("Time: %s, processed %d/%d pairs, %.2f%%" % (time.strftime(
                            "%Y-%m-%d %H:%M:%S", time.localtime()), num, total_num, num/total_num*100))
            if i == j:
                results[i, j] = 1.0
    return results


def parallel_ibs_computation(genotype_matrix, pair):
    i, j = pair
    gi = genotype_matrix[:, i]
    gj = genotype_matrix[:, j]
    return i, j, identity_by_state(gi, gj)


def get_IBS_matrix_parallel(genotype_matrix, num_threads=-1):
    genotype_matrix = np.sum(genotype_matrix, axis=2)

    num_individuals = genotype_matrix.shape[1]
    pairs = list(combinations(range(num_individuals), 2))

    # 使用 Joblib 进行并行计算
    results_list = Parallel(n_jobs=num_threads)(
        delayed(parallel_ibs_computation)(genotype_matrix, pair) for pair in pairs)

    results = np.zeros((num_individuals, num_individuals))
    for i, j, value in results_list:
        results[i, j] = value
        results[j, i] = value

    # 对角线上的值设置为 1
    np.fill_diagonal(results, 1.0)

    return results


def get_IBS_matrix_broadcasting(genotype_matrix):
    genotype_matrix = np.sum(genotype_matrix, axis=2)
    n_samples = genotype_matrix.shape[1]

    # 创建一个三维数组，其中第一维表示样本对，第二维和第三维表示每个样本的基因型
    genotype_pairs = np.abs(
        genotype_matrix[:, :, None] - genotype_matrix[:, None, :])

    # 一次性计算所有样本对的 IBS 值
    results = (1 - genotype_pairs / 2).mean(axis=0)

    return results


def get_IBS_matrix_broadcasting_chunk(genotype_matrix, chunk_size=1000):
    # 首先对genotype_matrix进行求和处理以简化数据

    genotype_matrix = np.sum(genotype_matrix, axis=2)
    n_samples = genotype_matrix.shape[1]

    results = np.zeros((n_samples, n_samples))

    # 对第一维度也进行分块处理
    chunk_number = 0
    for k in range(0, genotype_matrix.shape[0], chunk_size):
        k_end = min(k + chunk_size, genotype_matrix.shape[0])

        # 处理当前块的数据
        genotype_pairs = np.abs(
            genotype_matrix[k:k_end, :, None] -
            genotype_matrix[k:k_end, None, :]
        )

        # 因为是分块处理，所以这里需要累加结果，而不是直接赋值
        block_result = 1 - genotype_pairs / 2
        if k == 0:
            results = block_result
        else:
            results += block_result

        chunk_number += 1

    results = np.mean(results, axis=0)/chunk_number

    return results


def process_chunk_for_IBS(k, chunk_size, genotype_matrix):
    """
    处理单个数据块的函数。
    """
    k_end = min(k + chunk_size, genotype_matrix.shape[0])
    genotype_pairs = np.abs(
        genotype_matrix[k:k_end, :, None] - genotype_matrix[k:k_end, None, :]
    )
    block_result = 1 - genotype_pairs / 2
    return block_result.sum(axis=0)  # 返回当前块处理后的结果之和，减少内存使用


def get_IBS_matrix_broadcasting_chunk_parallel(genotype_matrix, chunk_size=100, n_jobs=8):
    """
    优化内存使用的并行处理计算IBS矩阵。
    """
    # 首先对genotype_matrix进行求和处理以简化数据
    genotype_matrix = np.sum(genotype_matrix, axis=2)
    n_samples = genotype_matrix.shape[1]

    # 初始化结果矩阵
    results = np.zeros((n_samples, n_samples))

    total_chunks = genotype_matrix.shape[0] // chunk_size + \
        (1 if genotype_matrix.shape[0] % chunk_size != 0 else 0)
    batch_size = n_jobs * 100  # 每次处理20个chunk
    processed_chunks = 0

    for start_chunk in range(0, total_chunks, batch_size):
        end_chunk = min(start_chunk + batch_size, total_chunks)
        all_results = Parallel(n_jobs=n_jobs)(
            delayed(process_chunk_for_IBS)(k * chunk_size, chunk_size, genotype_matrix)
            for k in range(start_chunk, end_chunk)
        )

        # 将这批次的结果累加到最终结果中
        for block_result in all_results:
            results += block_result

        processed_chunks += len(all_results)
        print("Time: %s, processed %d/%d chunks, %.2f%%" % (time.strftime(
            "%Y-%m-%d %H:%M:%S", time.localtime()), processed_chunks, total_chunks, processed_chunks/total_chunks*100))

    # 计算最终结果
    results /= genotype_matrix.shape[0]

    return results


def process_chunk_for_IBS_loop(k, chunk_size, genotype_matrix):
    """
    处理单个数据块的函数。
    """
    k_end = min(k + chunk_size, genotype_matrix.shape[0])
    m = genotype_matrix[k:k_end]
    n = m.shape[1]
    results = np.zeros((n, n))
    for i in range(n):
        gi = m[:, i]
        for j in range(n):
            if i > j:
                gj = m[:, j]
                results[i, j] += identity_by_state(gi, gj)
                results[j, i] = results[i, j]
            if i == j:
                results[i, j] += 1.0

    return results

def get_IBS_matrix_loop_chunk_parallel(genotype_matrix, chunk_size=200, n_jobs=80):
    """
    优化内存使用的并行处理计算IBS矩阵。
    """
    # 首先对genotype_matrix进行求和处理以简化数据
    genotype_matrix = np.sum(genotype_matrix, axis=2)
    n_samples = genotype_matrix.shape[1]

    # 初始化结果矩阵
    results = np.zeros((n_samples, n_samples))

    total_chunks = genotype_matrix.shape[0] // chunk_size + \
        (1 if genotype_matrix.shape[0] % chunk_size != 0 else 0)
    batch_size = n_jobs * 100  # 每次处理20个chunk
    processed_chunks = 0
    
    for start_chunk in range(0, total_chunks, batch_size):
        
        end_chunk = min(start_chunk + batch_size, total_chunks)
        all_results = Parallel(n_jobs=n_jobs)(
            delayed(process_chunk_for_IBS_loop)(k * chunk_size, chunk_size, genotype_matrix)
            for k in range(start_chunk, end_chunk)
        )

        # 将这批次的结果累加到最终结果中
        for block_result in all_results:
            results += block_result

        processed_chunks += len(all_results)
        print("Time: %s, processed %d/%d chunks, %.2f%%" % (time.strftime(
            "%Y-%m-%d %H:%M:%S", time.localtime()), processed_chunks, total_chunks, processed_chunks/total_chunks*100))

    # 计算最终结果
    results /= total_chunks

    return results

# dask
# import numpy as np
# import dask.array as da
# from dask import delayed, compute
# from itertools import combinations
# from dask.distributed import Client


# def get_IBS_matrix_dask(genotype_matrix, chunk_size=200, n_jobs=80):
#     client = Client(n_workers=1, threads_per_worker=n_jobs)  # 举例：启动1个worker，每个worker使用4个线程

#     genotype_matrix = np.sum(genotype_matrix, axis=2)
#     y = genotype_matrix.shape[1]

#     # 将NumPy memmap对象转换为Dask数组
#     # 注意：这里的chunks参数需要根据您的具体需求调整，以达到最佳性能
#     data = da.from_array(genotype_matrix, chunks=(chunk_size, y))

#     # 使用Dask计算每对列之间的自定义函数结果
#     results = []

#     for i, j in combinations(range(y), 2):
#         col1, col2 = data[:, i], data[:, j]
#         result = delayed(identity_by_state)(col1, col2)
#         results.append(result)

#     # 并行计算所有结果
#     distance_matrix = np.zeros((y, y))

#     # 计算结果
#     computed_distances = compute(*results)

#     # 填充距离矩阵
#     for k, (i, j) in enumerate(combinations(range(y), 2)):
#         distance_matrix[i, j] = computed_distances[k]
#         distance_matrix[j, i] = computed_distances[k]  # 因为是距离矩阵，所以它是对称的

#     # 对称轴上的值设置为0
#     np.fill_diagonal(distance_matrix, 1)

#     return distance_matrix

if __name__ == '__main__':
    from pyvmo import VMO

    vmo_dir = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/1.raw_data/beagle/reseq_landrace_passed_vmo"
    vmo = VMO(vmo_dir)
    m = vmo.get_matrix()

    # # 
    # t = m[:100000,:200]
    # %time ibs_m = get_IBS_matrix_loop(t)
    # %time ibs_m = get_IBS_matrix_parallel(t, 8)
    # %time ibs_m = get_IBS_matrix_broadcasting(t)
    # %time ibs_m = get_IBS_matrix_broadcasting_chunk(t, 1000)
    # %time ibs_m = get_IBS_matrix_broadcasting_chunk_parallel(t,1000,8)

    # # 
    # t = m[:10000,:200]
    # %time ibs_m = get_IBS_matrix_loop(t)
    # %time ibs_m = get_IBS_matrix_parallel(t, 8)
    # %time ibs_m = get_IBS_matrix_broadcasting_chunk_parallel(t,1000,8)
    # %time ibs_m = get_IBS_matrix_loop_chunk_parallel(t,100,8)
    # %time ibs_m = get_IBS_matrix_broadcasting_chunk_parallel(m,200,80)

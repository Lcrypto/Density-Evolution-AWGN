# Illustration of how multiprocessing will be used in algorithms
import time
from multiprocessing import shared_memory, Pool
import numpy as np
from random import random


def random_num(i):
    global shr_name, dim
    existing_shm = shared_memory.SharedMemory(name=shr_name)
    np_array = np.ndarray(dim, dtype=np.float64, buffer=existing_shm.buf)

    np_array[i] = random()

    existing_shm.close()

def increment(i):
    global shr_name, dim
    existing_shm = shared_memory.SharedMemory(name=shr_name)
    np_arr = np.ndarray(dim, dtype=np.float64, buffer=existing_shm.buf)

    np_arr[i] += 1

    existing_shm.close()


# Use initializer to share data that is fixed throughout computations
def initializer():
    global dim, shr_name


if __name__ == '__main__':
    global dim, shr_name
    dim = 10
    n_processes = 4

    # Use shared_memory to share data between processes that operate on the same array
    # but never on the same part of the array. In this case we are interested in sharing
    # numpy arrays
    a = np.ones(dim, dtype=np.float64)
    shm = shared_memory.SharedMemory(create=True, size=a.nbytes)
    # # Now create a NumPy array backed by shared memory
    np_array = np.ndarray(a.shape, dtype=np.float64, buffer=shm.buf)
    np_array[:] = a[:]  # Copy the original data into shared memory
    shr_name = shm.name
    print(np_array)

    # Use pool to split up processes
    pool = Pool(processes=n_processes, initializer=initializer)
    idx = range(dim)
    pool.map(random_num, idx)
    pool.map(increment, idx)

    print(np_array)
    shm.close()
    shm.unlink()

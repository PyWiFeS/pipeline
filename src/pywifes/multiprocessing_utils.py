import multiprocessing
import os
import math

# Consider using 
# with start_pool(max_processes=max_processes) as pool:
#     results = run_tasks(tasks, pool)

def start_pool(max_processes=-1):
    '''
    Start a multiprocessing pool with as many workers as are likely to be efficient for parallel tasks.

    If `max_processes` is less than 1, it will simply use os.cpu_count(),
    the number of hardware & logical cores.

    However, this value can differ from the number of *available* cores.
    e.g. for Slurm users.

    In such cases, `max_processes` is used as an upper limit on the number of workers.

    ***NOTE:*** The pool should be closed with `pool.close()` once you have finished your
    parallel tasks to free up hardware resources.
    '''
    num_processes = os.cpu_count()
    if max_processes > 0:
        num_processes = min(max_processes, num_processes)

    return multiprocessing.Pool(num_processes)


def get_task(func, *args, **kwargs):
    '''
    Convert the arguments provided into a 'task' tuple, suitable for `run_tasks`
    '''
    return (func, args, kwargs)


def _unwrap_and_run(task):
    func, args, kwargs = task
    return func(*args, **kwargs)


def run_tasks_singlethreaded(tasks):
    '''
    Run the `tasks` in a single thread.

    Each task should follow the pattern in `get_task`, storing the function, args and kwargs to run.

    The results will be returned in order.
    '''
    return [_unwrap_and_run(task) for task in tasks]


def run_tasks(tasks, pool, chunksize=-1):
    '''
    Run the `tasks` using the pool.

    Each task should follow the pattern in `get_task`, storing the function, args and kwargs to run.

    The results will be returned in order.
    '''
    if chunksize < 1:
        # By default, divide tasks evenly between processes.
        # For large numbers of tasks, consider setting a smaller chunksize.
        chunksize = math.ceil(len(tasks) / pool._processes)

    # `lazy_results` must be enumerated within the pool's scope or the threads will not complete.
    # To avoid confusion during development, this is done within this function.
    lazy_results = pool.imap(_unwrap_and_run, tasks, chunksize=chunksize)
    return list(lazy_results)

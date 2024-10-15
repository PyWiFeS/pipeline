import math
import multiprocessing
import os


def get_num_processes(max_processes=-1):
    """
    This function returns the number of processes that are likely to be efficient for
    parallel tasks.

    If `max_processes` is less than 1, it will simply return os.cpu_count(),
    the number of hardware & logical cores.

    However, this value can differ from the number of *available* cores.
    e.g. for Slurm users.

    In such cases, `max_processes` is used as an upper limit on the returned value.
    """
    num_processes = os.cpu_count()
    if max_processes > 0:
        num_processes = min(max_processes, num_processes)
    return num_processes


def _unwrap_and_run(task):
    func, args, kwargs = task
    return func(*args, **kwargs)


def map_tasks(tasks, max_processes=-1, chunksize=-1):
    """
    Run the `tasks`, divided between up to `max_processes` processes in chunks of
    `chunksize`.

    Each task should follow the pattern in `get_task`, storing the function, args and
    kwargs to run.

    The results will be returned in order.

    If `chunksize` is less than 1, then tasks will be run in a single batch.
    """
    num_processes = get_num_processes(max_processes)

    if chunksize < 1:
        # By default, divide tasks evenly between processes.
        # For extremely large numbers of tasks, consider setting a smaller chunksize.
        chunksize = math.ceil(len(tasks) / num_processes)

    results = []
    with multiprocessing.Pool(num_processes) as pool:
        # `lazy_results` must be enumerated within the pool's scope or the threads will
        # not complete.
        lazy_results = pool.imap(_unwrap_and_run, tasks, chunksize=chunksize)
        results = list(lazy_results)

    return results


def get_task(func, *args, **kwargs):
    """
    Convert the arguments provided into a 'task' tuple, suitable for `map_tasks`
    """
    return (func, args, kwargs)

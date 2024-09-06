import datetime
import logging
from functools import wraps


def timeit(func):
    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        start_time = datetime.datetime.now()
        result = func(*args, **kwargs)
        end_time = datetime.datetime.now()
        total_time = str(end_time - start_time).split(".")[0]
        logging.info(f"Total execution time {total_time}")
        return result

    return timeit_wrapper

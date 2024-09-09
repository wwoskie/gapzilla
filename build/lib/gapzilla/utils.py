"""
Utils for time measurement and args splitting
"""

import datetime
import logging
from functools import wraps
from typing import Callable, Any


def timeit(func: Callable) -> Callable:
    """
    A decorator that logs the execution time of the wrapped function.

    Parameters
    ----------
    func : Callable
        The function to be wrapped and timed.

    Returns
    -------
    Callable
        The wrapped function with added timing and logging functionality.

    Example
    -------
    >>> @timeit
    ... def example_function():
    ...     for _ in range(1000000):
    ...         pass
    >>> example_function()
    Total execution time 0:00:00
    """

    @wraps(func)
    def timeit_wrapper(*args: Any, **kwargs: Any) -> Any:
        start_time = datetime.datetime.now()
        result = func(*args, *kwargs)
        end_time = datetime.datetime.now()
        total_time = str(end_time - start_time).split(".")[0]
        logging.info(f"Total execution time {total_time}")
        return result

    return timeit_wrapper


def list_of_strings(arg: str) -> list[str]:
    """
    Function to split string into list of strings

    Parameters
    ----------
    arg : str
        The input string to be split.

    Returns
    -------
    list of str
        A list of substrings obtained by splitting the input string at each ", " delimiter.
    """
    return arg.split(", ")

"""Timeout decorator for pipeline function execution limits."""
from __future__ import annotations
import os
import errno
import signal
import functools


# pylint: disable=unused-argument
def timeout(seconds: int = 10, error_message: str = os.strerror(errno.ETIME)):
    """Decorator to raise a TimeoutError if the process runs over time.

    Args:
        seconds: Default timeout in seconds. Can be overridden by passing
                 'timeout' kwarg to the decorated function.
        error_message: Error message for the TimeoutError.

    Example:
        @timeout(seconds=30)
        def slow_function(x, timeout=None):
            # timeout kwarg will override the decorator default
            ...

    Note:
        Uses signal.SIGALRM, so this only works on Unix-like systems and
        only in the main thread.
    """
    def decorator(func):
        def _handle_timeout(signum, frame):
            raise TimeoutError(error_message)

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            timeout_seconds = kwargs.get('timeout', seconds)
            if timeout_seconds:
                signal.signal(signal.SIGALRM, _handle_timeout)
                signal.alarm(timeout_seconds)
            try:
                result = func(*args, **kwargs)
            finally:
                if timeout_seconds:
                    signal.alarm(0)
            return result

        return wrapper

    return decorator

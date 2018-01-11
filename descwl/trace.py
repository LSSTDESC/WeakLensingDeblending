"""Trace program resource usage.
"""
from __future__ import print_function
import os

class Memory(object):
    """Trace memory usage for the current program.

    Args:
        enabled(bool): Enable memory tracing.
    """
    def __init__(self,enabled):
        self.enabled = enabled
        if self.enabled:
            # Defer the psutil import to here so that it does not need
            # to be installed unless tracing is enabled.
            try:
                import psutil
            except ImportError:
                raise RuntimeError(
                    'Missing required psutil import for memory tracing.')
            self.this_process = psutil.Process(os.getpid())
            self.last_usage = 0

    def __call__(self,label):
        """Register a memory usage checkpoint.

        This method does nothing if this object was initialized with
        enabled = False.

        Args:
            label(str): A brief description of this checkpoint.
        """
        if not self.enabled:
            return
        usage = self.this_process.get_memory_info()[0]
        print('%s memory usage: %.3f Mb (%+d bytes)' % (label,
            usage/float(2**20),usage-self.last_usage))
        self.last_usage = usage

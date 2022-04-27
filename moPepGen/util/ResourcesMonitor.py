""" ResourcesMonitor monitors the resources usage including CPU and memory """
from __future__ import annotations
import datetime
import math
import psutil


class MemoryUnit():
    """ Memory Unit """
    _suffixes = ['bytes', 'KiB', 'MiB', 'GiB', 'TiB', 'PiB', 'EiB', 'ZiB', 'YiB']
    def __init__(self, size:int):
        """ constructor """
        self._data = size

    def __str__(self) -> str:
        """ str """
        return self.human_readable(self._data)

    def human_readable(self, x:int) -> str:
        """ Get humand readable representation of the memroy unit. """
        # https://stackoverflow.com/a/25613067/11081630
        order = int(math.log2(x) / 10) if x else 0
        val = x/(1 << (order * 10))
        unit = self._suffixes[order]
        return f"{val:.4g} {unit}"

class ResourcesMonitor():
    """ Monitor the resources usage of the current process. """
    def __init__(self):
        """ constructor """
        self.process = psutil.Process()
        # According to the documentation, this function returns the cpu percent
        # since the last call. So calling it at the beginning.
        self.process.cpu_percent()

    def get_resource_usage(self) -> ResourcesUsage:
        """ Get resources usage """
        mem = self.process.memory_info()
        rss = MemoryUnit(mem.rss)
        vms = MemoryUnit(mem.vms)

        cpu_times = self.process.cpu_times()
        wallclock = datetime.timedelta(seconds=cpu_times.user)
        system = datetime.timedelta(seconds=cpu_times.system)

        cpu_percent = self.process.cpu_percent()

        return ResourcesUsage(
            rss=rss, vms=vms, wallclock=wallclock, system=system,
            cpu_percent=cpu_percent
        )

class ResourcesUsage():
    """ Resources usage """
    def __init__(self, rss:MemoryUnit=None, vms:MemoryUnit=None,
            wallclock:datetime.timedelta=None, system:datetime.timedelta=None,
            cpu_percent:float=None):
        """ """
        self.rss = rss
        self.vms = vms
        self.wallclock = wallclock
        self.system = system
        self.cpu_percent = cpu_percent

    def __str__(self) -> str:
        """ str """
        return  f"rss={self.rss}, vms={self.vms}, wallclock={self.wallclock}, " +\
            f"system={self.system}, cpu_usage={self.cpu_percent}%"

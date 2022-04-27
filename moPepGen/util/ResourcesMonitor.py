""" ResourcesMonitor monitors the resources usage including CPU and memory """
import datetime
import math
import psutil


class MemoryUnit():
    """ Memory Unit """
    _suffixes = ['bytes', 'KiB', 'MiB', 'GiB', 'TiB', 'PiB', 'EiB', 'ZiB', 'YiB']
    def __init__(self, size:int):
        """ """
        self._data = size

    def __str__(self) -> str:
        """ """
        return self.human_readable(self._data)

    def human_readable(self, x:int) -> str:
        """ """
        # https://stackoverflow.com/a/25613067/11081630
        order = int(math.log2(x) / 10) if x else 0
        val = x/(1 << (order * 10))
        unit = self._suffixes[order]
        return f"{val:.4g} {unit}"

class ResourcesMonitor():
    """ """
    def __init__(self):
        """ """
        self.process = psutil.Process()
        # According to the documentation, this function returns the cpu percent
        # since the last call. So calling it at the beginning.
        self.process.cpu_percent()
        self.rss:MemoryUnit = None
        self.vms:MemoryUnit = None
        self.wallclock:datetime.timedelta = None
        self.system:datetime.timedelta = None
        self.cpu_percent:float = None

    def get_resource_usage(self):
        """ Get resources rsage """
        mem = self.process.memory_info()
        self.rss = MemoryUnit(mem.rss)
        self.vms = MemoryUnit(mem.vms)

        cpu_times = self.process.cpu_times()
        self.wallclock = datetime.timedelta(seconds=cpu_times.user)
        self.system = datetime.timedelta(seconds=cpu_times.system)

        self.cpu_percent = self.process.cpu_percent()

    def __str__(self) -> str:
        return  f"rss={self.rss}, vms={self.vms}, wallclock={self.wallclock}, " +\
            f"system={self.system}, cpu_usage={self.cpu_percent}%"

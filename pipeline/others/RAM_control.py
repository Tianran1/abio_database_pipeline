"""
script to control memory taken up by jupyter notebooks,

paramerters to configure: 
    see buttom, if __name__ == '__main__':
    MEMORY_CUTOFF = 30 # in GB
    INTERVAL = 120 # in seconds

how this work:
    psutil is a python package monitoring system resource,
    this code is a simply a wrapper of psutil functions, customized for our need,

    every time INTERVAL, perform a check on total memory usage,
    once reaches the MEMORY_CUTOFF, kill from the largest process, until memory if below MEMORY_CUTOFF

how to use:

    nohup python3 -u RAM_control.py > RAM_control.log &
"""


import os
import signal
import time
from datetime import datetime
from typing import Tuple
from pprint import pprint

import psutil


def get_memory_usage(cutoff :int = 100) -> Tuple[bool, float]:
    """
    get total memory usage, 
    cutoff is integer indicating GB of memory usage
    """
    
    overflow = False
    total_memory = psutil.virtual_memory().used / 1014 / 1024 / 1024 # in GB
    if total_memory > cutoff:
        overflow = True
        
    return overflow, total_memory


def get_large_processes(percent_memory :int = 1) -> list:
    """
    get large jupyter notebook processes excedding percent_memory percent of memory
    returns a dictionary with pid as keys.
    
    e.g.:
    
    get_large_processes(percent_memory = 1) returns:
    
    {19599: '--port=9024 taking up 2.968020787359343 percent memory, total memory '
            'is 128GB',
     21532: '--port=9019 taking up 5.857101495118028 percent memory, total memory '
            'is 128GB',
     21996: '--port=9018 taking up 12.677567467218703 percent memory, total memory '
            'is 128GB',
     22733: '--port=9003 taking up 4.495330484205427 percent memory, total memory '
            'is 128GB',
     24068: '--port=9014 taking up 1.5671722053595774 percent memory, total memory '
            'is 128GB'}
    """
    
    # select parent jupyter processes from all processes
    pids = psutil.pids()
    jupyter_processes = [psutil.Process(pid) for pid in pids if psutil.Process(pid).name() == 'jupyter-noteboo']

    # select large parent jupyter processes
    large_processes = []
    for jupyter_process in jupyter_processes:
        children_processes = jupyter_process.children()
        for children_process in children_processes:
            if children_process.memory_percent() > percent_memory:
                pid = jupyter_process.pid
                pid_memory_percent = children_process.memory_percent()
                pid_info = '%s taking up %s percent memory, total memory is 128GB'%(
                                jupyter_process.as_dict()['cmdline'][-4:], pid_memory_percent)
                large_processes.append((pid, pid_memory_percent, pid_info))
    
    # sort large parent jupyter processes by memory usage
    large_processes = sorted(large_processes, key=lambda p: p[1], reverse=True)
    
    return large_processes


# taken from https://psutil.readthedocs.io/en/latest/#recipes
def kill_proc_tree(pid, sig=signal.SIGTERM, include_parent=True,
                   timeout=None, on_terminate=None):
    """
    Kill a process tree (including grandchildren) with signal
    "sig" and return a (gone, still_alive) tuple.
    "on_terminate", if specified, is a callabck function which is
    called as soon as a child terminates.
    """
    parent = psutil.Process(pid)
    children = parent.children(recursive=True)
    if include_parent:
        children.append(parent)
    for p in children:
        p.send_signal(sig)
    gone, alive = psutil.wait_procs(children, timeout=timeout,
                                    callback=on_terminate)
    return (gone, alive)
      


if __name__ == '__main__':
    
    MEMORY_CUTOFF = 105 # in GB
    INTERVAL = 10 # in seconds

    while True:
        time.sleep(INTERVAL)
        try:
            overflow, total_memory = get_memory_usage(cutoff = MEMORY_CUTOFF)
            
            if overflow:
                print('\n'*3, datetime.now(), 'encountered an overflow')
                print('overflow: ', overflow, 'total memory in GB: ', total_memory)
                large_processes = get_large_processes()
                pprint(large_processes)
                
                for p in large_processes:
                    if overflow:
                        kill_proc_tree(p[0])
                        print('\n', datetime.now())
                        print('killed a jupyter notebook to release memory, info: %s' %p[2])
                        overflow, total_memory = get_memory_usage(cutoff=MEMORY_CUTOFF)
                    else:
                        print('memory goes back to normal')
                        break
                        
        except Exception as e:
            print('encountered exception %s' %e)

import threading

class LockedIterator(object):
    '''make generator/iterator thread safe'''
    def __init__(self, it):
        self.lock = threading.Lock()
        self.it = iter(it)

    def __iter__(self): 
        return self

    def __next__(self):
        with self.lock:
            return self.it.__next__()

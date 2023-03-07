import threading

class LockedIterator:
    '''
    make generator/iterator thread safe
    from https://gist.github.com/platdrag/e755f3947552804c42633a99ffd325d4
    '''
    def __init__(self, it):
        self.lock = threading.Lock()
        #self.it = iter(it)
        self.it = it

    def __iter__(self): 
        return self

    def __next__(self):
        with self.lock:
            return next(self.it)


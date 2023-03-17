import threading
from multiprocessing import Queue

class LockedIterator:
    '''
    make generator/iterator thread safe
    from https://gist.github.com/platdrag/e755f3947552804c42633a99ffd325d4
    '''

    def __init__(self, it):
        self.it = iter(it) #otherwise it's not iterable
        self.queue = Queue()
        self.thread = threading.Thread(target=self._run)
        self.thread.start()

    def _run(self):
        for item in self.it:
            self.queue.put(item)
        self.queue.put(None)
    
    def __iter__(self): 
        return self

    def __next__(self):
        item = self.queue.get()
        if item is None:
            raise StopIteration
        return item


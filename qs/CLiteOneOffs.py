import queue

from .CLiteDbThread import CLiteDbThread
from .CLiteHTTPThread import CLiteHTTPThread


def cLiteDbThread(*argv, **kw):
    if CLiteDbThread.insts is None:
        hasStartedQueue = queue.Queue()
        kw['hasStartedQueue'] = hasStartedQueue
        newCLiteDbThread = CLiteDbThread(*argv, **kw)
        newCLiteDbThread.start()
        print(hasStartedQueue.get(True, 10))
    return CLiteDbThread.insts


def cLiteHTTPThread(*argv, **kw):
    if CLiteHTTPThread.insts is None:
        hasStartedQueue = queue.Queue()
        kw['hasStartedQueue'] = hasStartedQueue
        newCLiteHTTPThread = CLiteHTTPThread(*argv, **kw)
        newCLiteHTTPThread.start()
        print(hasStartedQueue.get(True, 10))
    return CLiteHTTPThread.insts

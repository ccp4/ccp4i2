from __future__ import print_function

import queue

def cLiteDbThread(*argv, **kw):
    from CLiteDbThread import CLiteDbThread
    if CLiteDbThread.insts is None:
        hasStartedQueue = queue.Queue()
        kw['hasStartedQueue'] = hasStartedQueue
        newCLiteDbThread = CLiteDbThread(*argv, **kw)
        newCLiteDbThread.start()
        print(hasStartedQueue.get(True, 10))
    return CLiteDbThread.insts

def cLiteHTTPThread(*argv, **kw):
    from CLiteHTTPThread import CLiteHTTPThread
    if CLiteHTTPThread.insts is None:
        hasStartedQueue = queue.Queue()
        kw['hasStartedQueue'] = hasStartedQueue
        newCLiteHTTPThread = CLiteHTTPThread(*argv, **kw)
        newCLiteHTTPThread.start()
        print(hasStartedQueue.get(True, 10))
    return CLiteHTTPThread.insts


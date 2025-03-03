from __future__ import print_function

import sys
import multiprocessing
import os
import traceback


def QuarantinedProcessing_decorate(func):
    def func_wrapper(*args, **kwargs):
        try:
            resultQueue = args[0]
            result = func(*args, **kwargs)
            resultQueue.put(result)
        except Exception as err:
            traceback.print_exc(10)
            print("Process failed with exception ", err)
            resultQueue.put("FailedProcess")
    return func_wrapper

@QuarantinedProcessing_decorate
def remakeReport(resultQueue, jobId, jobStatus, dbFileName):
    print(jobId, jobStatus, dbFileName)
    if jobId is not None:
        from core import CCP4ProjectsManager
        from utils.startup import startProjectsManager
        pm = startProjectsManager(dbFileName=dbFileName)
        from report import CCP4ReportGenerator
        generator = CCP4ReportGenerator.CReportGenerator(jobId=jobId, jobStatus=jobStatus)
        if jobStatus == 'Failed':
            reportFile = generator.makeFailedReportFile(redo=True)
            pm.db().close()
            CCP4ProjectsManager.CProjectsManager.insts = None
            return reportFile
        reportFile, newPageOrNewData = generator.makeReportFile(redo=True, makePictures=False)
        resultQueue.put((reportFile, newPageOrNewData))
        if hasattr(generator,"report") and generator.report.containsPictures():
            for picture in generator.report.pictureQueue:
                pass
        pm.db().close()
        CCP4ProjectsManager.CProjectsManager.insts = None
        return "Completed"

if __name__ == "__main__":
    print('sys.argv', sys.argv)
    if sys.argv[1] == "remakeReport":
        resultQueue = multiprocessing.Queue()
        reportProcess = multiprocessing.Process(target=remakeReport, args=(resultQueue,sys.argv[2],sys.argv[3], sys.argv[4],))
        reportProcess.start()
        print(resultQueue.get())
        reportProcess.join()

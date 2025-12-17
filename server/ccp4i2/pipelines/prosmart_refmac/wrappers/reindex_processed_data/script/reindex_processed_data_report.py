from ccp4i2.report import Report


class reindex_processed_data_report(Report):
    TASKNAME = 'reindex_processed_data'
    RUNNING = False
    USEPROGRAMXML = False
    
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)

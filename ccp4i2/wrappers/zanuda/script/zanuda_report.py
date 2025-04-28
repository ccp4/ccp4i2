import os
import re

from ....report.CCP4ReportParser import Report


class zanuda_report(Report):
    TASKNAME = 'zanuda'
    RUNNING = True
    USEPROGRAMXML = False
    WATCHED_FILE = 'log.txt'

    use_folds = False  # folds do not look good

    def __init__(self, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        self.addDiv(style='clear:both;')
        if jobStatus in ['Running', 'Running remotely']:
            self.append('<p><b>The job is currently running.</b></p>')
        log_path = os.path.join(jobInfo['fileroot'], self.WATCHED_FILE)
        if os.path.exists(log_path):
            with open(log_path) as istream:
                txt = istream.read()
                p1 = '\n *Step [1-3][.] *\n.*\n.*\n\n'
                p2 = '^ *[|] *([<> ]{2}) *([0-9]+) *[|]'
                p2 += ' *([A-Z][0-9 ]*[0-9]) *[|]'
                p2 += 4* ' *(--|[0-9.]+|error) *[|]' + ' *$'
                p3 = '\n *-+ *\n([^|]*)end of job:'
                hh = [mo.span() for mo in re.finditer(p1, txt)]
                if hh:
                    hhz = list(zip(*hh))
                    bb = list(zip(hhz[1], hhz[0][1:] + (len(txt),)))
                    b = bb[-1]
                    mobj = re.search(p3, txt[b[0]:b[1]], re.M)
                    if mobj:
                        self.append(' '.join(mobj.group(1).split()))
                    oo = (len(hh) - 1)* [False] + [True]
                    for h, b, o in zip(hh, bb, oo):
                        msg = ' '.join(txt[h[0]:h[1]].strip().split())
                        if self.use_folds:
                            label, sep, msg = msg.partition('. ')
                            fold = self.addFold(label = label, initiallyOpen = o)
                        else:
                            fold = self

                        fold.append(msg)
                        rr = re.findall(p2, txt[b[0]:b[1]], re.M)
                        if rr:
                            table = fold.addTable(transpose=False)
                            rrz = list(zip(*rr))
                            table.addData(title='In/Out', data=rrz[0])
                            table.addData(title='Subgroup No', data=rrz[1])
                            table.addData(title='Spacegroup', data=rrz[2])
                            table.addData(title='rmsd', data=rrz[3])
                            table.addData(title='R-work (RB)', data=rrz[4])
                            table.addData(title='R-work', data=rrz[5])
                            table.addData(title='R-free', data=rrz[6])

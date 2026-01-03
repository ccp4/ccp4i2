import os
import sys
import time


def getCCP4I2Dir(up=1):
    target = os.path.join(os.path.realpath(sys.argv[0]),"..")
    abstarget = os.path.abspath(target)
    splittarget = abstarget.split()
    if splittarget.count('ccp4i2'):
        splittarget.reverse()
        up = splittarget.index('ccp4i2')
    while up>0:
        abstarget = os.path.dirname(abstarget)
        up = up -1
    return abstarget

def getUserId():
    name = os.environ.get('LOGNAME',None)
    if name is not None:
        return name
    name = os.environ.get('USERNAME',None)
    if name is not None:
        return name
    import getpass
    name = getpass.getuser()
    if name is not None:
        return name
    try:
        return os.getlogin()
    except:
        return None


if __name__ == '__main__':

    import psutil
    containsList = ['ccp4']

    def contains(exe,containsList):
        if exe is None:
            return False
        for item in containsList:
            if exe.count(item):
                return True
        return False

    pInfoDict = {}
    me = getUserId()
    for proc in psutil.process_iter():
        try:
            pinfo = proc.as_dict(attrs=['pid', 'name', 'username','exe','create_time'])
        except psutil.NoSuchProcess:
            pass
        else:
            if pinfo['username'] == me and (len(containsList)==0 or contains(pinfo['exe'],containsList)):
                try:
                    pinfo['parent'] = proc.parent().as_dict(attrs=['pid'])['pid']
                except:
                    pinfo['parent'] = -1
                pinfo['children'] = []
                for p in proc.children(recursive=True):
                    try:
                        pinfo['children'].append(p.pid)
                    except:
                        pass
                #print(pinfo)
                pInfoDict[pinfo['pid']] = pinfo
    print('processes=',pInfoDict)
    print('atTime='+time.time())



from __future__ import print_function


class CCP4I2_XXX_Interface:

  SHOW_USER_PARAMS = [ 'username' , 'password' ]
  def __init__(self,serverParams={}):
    print('CCP4I2_XXX_Interface')
    #  self.serverParams is a reference to the CJobServer.serverParams instance with data for this job run   
    self.serverParams = serverParams
    self.serverToken = None
    # This is complicated -- keep a list of the jobIds for jobs that are running
    # and must be polled to see if they have finished.  If the run mechanism has
    # build in mechanism to report job finish then this is unnecessary.. or it may be
    # necessary only after i2 is closed down and lost connections.
    self.pollList = []

  def openConnection(self,jobId,args):
    if self.serverToken is None:
      # Open connection and get server token
      pass

  def setup(self,jobId):
    # Setup the job by:
    # 1) transporting the tarball self.serverParams[jobId].local_tarball  to remote machine
    # 2) creating/transporting scripts
    # 3) The contents of self.serverParams[jobId] should be updated appropriately
    #   If the tarball path on remote machine is remote_path/projectname_jobnumber_setup.ccp4db.zip
    #   then need to set in self.serverParams[jobId] path to files that are going to be used later:
    #   sP.remote_finished_tarball=remote_path/projectname_jobnumber_finished.ccp4db.zip
    #   sP.remote_report_file = remote_path/projectname_jobnumber_work/project/CCP4_JOBS/job_jobNumber/program.xml'
    pass

  def submit(self,jobId):
    # Actually start the remote job
    pass

  def transportFiles(self,jobId,copyList=[],mode='put'):
    # Each item of copyList is list of: local_file,remote_file
    # Mode is 'put'/'get'
    pass
    
  def pollForFinishedJob(self):
    # For each job in self.pollList check if it has finished
    # Remove any finished job from the pollList
    # Return a list of jobId for jobs that have finished   
    return []

  def handleFinishedJob(self,jobId):
    # Called after remote job finished needs to:
    # 1) Retrieve remote tarball to self.serverParams[jobId].local_finished_tarball
    pass


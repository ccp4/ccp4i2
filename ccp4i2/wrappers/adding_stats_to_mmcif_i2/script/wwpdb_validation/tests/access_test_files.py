import os


class TestFiles:

    def __init__(self):
        self.FILE_ROOT = os.path.dirname(os.path.realpath(__file__))
        self.path = os.path.dirname(os.path.join(self.FILE_ROOT, '..', '..'))
        self.test_data = os.path.join(self.path, 'test_data')

        self.xml = None
        self.cif = None
        self.sf = None

    def pdb3zt9(self):
        self.xml = os.path.join(self.test_data, "3zt9_validation.xml")
        self.cif = os.path.join(self.test_data,'pdb3zt9_refmac1_output.cif')
        self.sf = os.path.join(self.test_data,'r3zt9sf.ent')

    def pdb1cbs(self):
        self.xml = os.path.join(self.test_data, "1cbs_validation.xml")

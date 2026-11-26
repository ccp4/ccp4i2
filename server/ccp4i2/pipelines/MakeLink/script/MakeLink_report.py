import os

from ccp4i2.pipelines.MakeLink.script.link_dictionary import describe_bond, link_modifications
from ccp4i2.report.CCP4ReportParser import Report


class MakeLink_report(Report):
    TASKNAME = 'MakeLink'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        clearingDiv = self.addDiv()
        self.addModifications(self)
        self.addDefaultReport(self)
        clearingDiv = self.addDiv()

    def linkDictionaryPath(self):
        """The link dictionary this job wrote, or None."""
        fileroot = self.jobInfo.get('fileroot') if self.jobInfo else None
        if not fileroot:
            return None
        for outfile in self.jobInfo.get('outputfiles') or []:
            if outfile.get('jobparamname') == 'CIF_OUT' and outfile.get('filename'):
                path = os.path.join(fileroot, outfile['filename'])
                if os.path.isfile(path):
                    return path
        try:
            names = sorted(name for name in os.listdir(fileroot) if name.endswith('_link.cif'))
        except OSError:
            return None
        return os.path.join(fileroot, names[0]) if names else None

    def addModifications(self, parent=None):
        """How each monomer changes when the link forms, from the dictionary.

        The monomer's own block in the link dictionary is left as it was;
        the change is in a modification block that refinement applies to the
        linked residues. Laid out here, so that a bond made single by the
        link is seen to be single without reading the _chem_mod loops.
        """
        if parent is None: parent = self
        path = self.linkDictionaryPath()
        if path is None:
            return
        try:
            mods = link_modifications(path)
        except Exception as e:
            parent.addText(text='The link dictionary could not be read: %s' % e, tone='warning')
            return
        if not mods:
            return

        parent.addText(text='How each monomer changes when the link forms', tone='heading')
        parent.addText(
            text='AceDRG leaves each monomer\'s own description as it was and writes the '
                 'changes as a modification that refinement applies to the residues the '
                 'link joins. A viewer drawing the monomer alone shows the unlinked chemistry.',
            tone='muted')
        for mod in mods:
            parent.addText(text='%s (modification %s)' % (mod.comp_id, mod.mod_id), tone='emphasis')
            changes, atoms, details = [], [], []
            for bond in mod.bonds:
                changes.append('Bond')
                atoms.append('%s - %s' % (bond.atom_1, bond.atom_2))
                details.append(describe_bond(bond))
            for atom in mod.atoms:
                changes.append('Atom added' if atom.function == 'add' else 'Atom deleted')
                atoms.append(atom.name)
                what = atom.element or ''
                if atom.bonded_to:
                    what = ('%s on %s' % (what, atom.bonded_to)).strip()
                details.append(what)
            if changes:
                table = parent.addTable()
                table.addData(title='Change', data=changes)
                table.addData(title='Atoms', data=atoms)
                table.addData(title='Detail', data=details)
            also = ['%d %s' % (count, label if count != 1 else label[:-1])
                    for label, count in mod.revised.items()]
            if mod.retyped:
                also.append('the atom types of ' + ', '.join(mod.retyped))
            if also:
                parent.addText(text='Also revised: ' + '; '.join(also) + '.', tone='muted')

    def addDefaultReport(self, parent=None):
        if parent is None: parent=self
        for AcedrgLinkNode in self.xmlnode.findall("."):
            try:
                cycleNode = AcedrgLinkNode.findall("Cycle")[0]
            except:
                print("Missing cycle")
            try:
                logTextNode = AcedrgLinkNode.findall("LogText")[0]
            except:
                print("Missing logText")
            try:
                newFold = parent.addFold(label="Log text for iteration "+cycleNode.text, initiallyOpen=True)
            except:
                print("Unable to make fold")
            try:
                newFold.addPre(text = logTextNode.text)
            except:
                print("Unable to add Pre")

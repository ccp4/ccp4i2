import html

from ccp4i2.report import Report


def _text(xmlnode, tag):
    found = xmlnode.find('.//' + tag)
    if found is None or not found.text:
        return None
    return found.text.strip()


def safe(text):
    """Server-provided text, escaped: addText keeps inline markup."""
    return html.escape(text or '', quote=False)


def depui_url(site_url, dep_id):
    """The page a depositor opens: the site's /deposition root, not the API view."""
    if site_url and '/deposition' in site_url:
        return site_url.split('/deposition', 1)[0] + '/deposition'
    if dep_id:
        return 'https://deposit.wwpdb.org/deposition'
    return None


class onedep_submit_report(Report):
    TASKNAME = 'onedep_submit'
    RUNNING = True

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        dep_id = safe(_text(xmlnode, 'DEP_ID'))
        site_url = _text(xmlnode, 'SITE_URL')
        status = safe(_text(xmlnode, 'STATUS'))
        details = safe(_text(xmlnode, 'DETAILS'))
        error = safe(_text(xmlnode, 'ERROR'))
        replaced = _text(xmlnode, 'REPLACED') == '1'
        link = depui_url(site_url, dep_id)
        link = html.escape(link, quote=True) if link else None

        if jobStatus in ('Running', 'Running remotely'):
            if dep_id:
                self.addText(text=f'OneDep deposition {dep_id} created; waiting for OneDep '
                                  'to process the uploaded files.', tone='emphasis')
            else:
                self.addText(text='Creating the OneDep deposition and uploading the files.',
                             tone='emphasis')
            return

        if dep_id:
            self.addText(text=f'OneDep deposition {dep_id}', tone='heading')
            if link:
                self.addText(text=f'Open it at <a href="{link}" target="_blank">{link}</a> '
                                  'and sign in with the ORCID the key belongs to.', tone='info')
        if error:
            self.addText(text=error, tone='error')
            if replaced and dep_id:
                self.addText(text='A replacement removes the server\'s copy of a file before '
                                  'uploading the new one, so the deposition may now be '
                                  'without the file that failed until this task is re-run.',
                             tone='warning')
        elif status == 'finished':
            self.addText(text=f'OneDep processed the upload: {details or "done"}.',
                         tone='success')
        elif status:
            self.addText(text=f'OneDep reports {status}: {details or ""}', tone='warning')

        files = xmlnode.findall('.//FILES/FILE')
        if files:
            table = self.addTable(title='Files on the deposition')
            table.addData(title='File', data=[safe(_text(f, 'NAME')) for f in files])
            table.addData(title='Type', data=[safe(_text(f, 'TYPE')) for f in files])
            table.addData(title='OneDep file id', data=[safe(_text(f, 'REMOTE_ID')) for f in files])

        if dep_id and not error:
            fold = self.addFold(label='Next steps', brief='Next', initiallyOpen=True)
            if replaced:
                fold.addText(text='The files were replaced and reprocessed, but OneDep reads '
                                  'the descriptive metadata (title, authors, sources, data '
                                  'collection) from the model only when a deposition is first '
                                  'populated: the pages already filled in the interface are '
                                  'not updated from the new file. Edit them there, or create '
                                  'a new deposition if the metadata changed substantially.',
                             tone='warning')
            else:
                fold.addText(text='The title, authors, sources, sequence references and data '
                                  'collection details were read from the model file; check them '
                                  'on the Admin, Macromolecules and Data collection pages.',
                             tone='info')
            fold.addText(text='Still to enter in the OneDep interface: contact phone number, '
                              'grant information, release status, crystallisation condition '
                              'if it was not in the file, ligand descriptions, and the terms.',
                         tone='info')
            fold.addText(text='To replace the files after a change, re-run this task with the '
                              'deposition id in "Replace files on existing deposition".',
                         tone='muted')

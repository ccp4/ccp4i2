# ## @package arcimboldo_report
# # The initialisation method creates the whole report.
# # Offers functionalised report components for use on
# # pipeline reports.

import os
import time

from ccp4i2.core.CCP4ErrorHandling import SEVERITY_WARNING
from ccp4i2.report import Report


class arcimboldo_report ( Report ) :
	TASKNAME="arcimboldo"
	RUNNING = True
	CSS_VERSION = '0.1.0'
	WATCHED_FILE = 'program.xml'

	def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,taskVersion=None,**kw):
		Report.__init__(self,xmlnode=xmlnode,jobInfo=jobInfo,cssVersion=self.CSS_VERSION,**kw)

		if self.errorReport().maxSeverity()>SEVERITY_WARNING:
			print('FAILED instantiating MyWrapper report generator')
			self.errorReport().report()
			return

		# ARCIMBOLDO exits 0 even when it aborts, so the wrapper records any
		# FATAL lines it found in the log into PROGRAMXML. Show them first --
		# without this the report for a failed run is an empty "job finished".
		self.fatalErrorReport()

		fileXml = self.xmlnode.findall('pathXml')[0].text

		if os.path.exists(fileXml):
			self.arcimboldoXml = self.loadXmlFile(fileXml)

			div = self.addDiv(style='clear:both;')
			if jobStatus is not None and jobStatus.lower() == 'nooutput':
				return
			elif jobStatus is not None and jobStatus.lower() == 'running':
				self.defaultReport('running', parent=div)
			else:
				self.defaultReport('finished', parent=div)

			self.summaryReport()
			self.statusReport()
			self.borReport(jobInfo)

		elif jobStatus is not None and jobStatus.lower() != 'nooutput' and jobStatus.lower() != 'running':
			self.borReport(jobInfo)
		else:
			self.append('<p><b>Starting job. Waiting for results.</b></p>')

	def fatalErrorReport(self, parent=None):
		"""Render the FATAL messages ARCIMBOLDO wrote to its log, if any."""
		if parent is None:
			parent = self

		fatalErrors = [
			node.text for node in self.xmlnode.findall('fatalError') if node.text
		]
		if not fatalErrors:
			return

		div = parent.addDiv(style='clear:both;')
		html = ('<div style="border:1px solid #d32f2f; border-left:6px solid #d32f2f; '
				'background:#fdecea; padding:0.6em 1em; margin:0.5em 0;">'
				'<p style="color:#d32f2f; margin:0 0 0.4em 0;"><b>ARCIMBOLDO aborted</b></p>')
		for message in fatalErrors:
			html += '<p style="margin:0.2em 0;">{0}</p>'.format(
				self.escapeHtml(message))
		html += ('<p style="margin:0.6em 0 0 0;">See the log file for full details. '
				 'If a program could not be found, set its location in '
				 'Preferences &#8594; Program locations.</p></div>')
		div.append(html)

	@staticmethod
	def escapeHtml(text):
		return (str(text).replace('&', '&amp;')
						 .replace('<', '&lt;')
						 .replace('>', '&gt;'))

	def defaultReport(self, status, parent=None):
		if parent is None:
			parent = self

		resultsHtml = os.path.normpath(self.xmlnode.findall('pathHtml')[0].text)

		if status == 'running':
			htmlText = '<div><p><b>The job is currently running.</b></p>'
		else:
			htmlText = '<div><p><b>The job is finished.</b></p>'
			try:
				if len(self.arcimboldoXml.findall('backtracing/finalcc'))>0:
					finalcc = self.arcimboldoXml.findall('backtracing/finalcc')[0].text
					if float(finalcc) >= 30:
						htmlText += '<p style=\"color:green\"><b>It seems you have a good solution! (Final CC: {0}%)</b></p>'.format(finalcc)
			except:
				pass
		htmlText += '</div>'
		parent.append(htmlText)
		if os.path.exists(resultsHtml):
			parent.addFileLink(
				label='ARCIMBOLDO results (Last modified: {0})'.format(time.ctime(os.path.getmtime(resultsHtml))),
				relativePath=os.path.basename(resultsHtml),
				fileType='html',
			)

	def summaryReport(self, parent=None):
		if parent is None:
			parent = self

		summaryDiv = parent.addDiv(style='clear:both;')
		summaryFold = summaryDiv.addFold(label='Summary', brief='Summary', initiallyOpen=True)

		completeness = []
		try:
			completeness.append(self.arcimboldoXml.findall('data/completeness')[0].text)
		except:
			completeness.append('Not found')

		spaceGroup = []
		try:
			spaceGroup.append(self.arcimboldoXml.findall('data/spacegroup')[0].text)
		except:
			spaceGroup.append('Not found')

		cells = []
		try:
			cell = self.arcimboldoXml.findall('data/cell_dim/A')[0].text + ', ' + \
		 				self.arcimboldoXml.findall('data/cell_dim/B')[0].text + ', ' + \
		 				self.arcimboldoXml.findall('data/cell_dim/C')[0].text + ', ' + \
		 				self.arcimboldoXml.findall('data/cell_dim/alpha')[0].text + ', ' + \
		 				self.arcimboldoXml.findall('data/cell_dim/beta')[0].text +  ', ' + \
		 				self.arcimboldoXml.findall('data/cell_dim/gamma')[0].text
			cells.append(cell)
		except:
			cells.append('Not found')

		resolution = []
		try:
			resolution.append(self.arcimboldoXml.findall('data/resolution')[0].text)
		except:
			resolution.append('Not found')

		refl = []
		try:
			refl.append(self.arcimboldoXml.findall('data/unique_refl')[0].text)
		except:
			refl.append('Not found')

		table = summaryFold.addTable(transpose=True)
		table.addData(title = 'Completeness', data = completeness)
		table.addData(title = 'Space group', data = spaceGroup)
		table.addData(title = 'Unit cell', data = cells)
		table.addData(title = 'Resolution', data = resolution)
		table.addData(title = 'Number of unique reflections', data = refl)

	def statusReport(self, parent=None):
		if parent is None:
			parent = self

		statusDiv = parent.addDiv(style='clear:both;')
		statusFold = statusDiv.addFold(label='Status', brief='Status', initiallyOpen=True)

		try:
			if len(self.arcimboldoXml.findall('TIME'))>0:
				logInfo = self.arcimboldoXml.findall('TIME')[0].text.replace('[','').replace(']]','').split('],')

				timing = []
				step = []
				for x in logInfo:
					aux = x.split(',')
					step.append(aux[1])
					timing.append(aux[2])

				table = statusFold.addTable()
				table.addData(title = 'Step', data = step)
				table.addData(title = 'Time', data = timing)

			else:
				raise Exception
		except:
			statusFold.append('Not found')

	def borReport(self, jobInfo, parent=None):
		if parent is None:
			parent = self

		borDiv = parent.addDiv(style='clear:both;')
		borFold = borDiv.addFold(label='Setup file', brief='Setup', initiallyOpen=False)

		borFile = os.path.join(jobInfo['fileroot'],'setup.bor')
		try:
			if os.path.exists(borFile):
				with open(borFile, 'r') as bor:
					for line in bor:
						borFold.append(line)
			else:
				raise Exception
		except:
			borFold.append('Not found')

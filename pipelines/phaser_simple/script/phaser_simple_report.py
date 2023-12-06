from report.CCP4ReportParser import *
import sys
from lxml import etree
from pipelines.phaser_pipeline.script.phaser_pipeline_report import phaser_pipeline_report

class phaser_simple_report(phaser_pipeline_report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'phaser_simple'
    



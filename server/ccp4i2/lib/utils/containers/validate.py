import logging
from xml.etree import ElementTree as ET
from ccp4i2.core import CCP4Container
from ccp4i2.core import CCP4ErrorHandling


logger = logging.getLogger(f"ccp4x:{__name__}")

# Mapping from severity codes to text
SEVERITY_TEXT = {
    CCP4ErrorHandling.SEVERITY_OK: "OK",
    CCP4ErrorHandling.SEVERITY_UNDEFINED: "UNDEFINED",
    CCP4ErrorHandling.SEVERITY_WARNING: "WARNING",
    CCP4ErrorHandling.SEVERITY_UNDEFINED_ERROR: "UNDEFINED_ERROR",
    CCP4ErrorHandling.SEVERITY_ERROR: "ERROR"
}


def getEtree(error_report: CCP4ErrorHandling.CErrorReport):
    element = ET.Element("errorReportList")
    # Use getErrors() public API instead of accessing _reports directly
    for item in error_report.getErrors():
        try:
            ele = ET.Element("errorReport")
            e = ET.Element("className")
            # In new API, 'class' is a string, not a class object
            class_name = item["class"] if isinstance(item["class"], str) else item["class"].__name__
            e.text = class_name
            ele.append(e)
            e = ET.Element("code")
            e.text = str(item["code"])
            ele.append(e)
            e = ET.Element("description")
            # Description is in 'details' field in new API
            e.text = item["details"]
            ele.append(e)
            e = ET.Element("severity")
            severity = item["severity"]
            e.text = SEVERITY_TEXT.get(severity, f"UNKNOWN({severity})")
            ele.append(e)
            if item["details"] is not None:
                e = ET.Element("details")
                e.text = str(item["details"])
                ele.append(e)
            if item.get("time", None) is not None:
                e = ET.Element("time")
                e.text = str(item["time"])
                ele.append(e)
            if item.get("stack", None) is not None:
                e = ET.Element("stack")
                # print 'CErrorReport.getEtree stack',item['stack'],type(item['stack']),type(item['stack'][0])
                text = ""
                for line in item["stack"]:
                    text = text + line
                e.text = text
                ele.append(e)
            element.append(ele)
        except Exception as e:
            logger.exception("Error in validate_container", exc_info=e)
    return element


def validate_container(
    container: CCP4Container.CContainer,
) -> ET.Element:
    error_report: CCP4ErrorHandling.CErrorReport = container.validity()
    error_etree: ET.Element = getEtree(error_report)

    # Remove stacks (too nasty to live with)
    namespace = ""
    error_reports = error_etree.findall(".//{0}errorReport".format(namespace))
    for error_report in error_reports:
        stack_children = error_report.findall("./stack")
        for stack_child in stack_children:
            error_report.remove(stack_child)
        description_children = error_report.findall("./description")

        for description_child in description_children:
            description_text = description_child.text
            broken_text = description_text.split(":")
            if len(broken_text) > 1:
                object_element = ET.Element("objectPath")
                object_element.text = broken_text[0].strip()
                error_report.append(object_element)

    ET.indent(error_etree, " ")
    # print(error_log)
    return error_etree

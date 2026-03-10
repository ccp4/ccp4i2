import os
import re
from collections import OrderedDict

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core import CCP4Utils


class ProvideTLS(CPluginScript):

    TASKNAME = 'ProvideTLS'

    def process(self):
        invalidFiles = self.checkInputData()
        if len(invalidFiles) > 0:
            self.reportStatus(CPluginScript.FAILED)

        self.checkOutputData()

        # Use guiParameters.EDIT_MODE to decide which source to use
        edit_mode = 'table'
        if hasattr(self.container, 'guiParameters') and \
                hasattr(self.container.guiParameters, 'EDIT_MODE') and \
                self.container.guiParameters.EDIT_MODE.isSet():
            edit_mode = str(self.container.guiParameters.EDIT_MODE)

        if edit_mode == 'table':
            tls_text = self._generate_tls_text_from_groups()
            # Fall back chain: TLSGROUPS → TLSIN file → TLSTEXT
            if not tls_text.strip():
                tls_text = self._read_tlsin_file()
            if not tls_text.strip():
                tls_text = self.container.controlParameters.TLSTEXT.__str__()
        else:
            tls_text = self.container.controlParameters.TLSTEXT.__str__()

        with open(self.container.outputData.TLSFILE.fullPath.__str__(), "w") as myFile:
            myFile.write(tls_text)

        from lxml import etree
        root = etree.Element('ProvideTLSOutput')
        tlsElement = etree.SubElement(root, 'TLSProvided')
        tlsElement.text = tls_text
        with open(self.makeFileName('PROGRAMXML'), 'w') as xmlFile:
            CCP4Utils.writeXML(xmlFile, etree.tostring(root, pretty_print=True))

        self.reportStatus(CPluginScript.SUCCEEDED)

    def _read_tlsin_file(self):
        """Read TLSIN input file if provided, return contents or empty string."""
        tlsin = self.container.inputData.TLSIN
        if tlsin is None or not tlsin.isSet():
            return ""
        file_path = str(tlsin.fullPath)
        if not file_path or not os.path.exists(file_path):
            return ""
        with open(file_path, 'r') as f:
            return f.read()

    def _generate_tls_text_from_groups(self):
        """Generate TLS file text from structured TLSGROUPS parameter."""
        groups_list = self.container.controlParameters.TLSGROUPS
        if groups_list is None:
            return ""

        items = list(groups_list)
        if not items:
            return ""

        # Group ranges by groupId
        groups = OrderedDict()
        for item in items:
            data = item.get()
            gid = int(data.get('groupId', 0) or 0)
            if gid not in groups:
                groups[gid] = []
            groups[gid].append(data)

        lines = []
        for gid, ranges in groups.items():
            lines.append("TLS    Group %d" % gid)
            for r in ranges:
                chain = str(r.get('chainId', 'A') or 'A')
                first = str(r.get('firstRes', '0') or '0')
                last = str(r.get('lastRes', '999') or '999')
                sel = str(r.get('selection', 'ALL') or 'ALL')
                lines.append(
                    "RANGE  '%s%4s.' '%s%4s.' %s" % (chain, first, chain, last, sel)
                )
            lines.append("")  # blank line between groups

        return "\n".join(lines)

    def suggest_tls_groups(self):
        """
        Suggest TLS groups from the XYZIN coordinate file.

        Creates one TLS group per polymer chain, covering the full residue range.
        Called via the plugin_method API endpoint from the frontend.

        Returns:
            list of dicts: [{groupId, chainId, firstRes, lastRes, selection}, ...]
        """
        import gemmi
        from ccp4i2.core.CCP4ModelData import CPdbDataComposition

        xyzin = self.container.inputData.XYZIN
        if xyzin is None or not xyzin.isSet():
            return {"error": "No coordinate file (XYZIN) provided"}

        file_path = str(xyzin.fullPath)
        if not file_path or not os.path.exists(file_path):
            return {"error": "Coordinate file not found: %s" % file_path}

        try:
            structure = gemmi.read_structure(file_path)
            comp = CPdbDataComposition(structure)
        except Exception as e:
            return {"error": "Failed to read coordinate file: %s" % str(e)}

        suggested = []
        group_id = 1

        for detail in comp.chainDetails:
            # Only polymer chains (protein or nucleic)
            if detail["type"] not in ("protein", "nucleic"):
                continue

            first_res = detail["firstRes"]
            last_res = detail["lastRes"]
            try:
                first_res = int(first_res)
            except (ValueError, TypeError):
                first_res = 1
            try:
                last_res = int(last_res)
            except (ValueError, TypeError):
                last_res = 999

            suggested.append({
                "groupId": group_id,
                "chainId": detail["id"],
                "firstRes": first_res,
                "lastRes": last_res,
                "selection": "ALL",
            })
            group_id += 1

        return suggested

    @staticmethod
    def parse_tls_text(tls_text):
        """
        Parse TLS-format text into structured group data.

        Handles the standard TLS file format:
            TLS    <group name>
            RANGE  '<chain><resnum>.' '<chain><resnum>.' <selection>

        Returns:
            list of dicts: [{groupId, chainId, firstRes, lastRes, selection}, ...]
        """
        result = []
        current_group_id = 0

        for line in tls_text.splitlines():
            stripped = line.strip()
            if not stripped:
                continue

            if stripped.upper().startswith('TLS'):
                current_group_id += 1
                continue

            range_match = re.match(
                r"RANGE\s+'(\w)\s*([^']+)'\s+'(\w)\s*([^']+)'\s+(\w+)",
                stripped, re.IGNORECASE
            )
            if range_match:
                chain = range_match.group(1)
                first_res = range_match.group(2).strip().rstrip('.')
                last_res = range_match.group(4).strip().rstrip('.')

                try:
                    first_res_int = int(first_res)
                except (ValueError, TypeError):
                    first_res_int = 0
                try:
                    last_res_int = int(last_res)
                except (ValueError, TypeError):
                    last_res_int = 999

                result.append({
                    "groupId": current_group_id if current_group_id > 0 else 1,
                    "chainId": chain,
                    "firstRes": first_res_int,
                    "lastRes": last_res_int,
                    "selection": range_match.group(5).upper(),
                })

        return result

    def import_tls_from_text(self):
        """
        Parse TLSTEXT into structured groups. Called via plugin_method.

        Returns:
            list of dicts: [{groupId, chainId, firstRes, lastRes, selection}, ...]
        """
        tls_text = str(self.container.controlParameters.TLSTEXT)
        return self.parse_tls_text(tls_text)

    def import_tls_from_file(self):
        """
        Parse TLSIN file into structured groups. Called via plugin_method.

        Returns:
            list of dicts or error dict
        """
        tlsin = self.container.inputData.TLSIN
        if tlsin is None or not tlsin.isSet():
            return {"error": "No TLS file (TLSIN) provided"}

        file_path = str(tlsin.fullPath)
        if not file_path or not os.path.exists(file_path):
            return {"error": "TLS file not found: %s" % file_path}

        with open(file_path, 'r') as f:
            tls_text = f.read()

        return self.parse_tls_text(tls_text)

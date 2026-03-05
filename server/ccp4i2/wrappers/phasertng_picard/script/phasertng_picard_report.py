"""
phasertng_picard report — builds a CCP4i2 report from PhaserTNG Picard output.

Picard writes structured text files rather than JSON:
  - result.cards: key-value metadata (resolution, spacegroup, R-factor, etc.)
  - best.N.dag.cards: solution blocks separated by '===='
  - *.dag.html: per-node pathway files with parent/child identifiers and Z-scores

This report follows the ModelCraft pattern (USEPROGRAMXML=False, RUNNING=True)
and adds a vis-network solution pathway tree matching Airlie McCoy's pyvis output.
"""

import json
import logging
import re
from pathlib import Path

from ccp4i2.report.CCP4ReportParser import Report

logger = logging.getLogger(__name__)


def _dag_node_colour(tag, rfactor=0.0):
    """Return a colour string matching runTree.py's tag-prefix colour scheme."""
    t = tag.lower()
    if t.startswith(("root", "find")):
        return "green"
    elif t.startswith("join"):
        return "gold"
    elif t.startswith("put"):
        return "grey"
    elif t.startswith("fuse"):
        return "darkkhaki"
    elif t.startswith(("rfac", "xref")):
        if rfactor is not None and rfactor > 0:
            if rfactor > 45:
                return "#ff0000"
            if rfactor > 40:
                return "#d5002a"
            if rfactor > 35:
                return "#aa0055"
            if rfactor > 31:
                return "#800080"
            if rfactor > 27:
                return "#5500aa"
            if rfactor > 24:
                return "#2b00d5"
            return "#0000ff"
        return "skyblue"
    return "orange"


def _dag_node_shape(tag):
    """Return a vis-network shape matching runTree.py's tag-prefix mapping."""
    t = tag.lower()
    if t.startswith("find") or t.startswith("fuse"):
        return "square"
    elif t.startswith("join"):
        return "triangle"
    return "dot"


class phasertng_picard_report(Report):
    TASKNAME = "phasertng_picard"
    RUNNING = True
    USEPROGRAMXML = False
    WATCHED_FILE = "phasertng_picard/best.1.dag.cards"

    @staticmethod
    def _find_db_dir(fileroot):
        """Find the phasertng database subdirectory (picard or riker)."""
        for name in ("phasertng_picard", "phasertng_riker"):
            candidate = Path(fileroot, name)
            if candidate.is_dir():
                return candidate
        return None

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(
            self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw
        )
        self.addDiv(style="clear:both;")
        if jobStatus in ["Running", "Running remotely"]:
            self.append("<p><b>The job is currently running.</b></p>")

        db_dir = self._find_db_dir(jobInfo.get("fileroot", ""))
        if db_dir is None:
            self.append("<p>No output directory found.</p>")
            return

        # Parse output files
        result = self._parse_result_cards(db_dir)
        solutions = self._parse_dag_cards(db_dir)
        dag_dot = self._parse_dag_html_tree(db_dir)

        # Build report sections
        if result or solutions:
            summaryFold = self.addFold(
                label="Summary", initiallyOpen=True, brief="Summary"
            )
            self._add_summary_table(result, solutions, parent=summaryFold)
            self._add_quality_message(result, solutions, jobStatus, parent=summaryFold)

        if len(solutions) > 1:
            solFold = self.addFold(
                label="Solutions", initiallyOpen=True, brief="Solutions"
            )
            self._add_solutions_table(solutions, parent=solFold)

        if dag_dot:
            dagFold = self.addFold(
                label="Solution Pathway", initiallyOpen=True, brief="Pathway"
            )
            dagFold.addDAGGraph(
                title="PhaserTNG Solution Pathway",
                elements=dag_dot,
                layout="hierarchical",
            )

        if result:
            self._add_notifications(result, parent=self)

        self.addDiv(style="clear:both;")

    # ------------------------------------------------------------------
    # Parsers
    # ------------------------------------------------------------------

    def _parse_result_cards(self, db_dir):
        """Parse result.cards from the highest-numbered rfac subdirectory.

        Returns a dict where keys are space-separated paths after 'phasertng '.
        Multi-value keys (e.g., notifications advisory) become lists.
        """
        rfac_dirs = sorted(db_dir.glob("*-rfac"))
        if not rfac_dirs:
            return {}

        cards_files = list(rfac_dirs[-1].glob("*.result.cards"))
        if not cards_files:
            return {}

        result = {}
        try:
            text = cards_files[0].read_text(encoding="utf-8", errors="replace")
        except OSError:
            return {}

        for line in text.splitlines():
            line = line.strip()
            if not line or not line.startswith("phasertng "):
                continue
            # Format: "phasertng <key>   <value>"
            # Split on 2+ spaces to separate key from value
            rest = line[len("phasertng "):]
            parts = re.split(r"  +", rest, maxsplit=1)
            if len(parts) == 2:
                key, value = parts[0].strip(), parts[1].strip()
            else:
                key, value = parts[0].strip(), ""

            # Strip surrounding quotes
            if value.startswith('"') and value.endswith('"'):
                value = value[1:-1]

            # Multi-value keys become lists
            if key in result:
                existing = result[key]
                if isinstance(existing, list):
                    existing.append(value)
                else:
                    result[key] = [existing, value]
            else:
                result[key] = value

        return result

    def _parse_dag_cards(self, db_dir):
        """Parse best.N.dag.cards into a list of solution dicts.

        Each solution block is separated by '====' lines. The first block
        is the best solution.
        """
        dag_files = sorted(db_dir.glob("best.*.dag.cards"))
        if not dag_files:
            return []

        try:
            text = dag_files[0].read_text(encoding="utf-8", errors="replace")
        except OSError:
            return []

        blocks = re.split(r"^={4,}\s*$", text, flags=re.MULTILINE)
        solutions = []

        for block in blocks:
            block = block.strip()
            if not block:
                continue

            sol = {}
            for line in block.splitlines():
                line = line.strip()
                if not line.startswith("phaserdag node "):
                    continue
                rest = line[len("phaserdag node "):]
                # Split on first space after the key
                parts = rest.split(None, 1)
                if len(parts) == 2:
                    key, value = parts
                else:
                    continue

                # Strip quotes
                if value.startswith('"') and value.endswith('"'):
                    value = value[1:-1]

                # Store first occurrence only for simple keys
                if key not in sol:
                    sol[key] = value

            if sol:
                solutions.append(sol)

        return solutions

    def _parse_dag_html_tree(self, db_dir):
        """Walk db_dir for *.dag.html files and build a solution pathway tree.

        Matches runTree.py / inner_voyager_pyvis: each dag.html has a 5-line
        HTML-comment header encoding parent→child pathway relationships,
        Z-scores, and R-factors.

        Returns a JSON string with {nodes: [...], edges: [...]},
        or empty string if no dag.html files found.
        """
        dag_files = list(db_dir.rglob("*.dag.html"))
        if not dag_files:
            return ""

        nodes_dict = {}   # identifier (str or int 0) -> node dict
        edges_set = set()  # (from_id, to_id) tuples for dedup

        def _make_label(tag, str_id):
            """Build label like 'rfac-28' from tag and zero-padded id string."""
            short = str(int(str_id))[:3]
            return "{}-{}".format(tag, short)

        for dag_file in dag_files:
            try:
                with open(dag_file, "r", encoding="utf-8", errors="replace") as f:
                    lines = [f.readline() for _ in range(5)]
            except OSError:
                continue

            if len(lines) < 5 or not lines[0].strip():
                continue

            try:
                # Line 0: <!-- identifier tag -->
                w0 = lines[0].split()
                uid_id = w0[1]          # keep as string e.g. "0000000073"
                uid_tag = w0[2]

                # Line 1: <!-- identifier tag --> or <!-- --- -->
                w1 = lines[1].split()
                if w1[1] == "---":
                    pid_id = 0           # root is int 0
                    pid_tag = "root"
                else:
                    pid_id = w1[1]       # keep as string
                    pid_tag = w1[2]

                # Line 4: <!-- tracker_id zscore rfactor -->
                w4 = lines[4].split()
                zscore = float(w4[1])
                rfactor = float(w4[2])
            except (IndexError, ValueError):
                continue

            # Ensure root node exists
            if 0 not in nodes_dict:
                nodes_dict[0] = {
                    "id": 0,
                    "label": db_dir.name or "root",
                    "color": "green",
                    "shape": "dot",
                    "value": 20,
                }

            # Size from Z-score (matching runTree.py)
            size = 20 * max(1, zscore) if zscore > 0 else 20

            # Add parent node if not already present
            if pid_id != 0 and pid_id not in nodes_dict:
                nodes_dict[pid_id] = {
                    "id": pid_id,
                    "label": _make_label(pid_tag, pid_id),
                    "color": _dag_node_colour(pid_tag),
                    "shape": _dag_node_shape(pid_tag),
                    "value": 20,
                }

            # Add/update current node (later files may refine the same node)
            nodes_dict[uid_id] = {
                "id": uid_id,
                "label": _make_label(uid_tag, uid_id),
                "color": _dag_node_colour(uid_tag, rfactor),
                "shape": _dag_node_shape(uid_tag),
                "value": size,
            }

            edges_set.add((pid_id, uid_id))

        if not nodes_dict:
            return ""

        nodes = list(nodes_dict.values())
        edges = [{"from": f, "to": t} for f, t in edges_set]
        return json.dumps({"nodes": nodes, "edges": edges})

    # ------------------------------------------------------------------
    # Report sections
    # ------------------------------------------------------------------

    def _add_summary_table(self, result, solutions, parent):
        """Add a transposed key-value summary table."""
        table = parent.addTable(transpose=True)

        best = solutions[0] if solutions else {}

        # Space group
        sg = best.get("hermann_mauguin", result.get("subgroup original hermann_mauguin", ""))
        if sg:
            table.addData(title="Space group", data=[sg])

        # Unit cell
        cell = best.get("unitcell", "")
        if cell:
            parts = cell.split()
            if len(parts) >= 6:
                cell_str = (
                    f"{float(parts[0]):.2f}  {float(parts[1]):.2f}  {float(parts[2]):.2f}  "
                    f"{float(parts[3]):.1f}  {float(parts[4]):.1f}  {float(parts[5]):.1f}"
                )
                table.addData(title="Unit cell", data=[cell_str])

        # Resolution
        res = result.get("data resolution_available", "")
        if res:
            parts = res.split()
            if parts:
                table.addData(title="Resolution (\u00C5)", data=[f"{float(parts[0]):.2f}"])

        # Wilson B-factor
        wilson_b = result.get("anisotropy wilson_bfactor", "")
        if wilson_b:
            table.addData(
                title="Wilson B-factor (\u00C5\u00B2)", data=[f"{float(wilson_b):.1f}"]
            )

        # R-factor
        rfac = best.get("rfactor", "")
        if rfac:
            try:
                rfac_val = float(rfac)
                if rfac_val <= 100:
                    table.addData(title="R-factor (%)", data=[f"{rfac_val:.1f}"])
            except ValueError:
                pass

        # TFZ score
        tfz = best.get("zscore", "")
        if tfz:
            table.addData(title="TFZ score", data=[f"{float(tfz):.1f}"])

        # LLG
        llg = best.get("llg", "")
        if llg:
            table.addData(title="LLG", data=[f"{float(llg):.1f}"])

        # Matthews coefficient
        matthews_lines = result.get("matthews", "")
        if isinstance(matthews_lines, str) and matthews_lines:
            # Parse "vm  2.10 probability  0.97 z  1"
            m = re.search(r"vm\s+([\d.]+)\s+probability\s+([\d.]+)\s+z\s+(\d+)", matthews_lines)
            if m:
                table.addData(title="V<sub>M</sub> (\u00C5\u00B3/Da)", data=[f"{float(m.group(1)):.2f}"])
                table.addData(title="Z (copies in ASU)", data=[m.group(3)])

        # Twinning
        twinned = best.get("twinned", result.get("twinning indicated", ""))
        if twinned:
            table.addData(title="Twinning", data=["Yes" if twinned.lower() == "true" else "No"])

        # tNCS
        tncs = best.get("tncs_indicated", "")
        if tncs:
            table.addData(title="tNCS", data=["Yes" if tncs.lower() == "true" else "No"])

        # Number of solutions
        if solutions:
            table.addData(title="Solutions found", data=[str(len(solutions))])

        # Runtime
        wall = result.get("time cumulative wall", "")
        if wall:
            try:
                secs = float(wall)
                if secs < 60:
                    table.addData(title="Runtime", data=[f"{secs:.0f} s"])
                elif secs < 3600:
                    table.addData(title="Runtime", data=[f"{secs / 60:.1f} min"])
                else:
                    table.addData(title="Runtime", data=[f"{secs / 3600:.1f} h"])
            except ValueError:
                pass

    def _add_quality_message(self, result, solutions, jobStatus, parent):
        """Add a quality assessment message based on R-factor and TFZ."""
        best = solutions[0] if solutions else {}
        rfac_str = best.get("rfactor", "")
        tfz_str = best.get("zscore", "")

        if not rfac_str or not tfz_str:
            return

        try:
            rfac = float(rfac_str)
            tfz = float(tfz_str)
        except ValueError:
            return

        if rfac > 60:
            # R-factor > 60% means the solution is likely wrong or incomplete
            return

        if jobStatus in ["Running", "Running remotely"]:
            message = "The best solution found so far"
        else:
            message = "The best solution"

        if rfac < 30 and tfz > 8:
            message += " is a strong molecular replacement solution."
        elif rfac < 40 and tfz > 6:
            message += " appears promising but may require further refinement."
        elif tfz > 8:
            message += " has a strong TFZ signal but a high R-factor, suggesting the model may need significant rebuilding."
        else:
            message += " has weak statistics and may not be correct."

        parent.append(f"<p>{message}</p>")

    def _add_solutions_table(self, solutions, parent):
        """Add a table comparing all solutions."""
        table = parent.addTable()

        nums = []
        sgs = []
        rfacs = []
        tfzs = []
        llgs = []
        twins = []
        packs_list = []

        for i, sol in enumerate(solutions):
            nums.append(i + 1)
            sgs.append(sol.get("hermann_mauguin", "?"))

            rfac = sol.get("rfactor", "")
            try:
                rfac_val = float(rfac)
                rfacs.append(f"{rfac_val:.1f}" if rfac_val <= 100 else "-")
            except (ValueError, TypeError):
                rfacs.append("-")

            try:
                tfzs.append(f"{float(sol.get('zscore', '')):.1f}")
            except (ValueError, TypeError):
                tfzs.append("-")

            try:
                llgs.append(f"{float(sol.get('llg', '')):.0f}")
            except (ValueError, TypeError):
                llgs.append("-")

            tw = sol.get("twinned", "")
            twins.append("Yes" if tw.lower() == "true" else "No" if tw else "-")

            packs_list.append(sol.get("packs", "-"))

        table.addData(title="#", data=nums)
        table.addData(title="Space group", data=sgs)
        table.addData(title="R-factor (%)", data=rfacs)
        table.addData(title="TFZ", data=tfzs)
        table.addData(title="LLG", data=llgs)
        table.addData(title="Twinned", data=twins)
        table.addData(title="Packs", data=packs_list)

    def _add_notifications(self, result, parent):
        """Add warnings and advisories from result.cards."""
        warnings = result.get("notifications warning", [])
        advisories = result.get("notifications advisory", [])

        if isinstance(warnings, str):
            warnings = [warnings]
        if isinstance(advisories, str):
            advisories = [advisories]

        if not warnings and not advisories:
            return

        fold = parent.addFold(
            label="Notifications", initiallyOpen=False, brief="Notifications"
        )

        if warnings:
            items = "".join(
                f"<li>{w.replace('_', ' ')}</li>" for w in warnings
            )
            fold.append(f"<div><p><b>Warnings:</b></p><ul>{items}</ul></div>")

        if advisories:
            items = "".join(
                f"<li>{a.replace('_', ' ')}</li>" for a in advisories
            )
            fold.append(f"<div><p><b>Advisories:</b></p><ul>{items}</ul></div>")


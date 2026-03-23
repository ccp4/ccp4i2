import json
import logging

from ccp4i2.report import Report

logger = logging.getLogger(f"ccp4i2:{__name__}")


def _parse_newick_to_vis(newick_str):
    """Parse a Newick-format string into vis-network nodes and edges.

    Returns a dict with 'nodes' and 'edges' lists suitable for
    the CCP4i2ReportDAGGraph component.
    """
    nodes = []
    edges = []
    node_id = [0]

    def _next_id():
        node_id[0] += 1
        return node_id[0]

    def _shorten_label(label):
        """Extract a short name from e.g. 'sp|Q01094|E2F1_HUMAN'."""
        parts = label.split('|')
        if len(parts) >= 3:
            return parts[2].replace('_HUMAN', '')
        return label

    def _parse(s, pos):
        """Recursive descent parser for Newick format.

        Returns (node_id, new_pos).
        """
        children = []
        if pos < len(s) and s[pos] == '(':
            # Internal node with children
            pos += 1  # skip '('
            # skip whitespace
            while pos < len(s) and s[pos] in ' \t\n\r':
                pos += 1
            child_id, pos = _parse(s, pos)
            children.append(child_id)
            while pos < len(s):
                # skip whitespace before comma or closing paren
                while pos < len(s) and s[pos] in ' \t\n\r':
                    pos += 1
                if pos >= len(s) or s[pos] != ',':
                    break
                pos += 1  # skip ','
                while pos < len(s) and s[pos] in ' \t\n\r':
                    pos += 1
                child_id, pos = _parse(s, pos)
                children.append(child_id)
            if pos < len(s) and s[pos] == ')':
                pos += 1  # skip ')'

        # Read label (may be empty for internal nodes)
        label = ''
        while pos < len(s) and s[pos] not in ':,);(\n\r':
            label += s[pos]
            pos += 1
        label = label.strip()

        # Skip whitespace before optional branch length
        while pos < len(s) and s[pos] in ' \t\n\r':
            pos += 1

        # Read branch length
        branch_length = None
        if pos < len(s) and s[pos] == ':':
            pos += 1
            bl_str = ''
            while pos < len(s) and s[pos] not in ',);(\n\r':
                bl_str += s[pos]
                pos += 1
            try:
                branch_length = float(bl_str.strip())
            except ValueError:
                pass

        nid = _next_id()
        if children:
            # Internal node
            nodes.append({
                'id': nid,
                'label': '',
                'shape': 'dot',
                'color': 'grey',
                'value': 5,
                '_bl': branch_length,
            })
            for child_id in children:
                child_bl = None
                for n in nodes:
                    if n['id'] == child_id:
                        child_bl = n.get('_bl')
                        break
                edge = {'from': nid, 'to': child_id}
                if child_bl is not None:
                    edge['label'] = f'{child_bl:.4f}'
                edges.append(edge)
        else:
            # Leaf node
            short_label = _shorten_label(label) if label else f'seq_{nid}'
            nodes.append({
                'id': nid,
                'label': short_label,
                'shape': 'box',
                'color': '#4fc3f7',
                'value': 15,
                '_bl': branch_length,
            })

        return nid, pos

    # Clean up the newick string
    clean = newick_str.strip()
    if clean.endswith(';'):
        clean = clean[:-1]

    if not clean:
        return {'nodes': [], 'edges': []}

    _parse(clean, 0)

    # Remove internal _bl keys (not needed by vis-network)
    for n in nodes:
        n.pop('_bl', None)

    return {'nodes': nodes, 'edges': edges}


def _extract_sequence_names(alignment_text):
    """Extract unique sequence names from CLUSTAL alignment text."""
    seq_names = []
    seen = set()
    for line in alignment_text.split('\n'):
        line = line.strip()
        if not line:
            continue
        if line.startswith('CLUSTAL'):
            continue
        # Skip conservation lines (only contain *, :, ., spaces)
        if all(c in '.:* ' for c in line):
            continue
        parts = line.split()
        if len(parts) >= 2:
            name = parts[0]
            if name not in seen:
                seen.add(name)
                seq_names.append(name)
    return seq_names


def _shorten_name(name):
    """Shorten a UniProt-style name for table display."""
    parts = name.split('|')
    if len(parts) >= 3:
        return parts[2].replace('_HUMAN', '')
    return name


class clustalw_report(Report):
    TASKNAME = 'clustalw'
    RUNNING = False

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo,
                        jobStatus=jobStatus, **kw)

        if jobStatus is not None and jobStatus.lower() == 'nooutput':
            return

        self.defaultReport(parent=self)

    def defaultReport(self, parent=None):
        self.addDiv(style='clear:both;')
        if parent is None:
            parent = self

        # --- Dendrogram ---
        try:
            dendogram_node = self.xmlnode.find('Dendogram')
            if dendogram_node is not None and dendogram_node.text and dendogram_node.text.strip():
                vis_data = _parse_newick_to_vis(dendogram_node.text)
                if vis_data['nodes']:
                    parent.addDAGGraph(
                        title='Guide Tree (Dendrogram)',
                        elements=json.dumps(vis_data),
                    )
        except Exception:
            logger.warning("Failed to parse dendrogram data", exc_info=True)

        # --- Pairwise Identity Matrix ---
        try:
            pairwise_nodes = self.xmlnode.findall('PairwiseScore')
            if pairwise_nodes:
                self._add_identity_matrix(parent, pairwise_nodes)
        except Exception:
            logger.warning("Failed to build identity matrix", exc_info=True)

        # --- Alignment ---
        for alignmentNode in self.xmlnode.findall('Alignment'):
            fold = parent.addFold(label='Sequence Alignment', initiallyOpen=True)
            fold.addAlignment(text=alignmentNode.text)

        # --- Statistics ---
        for node in self.xmlnode.findall('Statistics'):
            fold = parent.addFold(label='Alignment Statistics', initiallyOpen=False)
            fold.addPre(text=node.text)

    def _add_identity_matrix(self, parent, pairwise_nodes):
        """Build a pairwise identity score matrix as a report Table.

        PairwiseScore partner indices are 1-based and refer to the input
        sequence order (as listed in the log). We prefer the SequenceList
        element (added by the wrapper) which preserves this order. If
        unavailable, we fall back to extracting names from the alignment
        text (which may be reordered by ClustalW).
        """
        # Try SequenceList first (input order, matches PairwiseScore indices)
        seq_list_node = self.xmlnode.find('SequenceList')
        if seq_list_node is not None:
            seq_names = [s.text for s in seq_list_node.findall('Sequence')
                         if s.text]
        else:
            # Fallback: extract from alignment text (may differ in order)
            alignment_node = self.xmlnode.find('Alignment')
            if alignment_node is None or not alignment_node.text:
                return
            seq_names = _extract_sequence_names(alignment_node.text)

        if not seq_names:
            return

        n = len(seq_names)

        # Build the score lookup from PairwiseScore elements
        scores = {}
        for pw in pairwise_nodes:
            score_el = pw.find('Score')
            partners = pw.findall('Partner')
            if score_el is None or len(partners) != 2:
                continue
            try:
                i = int(partners[0].text) - 1
                j = int(partners[1].text) - 1
                score = int(score_el.text)
            except (ValueError, TypeError):
                continue
            scores[(i, j)] = score
            scores[(j, i)] = score

        short_names = [_shorten_name(s) for s in seq_names]

        table = parent.addTable(title='Pairwise Identity Scores (%)')
        table.addData(title='Sequence', data=short_names)
        for j in range(n):
            col_data = []
            for i in range(n):
                if i == j:
                    col_data.append('-')
                else:
                    val = scores.get((i, j))
                    col_data.append(str(val) if val is not None else '-')
            table.addData(title=short_names[j], data=col_data)

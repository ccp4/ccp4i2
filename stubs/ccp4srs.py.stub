"""
Minimal stub for ccp4srs (CCP4 Structure Refinement Suite) module.

This provides minimal stubs for the ccp4srs functionality used by acedrgNew's atomMatching.
Most functionality will raise NotImplementedError if actually called.
"""

# Constants
EXTTYPE_Ignore = 0


class Manager:
    """CCP4 SRS Manager class for structure database access (stub)."""

    def loadIndex(self, path):
        """Load SRS index from path (stub - does nothing)."""
        pass

    def loadStructure(self, name):
        """Load structure by name (stub - does nothing)."""
        pass


class Graph:
    """Graph representation for molecular structure matching (stub)."""

    def MakeGraph(self, residue):
        """Create graph from residue (stub - not implemented)."""
        raise NotImplementedError("ccp4srs.Graph.MakeGraph is not implemented in stub")

    def RemoveChirality(self):
        """Remove chirality information (stub - not implemented)."""
        raise NotImplementedError("ccp4srs.Graph.RemoveChirality is not implemented in stub")

    def Build(self, flag):
        """Build graph structure (stub - not implemented)."""
        raise NotImplementedError("ccp4srs.Graph.Build is not implemented in stub")

    def GetNofVertices(self):
        """Get number of vertices (stub - not implemented)."""
        raise NotImplementedError("ccp4srs.Graph.GetNofVertices is not implemented in stub")

    def Print(self):
        """Print graph (stub - does nothing)."""
        pass


class GraphMatch:
    """Graph matching class for substructure searches (stub)."""

    def SetTimeLimit(self, seconds):
        """Set time limit for matching (stub - does nothing)."""
        pass

    def MatchGraphs(self, g1, g2, maxAtoms, flag, exttype):
        """Match two graphs (stub - not implemented)."""
        raise NotImplementedError("ccp4srs.GraphMatch.MatchGraphs is not implemented in stub")

    def GetNofMatches(self):
        """Get number of matches found (stub - returns 0)."""
        return 0

    def GetMatch(self, index):
        """Get match by index (stub - not implemented)."""
        raise NotImplementedError("ccp4srs.GraphMatch.GetMatch is not implemented in stub")

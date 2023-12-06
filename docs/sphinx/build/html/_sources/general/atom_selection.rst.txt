#######################
CCP4i2 — Atom selection
#######################

`CCP4i projects <CCP4i1Projects.html>`__

`Experimental data <data_files.html>`__

`Model data <model_data.html>`__

`Atom selection <atom_selection.html>`__

`Servers <servers.html>`__

.. container:: content

   .. rubric:: Contents
      :name: contents

   `Coordinate Ids <#cids>`__
   `Selection operators <#selection_operators>`__
   `Examples of selection commands <#examples>`__

Coordinate IDs
^^^^^^^^^^^^^^

A coordinate ID is an text string that has the following format:

::

            /mdl/chn/seq(res).ic/atm[elm]:aloc

where

======== ==============================================
**Item** **Description**
mdl      the number of the model (usually an NMR model)
chn      chain ID
seq      residue sequence number
res      residue name
ic       residue insertion code
atm      atom name
elm      chemical element ID
aloc     alternative location indicator
======== ==============================================

The residue name and atom element type are usually redundant, that is
they are not essential to uniquely identify an atom, and do not need to
be included when entering a CID. But note that all spaces in the
coordinate ID are ignored and this means that certain atom names will be
interpreted correctly only if the chemical element name is supplied
(compare Calcium CA[CA] and Carbon in alpha-position CA[C]).

Any item in the coordinate ID may be replaced by a wildcard "*", which
means an indefinite value for that item. The wildcard value is
automatically implied for any missing item except the chain ID,
insertion code and alternative location indicator:

-  A missing chain ID is interpreted as an empty ID "" if there is a
   slash "/" after it. In all other cases, the missing chain ID is
   interpreted as a wildcard "*".
-  A missing insertion code is interpreted as a wildcard "*" if the
   residue sequence number is not specified or set to a wildcard. In all
   other cases, a missing insertion code is interpreted as "" ('no
   insertion code').
-  A missing alternative location indicator is defaulted to a wildcard
   "*" if no atom name or chemical element name is given. If any value
   for atom and/or chemical element name is given , even though a
   wildcard, then a missing alternative location indicator is defaulted
   to "" ('no alternative location indicator').

A coordinate ID may be incomplete. Below are the rules for
interpretation of incomplete IDs. Curling brackets {..} denote parts of
an ID string that may be omitted as a whole:

#. If the ID starts with slash "/" then it is interpreted as /mdl { /chn
   { /seq(res).ic { /atm[elm]:aloc }}}
#. If the ID starts with a letter: chn { /seq(res).ic { /atm[elm]:aloc
   }}
#. If ID starts with a (possibly negative) number, bracket "(" or dot
   ".": seq(res).ic { /atm[elm]:aloc }
#. If none of the previous cases apply but the ID string contains square
   bracket "[" or colon ":", it is interpreted as atm[elm]:aloc

Below are examples of valid coordinate IDs:

+------------------------+--------------------------------------------+
| **Coordinate ID**      | **Description**                            |
+------------------------+--------------------------------------------+
| /1/A/33(SER).B/CA[C].A | model 1, chain A, residue SER with         |
|                        | sequence number 33 and insertion code B,   |
|                        | C-alpha atom in alternative location A.    |
+------------------------+--------------------------------------------+
| /1/A/*(SER).\*         | any SER residue in chain A, model 1.       |
+------------------------+--------------------------------------------+
| /1/A/(SER)             | any SER residue with no insertion code in  |
|                        | chain A, model 1.                          |
+------------------------+--------------------------------------------+
| /1//(SER)              | any SER residue with no insertion code in  |
|                        | chain without a chain ID, model 1.         |
+------------------------+--------------------------------------------+
| /1/A/*.*/CA[C]         | any C-alpha atom with no alternative       |
|                        | location indicator in chain A, model 1, in |
|                        | residues with any sequence number and      |
|                        | insertion code.                            |
+------------------------+--------------------------------------------+
| /1/A/*/CA[C]           | any C-alpha atom with no alternative       |
|                        | location indicator in chain A, model 1, in |
|                        | residues with any sequence number and no   |
|                        | insertion code.                            |
+------------------------+--------------------------------------------+
| A/*/CA                 | any C-alpha or Calcium atom with no        |
|                        | alternative location indicator in chain A  |
|                        | of any model, in residues with any         |
|                        | sequence number and no insertion code.     |
+------------------------+--------------------------------------------+
| CA[C]                  | any C-alpha atom with no alternative       |
|                        | location indicator in any chain, any       |
|                        | model, in residues with any sequence       |
|                        | number and no insertion code.              |
+------------------------+--------------------------------------------+
| CA                     | any atom of chain CA with no alternative   |
|                        | location indicator, in any model, in       |
|                        | residues with any sequence number and no   |
|                        | insertion code.                            |
+------------------------+--------------------------------------------+

| 

The following coordinate IDs are incorrect:

::

     1.   /A/23/CA[C]      stating with a slash implies that the model number will be given
     2.   /1/CG[C]         the residue is not defined 
     3.   /-15             stating with a slash implies that the model number will be given
     4.   */*(*).*/*[*]:*

Selection Operators
^^^^^^^^^^^^^^^^^^^

Complex selection commands are built up from CIDs, aliases and commands
connected by the operators

+--------------------+------------------------------------------------+
| *a*\ **or**\ *b*   | select all atoms defined by *a* and *b*        |
+--------------------+------------------------------------------------+
| *a*\ **and**\ *b*  | select only the atoms in both *a* and *b*      |
+--------------------+------------------------------------------------+
| **not**\ *a*       | select atoms that are not in *a*               |
+--------------------+------------------------------------------------+
| *a*\ **xor**\ *b*  | select all atoms that are in *a* or *b* but    |
|                    | not in both                                    |
+--------------------+------------------------------------------------+
| *a*\ **excl**\ *b* | select the atoms that are in *a* but exclude   |
|                    | those in *b*. This is equivalent to 'and not'. |
+--------------------+------------------------------------------------+

The selection components can also be grouped by using curly brace. Note
that curly brace are used to be distinct from the brackets used to
denote residue types in the CID.

Examples of Selection Commands
''''''''''''''''''''''''''''''

| 
| **(GLU,ASP)/(O) and neighb (HIS)/NE2,ND1 5**
| Select oxygen atoms in glutamic acid or aspartic acid residues within
  5A of the ND1 or NE2 of a histidine residue
| **{ /CA,N,C,O and {A/10-25 or A/35-40} } or (MTX)**
| Select the main chain atoms for two given ranges of residues and any
  residue with residue name MTX

.. |CCP4 web logo| image:: ccp4i2.png
   :width: 75px

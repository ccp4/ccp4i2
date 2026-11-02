"""Building the files a wwPDB OneDep deposition takes, with gemmi.

The model mmCIF a refinement program writes carries coordinates, the refinement
statistics and one entity per chain. OneDep wants, in that same file, the data
reduction statistics, one entity per crystallised molecule with its expressed
sequence, and the descriptive metadata (title, authors, source organisms,
sequence references, data collection, assembly). It ingests all of that from
the mmCIF categories themselves: the web form is pre-filled from the file, and
a category that is missing is a form page the depositor types by hand.

Modules:

- :mod:`entities` -- one entity per AU-contents sequence, ``label_seq_id`` on
  the construct.
- :mod:`mmcif_edit` -- write a category exactly once and contiguously (OneDep
  rejects a split category as "Bad model file", with no further detail).
- :mod:`metadata` -- the descriptive categories, from a plain dict.
- :mod:`struct_ref` -- ``_struct_ref``/``_struct_ref_seq``/``_struct_ref_seq_dif``
  from an alignment of the construct to a UniProt entry.
- :mod:`statistics` -- merge the scaling statistics into the model.
- :mod:`assembly` -- Matthews coefficient and the biological assembly.
- :mod:`reflections` -- the reflections mmCIF.
- :mod:`lineage` -- trace a refined model back to its scaling job.
- :mod:`onedep_validation` -- the wwPDB validation web service.
"""

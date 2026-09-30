===============================
DIALS Reciprocal Lattice Viewer
===============================

DIALS Reciprocal Lattice Viewer displays the strong spots in 3D, after mapping them from their detector positions to reciprocal space. In a favorable case, you should be able to see the crystal’s reciprocal lattice by eye in the strong spot positions. Some practice may be needed in rotating the lattice to an orientation that shows off the periodicity in reciprocal lattice positions

`DIALS JSON <https://dials.github.io/documentation/data_files.html#experiment-list-files>`_ and `Xia2 or DIALS pickle data file <https://dials.github.io/documentation/data_files.html#reflection-files>`_ need to be provided for this task.


The ability to visualize the results of spot finding in reciprocal space allows the immediate diagnosis of many indexing problems, such as:

#. Problems with image headers, such as an incorrect beam center or an inverted rotation axis, may lead to an apparent distortion in the lattice. Depending on the severity of the distortion, autoindexing may identify an incorrect lattice or result in an offset in the assigned Miller indices. 

#. Visible features that are not part of the primary lattice, such as points arranged in a spherical shell, may indicate the presence of ice rings or low-quality powder samples. 

#. Split crystals or multiple lattices are visible as a set of two or more intersecting lattices. Unindexed reflections and reflections identified as belonging to distinct lattices are colored separately to aid visualization.  

#. Multiple sweeps from a single crystal on a multi-axis goniometer can be combined for display, with each sweep uniquely colored.

In case of incorrectly defined geometry the parameters may be adjusted within the GUI, allowing common causes of failure to be easily corrected. This is valuable when commissioning a new beamline, where an accurate description of the geometry may not be available.


**References**

`Winter, G., Waterman, D. G., Parkhurst, J. M., Brewster, A. S., Gildea, R. J., Gerstel, M., Fuentes-Montero, L., Vollmar, M., Michels-Clark, T., Young, I. D., Sauter, N. K. & Evans, G. (2018). Acta Cryst. D74, 85-97. <https://doi.org/10.1107/S2059798317017235>`_
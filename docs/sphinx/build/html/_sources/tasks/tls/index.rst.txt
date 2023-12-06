+++
TLS
+++
 
 **task generate or edit TLS set definitions**


`TLS Motion Determination (TLSMD) <http://skuld.bmsc.washington.edu/~tlsmd/>`_ analyzes a macromolecular crystal structure for evidence of flexibility, e.g. local or inter-domain motions. It does this by partitioning individual chains into multiple segments that are modeled as rigid bodies undergoing TLS (Translation/Libration/Screw) vibrational (anisotropic) motion. It generates all possible partitions up to a maximum number of segments. Each trial partition is scored by how well it explains the observed atomic displacement parameters (**B values**) that came out of crystallographic refinement.

Suitable for medium resolution, when full anisotropy is impossible
Per group (20 parameters):

• **Translation** – 6 parameters

• **Libration** – 6 parameters

• **Screw rotation** – 8 parameters

    ``TLS refinement is useful when there is NCS. It is often the case that different copies of a molecule in the asymmetric unit have different overall displacements. These can be accounted for by refining TLS parameters for each molecule. The residual atomic displacement parameters (B factors) should then be similar between molecules, and NCS restraints can be applied between them.``

*TLS rigid body model is highly dependent on the quality of the crystal lattice, and conformational flexibility unique to the protein.*
 
 .. toggle-header::
        :header:  **more**
        
                 TLS (Translation/Libration/Screw) is a mathematical model that predicts the local positional displacement of atoms in a crystal structure based on an underlying assumption that each atom acts as a member of a rigid body that is displaced normally about a mean position. This displacement is seen crystallographically as non-spherical electron density at the atomic positions. The net displacement results both from actual vibration of the molecule in the crystal and from the static disorder within the crystal lattice that results from trapping different microconformers of the molecule in different unit cells making up the crystal lattice. The TLS formalism was originally developed to confirm rigid body displacement for small molecule crystallography, which yields high resolution data and precise anisotropic thermal parameters. Its success in predicting the thermal parameters for small molecules makes it tempting for use in macromolecular crystallography where it could potentially be used to find domain and loop flexibility in proteins; however, the challenges in applying the TLS model to macromolecular structures are different than those encountered with small molecules. Most macromolecular structures do not diffract to a high enough resolution to solve for individual atomic anisotropic thermal parameters (six per atom), so isotropic thermal parameters (one per atom) are used instead.


TLS Refinement tutorial can be found `here <https://www2.mrc-lmb.cam.ac.uk/groups/murshudov/content/tutorials/refmac_tutorial/files/part_2.html>`_
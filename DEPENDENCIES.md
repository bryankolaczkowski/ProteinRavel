# Required External Software and other Dependencies

ProteinRavel uses a few non-standard Python3 packages, as well as a few third-party applications that execute various parts of the analysis pipeline.

To run ProteinRavel analyses, you will need the following Python3 packages and third-party applications available. In addition, all required third-party applications must be in your executable PATH.

## Python3 Packages

### modeller

  "modeller": "^9.24"

Modeller is a Python3 language interface to build structural models of protein sequences using structural homology modeling. Modeller is available free of charge for academic or research use. Documentation and package downloads are available at the [Modeller website]. ProteinRavel uses modeller to build structural models of ancestral and extant protein sequences.

[Modeller website]: https://salilab.org/modeller/

### dendropy

  "dendropy": "^4.4.0"

Dendropy is a Python3 package for reading, manipulating and writing phylogenetic trees and associated information. Dendropy is available free of charge. Documentation and package downloads are available at the [Dendropy website]. Dendropy is also available at the [Dendropy github repository].

[Dendropy website]: https://dendropy.org/
[Dendropy github repository]: https://github.com/jeetsukumaran/DendroPy/

## Third-party Programs

Some of the third-party programs required to run ProteinRavel may be available as binary packages for your operating system, which is probably the easiest way to install and configure them. Once installed, please make sure that the required executables are in your executable PATH, before attempting to run any ProteinRavel analyses.

### mafft

  "mafft": "^v7.453"

Mafft is a multiple-sequence alignment method. ProteinRavel doesn't use mafft to align protein sequences directly, but it does use mafft to align protein sequences to a structural alignment built by Modeller, leveraging the structural alignment to guide the placement of insertions and deletions. Mafft is available free of charge. Documentation and downloads are available at the [Mafft website].

[Mafft website]: https://mafft.cbrc.jp/alignment/software/

The "mafft" executable must be in your executable PATH.

## raxml-ng (coarse)

  "raxml-ng": "^0.9.0coarse"

Raxml-ng is a program that performs maximum-likelihood phylogenetic inference, ancestral sequence reconstruction and statistical analyses of phylogenetic trees. ProteinRavel uses raxml-ng to infer protein-family phylogenies, assess uncertainty in phylogenetic inferences, and reconstruct ancestral protein sequences. Raxml-ng is available free of charge. Documentation and downloads are available on the [raxml-ng github repository].

Note that ProteinRavel makes use of coarse-grained parallelism in raxml-ng, which is currently experimental and available only in the "coarse" branch of the raxml-ng repository. Please make sure you are using the raxml-ng "coarse" branch, if you want to compile raxml-ng for use with ProteinRavel.

[raxml-ng github repository]: https://github.com/amkozlov/raxml-ng

The "raxml-ng" executable must be in your executable PATH.

### gromacs

  "gromacs": "^v2020.1"

Gromacs is a software package for performing and analyzing molecular dynamics simulations. ProteinRavel uses gromacs to generate structural topologies for energy calculations using charmm. Gromacs is available free of charge. Documentation and downloads are available on the [Gromacs website]. Gromacs is also available at the [Gromacs github repository].

Although gromacs can utilize GPU resources during molecular dynamics runs, ProteinRavel does not directly use gromacs for molecular dynamics and does not benefit from having gromacs configured to use the GPU. You do not need gromacs with GPU for ProteinRavel. However, if you want to perform mdapi analyses, you will definitely want gromacs with GPU for performing molecular dynamics simulations, in which case we recommend compiling and optimizing gromacs for the specific hardware you will be using.

[Gromacs website]: http://www.gromacs.org/
[Gromacs github repository]: https://github.com/gromacs/gromacs

The "gmx" executable must be in your executable PATH.

### charmm

  "charmm": "^44b2"

Charmm is a software package for performing molecular dynamics and similar structural analyses. ProteinRavel uses charmm to calculate potential energy fields around a protein structure using an implicit solvent model. Charmm is available free of charge (as "charmm"; there is also a paid version, but this is not needed for ProteinRavel). Documentation and downloads are available at the [Charmm website].

Although charmm can utilize GPU resources during molecular dynamics runs, ProteinRavel does not use charmm for molecular dynamics and does not benefit from having charmm configured to use GPU calculations. You do not need charmm with GPU.

[Charmm website]: https://www.charmm.org/

The "charmm" executable must be in your executable PATH.

### pdb2pqr

  "pdb2pqr": "^2.1.1"

Pdb2pqr is an application that assigns hydrogen atoms to protein structures at a specified pH, optimizes amino-acid side-chain orientations and assigns charge+radius information based on a specified atomic force-field. ProteinRavel uses pdb2pqr for structure preparation prior to potential-energy field calculations. Pdb2pqr is available free of charge. Documentation and downloads are available at the [APBS-PDB2PQR website]. Pdb2pqr is also available at the [APBS-PDB2PQR github repository].

[ABPS-PDB2PQR website]: http://www.poissonboltzmann.org/
[APBS-PDB2PQR github repository]: https://github.com/Electrostatics/apbs-pdb2pqr

The "pdb2pqr" executable must be in your executable PATH.

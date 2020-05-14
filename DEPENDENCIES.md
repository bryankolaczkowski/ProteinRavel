# Required External Software and other Dependencies

ProteinRavel uses a few non-standard Python3 packages, as well as a few third-party applications that execute various parts of the analysis pipeline.

To run ProteinRavel analyses, you will need the following Python3 packages and third-party applications available. In addition, all required third-party applications must be in your executable PATH.

## Python3 Packages

  "modeller": "^9.24"

Modeller is a Python3 language interface to build structural models of protein sequences using structural homology modeling. Modeller is available free of charge for academic or research use. Documentation and package downloads are available at the [Modeller website]. ProteinRavel uses modeller to build structural models of ancestral and extant protein sequences.

[Modeller website]: https://salilab.org/modeller/

  "dendropy": "^4.4.0"

Dendropy is a Python3 package for reading, manipulating and writing phylogenetic trees and associated information. Dendropy is available free of charge. Documentation and package downloads are available at the [Dendropy website]. Dendropy is also available at the [Dendropy github repository].

[Dendropy website]: https://dendropy.org/
[Dendropy github repository]: https://github.com/jeetsukumaran/DendroPy/

## Third-party Programs

  "mafft": "^v7.453"

  "raxml-ng": "^0.9.0coarse"

  "gromacs": "^v2020.1"

  "charmm": "^44b2"

  "pdb2pqr": "^2.1.1"

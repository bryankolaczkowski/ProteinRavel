# Required External Software and other Dependencies

ProteinRavel uses a few non-standard Python3 packages, as well as a few third-party applications that execute various parts of the analysis pipeline. 

To run ProteinRavel analyses, you will need the following Python3 packages available.

## Python3 Packages

  "modeller": "^9.24"

Modeller is a Python3 language interface to build structural models of protein sequences using structural homology modeling. Modeller is available free of charge for academic or research use. Documentation and package downloads are available at the [modeller-reason]. ProteinRavel uses modeller to build structural models of ancestral and extant protein sequences.

[modeller-reason]: https://salilab.org/modeller/ "Modeller website"

  "dendropy": "^4.4.0"

Dendropy is a Python3 package for reading, manipulating and writing phylogenetic trees and associated information. Dendropy is available free of charge. Documentation and package downloads are available at the [dendropy-reason]. Dendropy also has a [dendropy-git].

[dendropy-reason]: https://dendropy.org/ "Dendropy website"
[dendropy-git]: https://github.com/jeetsukumaran/DendroPy/ "github repository"

## Third-party Programs

  mafft >= v7.453 (2019/Nov/8)
  raxml-ng >= 0.9.0coarse
  gromacs >= v2020.1
  charmm >= 44b2
  pdb2pqr >= 2.1.1
  

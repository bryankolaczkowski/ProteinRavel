# ProteinRavel
A collection of Python3 modules for characterizing the evolution of protein family sequence, structure and function

## Installing ProteinRavel

ProteinRavel has a few Python3 and external third-party dependencies, which are described in the [DEPENDENCIES.md] file. Please make sure you have all the required dependencies installed and functioning properly, before attempting to install and use ProteinRavel.

[DEPENDENCIES.md]: DEPENDENCIES.md

Other than the required Python3 and external third-party dependencies, ProteinRavel is pretty self-contained and easy to set up. The main executable program is:

bin/ravel

The required Python3 library is in proteinRavel/, and required data files and other information is in data/

So long as the directory hierarchy is maintained, you can put the ProteinRavel directory just about anywhere on your system, and it should work. You can also feel free to create a symbolic link to the bin/ravel executable, if you don't want to change your executable PATH.

For an easy system-wide installation on a typical 'NIX system, we recommend:

  cp -r ProteinRavel /usr/local/
  chown -r root:root /usr/local/ProteinRavel
  ln -s /usr/local/ProteinRavel/bin/ravel /usr/local/bin/

## Running ProteinRavel

Once installed, ProteinRavel can be executed from the command-line:

  ravel -h

Will print basic usage information to the screen.

### Common Options

Ravel uses multithreading in all analysis steps to speed up computation. By default, ravel is quite thread-greedy, using all available CPU cores. You can change the number of threads to a specified number using the "-t" option.

A ProteinRavel analysis produces lots of intermediate files. By default, ravel's output files are prefixed with "ravel-", but you can change the output file prefix using the "-p" option.

### Analysis Steps

A ProteinRavel analysis consists of multiple sequential steps, each of which should be thoroughly tested and examined before running the next step in the analysis. Although we have done our best to produce a tool that should be able to analyze most protein families reliably, default parameter values may not be applicable to your protein family of interest, and most real protein families will require some manual curation and parameter tuning to produce reliable results. Ravel is not intended as a 'fire and forget' type of software, and we do not recommend using it in this way.

The typical steps in a ProteinRavel analysis are (in order):

1. align - Align template structures, then align protein sequences to the structural alignment. This step also removes protein sequences with large insertions or deletions (relative to the aligned structures); large indels produce structural modeling artifacts that can ruin downstream analyses. Duplicate protein sequences are also removed to save computation time.

2. tree - Infer a maximum-likelihood phylogeny of the aligned protein family sequences. This step also assesses statistical confidence in each node on the phylogeny and generates a set of plausible candidate trees capturing uncertainty in the protein family tree.

3. asr - Reconstruct ancestral protein sequences. This step incorporates phylogenetic uncertainty to reconstruct probabilistic ancestral sequences at strongly-supported ancestral nodes on the protein family tree.

4. smodel - Infer replicate structural models of ancestral and extant protein family sequences. Uncertainty in ancestral sequence reconstruction (and phylogeny) is incorporated by randomly drawing replicate ancestral sequences from the (squashed and truncated) posterior probability distributions at each site. Uncertainty in structural modeling is incorporated by stochastic replicate structural modeling.

5. qmap - Calculate potential energy fields around replicate structural models of ancestral and extant protein sequences. This step generates grid-based electrostatic and hydrophobic potential-energy fields around each protein structural model. Although a simplification (compared to full molecular dynamics), this is probably the most computation-consuming step in the ProteinRavel analysis pipeline. As this step could take quite some time for large phylogenies, we recommend benchmarking your analysis using a small data set before launching a large-scale analysis.

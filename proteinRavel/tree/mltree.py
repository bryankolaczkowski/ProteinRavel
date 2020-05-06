################################################################################
# protein family phylogeny inference with bootstrap support calculation
#
# (C) 2020 Bryan Kolaczkowski, University of Florida, Gainesville, FL USA
# Released under GNU General Public License (GPL)
# bryank@ufl.edu
################################################################################

import sys
import os
import tempfile
import multiprocessing

################################################################################
# BEG VARIABLE DEFINITIONS

dnameprefix = 'mltree_'  # prefix to use for temporary directories
fnameprefix = 'mltree'   # prefix to use for temporary files

# END VARIABLE DEFINITIONS
################################################################################

################################################################################
# BEG CLASS DEFINITIONS

# END CLASS DEFINITIONS
################################################################################

################################################################################
# BEG HELPER FUNCTION DEFINITIONS

def run_raxml(raxml_exe, raxml_opts, raxml_threads, raxml_io, raxml_redirect,
              verbose):
    raxml_cmd = '{} {} {} {} {}'.format(raxml_exe, raxml_opts, raxml_threads,
                                        raxml_io, raxml_redirect)
    if verbose:
        sys.stdout.write('{}\n'.format(raxml_cmd))
    os.system(raxml_cmd)
    return

# END HELPER FUNCTION DEFINITIONS
################################################################################

################################################################################
# BEG MODULE INTERFACE

def mltree(alnfname, model='LG+G8+F', bootstraps=200, boot_metric='tbe',
           threads=0, verbose=True):
    """
    Infers a maximum-likelihood phylogeny.

    Infers the maximum-likelihood protein family phylogeny for aligned
    protein sequences, assuming the input model. Also calculates statistical
    support for the phylogeny.

    Parameters:
      alnfname    (string): aligned protein sequences in FASTA format
      model       (string): evolutionary model to use
      bootstraps (integer): number of bootstrap replicates
      boot_metric (string): bootstrap metric to use
      threads    (integer): Number of threads to run for inference
      verbose    (boolean): Show function execution to stdout
    Returns:
      Newick-formatted tree (string) with branch lengths and bootstrap
      proportions
    """

    final_tree_string = ''  # return value -> ML tree

    # set threads to max available if not set manually
    if threads <= 0:
        threads = multiprocessing.cpu_count()

    # set up temporary directory for raxml output
    with tempfile.TemporaryDirectory(prefix=dnameprefix) as tempdir:
        # set up raxml run
        raxml_prefix = tempdir + os.path.sep + fnameprefix

        if verbose:
            raxml_redirect = ''
        else:
            raxml_redirect = '> /dev/null'

        raxml_exe  = 'raxml-ng'
        raxml_io   = '--msa {} --model {} --prefix {}'.format(alnfname,
                                                              model,
                                                              raxml_prefix)
        raxml_thrd = '--threads {} --workers {}'.format(threads, threads)

        # execute raxml maximum-likelihood inference
        if verbose:
            sys.stdout.write('BEG Phylogeny Search\n')
        raxml_opts = '--search --tree pars{{{}}},rand{{{}}}'.format(threads,
                                                                    threads)
        run_raxml(raxml_exe, raxml_opts, raxml_thrd, raxml_io, raxml_redirect,
                  verbose)
        if verbose:
            sys.stdout.write('END Phylogeny Search\n')

        # execute raxml bootstrap run
        if verbose:
            sys.stdout.write('BEG Support Estimation\n')
        raxml_opts = '--bootstrap --bs-trees {}'.format(bootstraps)
        run_raxml(raxml_exe, raxml_opts, raxml_thrd, raxml_io, raxml_redirect,
                  verbose)
        if verbose:
            sys.stdout.write('END Support Estimation\n')

        # calculate bootstrap support
        if verbose:
            sys.stdout.write('BEG Support Integration\n')
        raxml_ml_fname = raxml_prefix + '.raxml.bestTree'
        raxml_bs_fname = raxml_prefix + '.raxml.bootstraps'
        raxml_io = '--tree {} --bs-trees {} --prefix {}'.format(raxml_ml_fname,
                                                                raxml_bs_fname,
                                                                raxml_prefix)
        raxml_opts = '--support --bs-metric {}'.format(boot_metric)
        run_raxml(raxml_exe, raxml_opts, raxml_thrd, raxml_io, raxml_redirect,
                  verbose)

        # read ML tree with bootstrap support
        raxml_supptree_fname = raxml_prefix + '.raxml.support'
        with open(raxml_supptree_fname, 'r') as handle:
            for line in handle:
                final_tree_string += line.strip()
        if verbose:
            sys.stdout.write('END Support Integration\n')

    return final_tree_string

# END MODULE INTERFACE
################################################################################

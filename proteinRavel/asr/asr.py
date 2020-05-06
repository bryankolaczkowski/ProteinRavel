################################################################################
# maximum-likelihood ancestral sequence reconstruction with indels
#
# (C) 2020 Bryan Kolaczkowski, University of Florida, Gainesville, FL USA
# Released under GNU General Public License (GPL)
# bryank@ufl.edu
################################################################################

import sys
import os
import math
import tempfile
import multiprocessing
import collections
import dendropy

import proteinRavel.align

################################################################################
# BEG VARIABLE DEFINITIONS

dnameprefix = 'asr_'  # prefix to use for temporary directories
fnameprefix = 'asr'   # prefix to use for temporary files

# END VARIABLE DEFINITIONS
################################################################################

################################################################################
# BEG CLASS DEFINITIONS

class AncestralSequence(proteinRavel.align.Sequence):
    """Probabilistic Ancestral Protein Sequence."""

    ambigchar = 'X'  # ambiguous amino-acid residue

    def __init__(self, filename, node_id):
        """Initialize a new AncestralSequence object from a raxml-ng file."""
        super().__init__(node_id, '')
        self.indelchars = [self.abschar, self.ambigchar]
        self.indelprobs = []
        self.aaprobs    = []
        with open(filename, 'r') as handle:
            header    = handle.readline().split()
            id_col    = header.index('Node')
            site_col  = header.index('Site')
            state_col = header.index('State')
            prob_cols = state_col+1
            self.aachars = [ x.split('p_')[-1] for x in header[prob_cols:] ]
            for line in handle:
                linearr = line.split()
                id = linearr[id_col]
                if id == node_id:
                    self.sequence += linearr[state_col]
                    self.indelprobs.append( [0.0, 1.0] )
                    self.aaprobs.append([ float(x) for x in linearr[prob_cols:] ])
        return

    def add_indels(self, filename):
        """Add indel probabilities and sequence gaps from a raxml-ng file."""
        all_missing_chars = '0-'
        with open(filename, 'r') as handle:
            header    = handle.readline().split()
            id_col    = header.index('Node')
            site_col  = header.index('Site')
            state_col = header.index('State')
            prob_cols = state_col+1
            rawchars  = [ x.split('p_')[-1] for x in header[prob_cols:] ]
            missing_index = rawchars.index('0')
            present_index = rawchars.index('1')
            for line in handle:
                linearr = line.split()
                id = linearr[id_col]
                if id == self.identifier:
                    site_idx = int(linearr[site_col]) - 1
                    if linearr[state_col] in all_missing_chars:
                        self.sequence = self.sequence[:site_idx]  \
                                        + self.abschar            \
                                        + self.sequence[site_idx+1:]
                    iprobs = [ float(x) for x in linearr[prob_cols:] ]
                    self.indelprobs[site_idx] = [ iprobs[missing_index],
                                                  iprobs[present_index] ]
        return

    def fullASR(self):
        """AncestralSequence as a full ASR report string."""
        retstr = self.PIR()
        retstr += str(self.indelchars) + ' ' + str(self.aachars) + '\n'
        for i in range(self.__len__()):
            retstr += str(i) + ' ' + str(self.indelprobs[i]) + ' ' + \
                                     str(self.aaprobs[i]) + '\n'
        return retstr

    def _aa_brier(self, idx):
        """returns a plausible_chars, probabilities pair for sequence position idx."""
        best_score = math.inf
        best_i     = 0
        aaprobs = []
        for i in range(len(self.aachars)):
            aaprobs.append((self.aaprobs[idx][i], i))
        aaprobs.sort(reverse=True)
        for i in range(len(aaprobs)):
            score = sum([math.pow(x[0]-(1.0/(i+1)), 2) for x in aaprobs[:i+1]]) + \
                    sum([math.pow(x[0],             2) for x in aaprobs[i+1:]])
            if score < best_score:
                best_score = score
                best_i     = i
        all_chars = [ self.aachars[x[1]] for x in aaprobs[:best_i+1] ]
        all_probs = [              x[0]  for x in aaprobs[:best_i+1] ]
        return (all_chars, all_probs)


    def brier(self):
        """Apply Brier-based ambiguity criterion to calculate sequence ambiguity."""
        retstr = '{}{}\n'.format(self.idprefix, self.identifier)
        for i in range(self.__len__()):
            retstr += '{} '.format(i)
            indelp = self.indelprobs[i]
            if indelp[0] > indelp[1]:
                # gap is ml state
                pchars = [ self.indelchars[0] ]
                pprobs = [ indelp[0] ]
                b_gap = math.pow(indelp[0]-1.0, 2) + \
                        math.pow(indelp[1]    , 2)
                b_amb = math.pow(indelp[0]-0.5, 2) + \
                        math.pow(indelp[1]-0.5, 2)
                if b_gap >= b_amb:
                    # anc state could be gap or an amino-acid state
                    chars,probs = self._aa_brier(i)
                    pchars += chars
                    pprobs += probs
            else:
                # gap is not ml state
                pchars, pprobs = self._aa_brier(i)
                b_prs = math.pow(indelp[1]-1.0, 2) + \
                        math.pow(indelp[0],     2)
                b_amb = math.pow(indelp[0]-0.5, 2) + \
                        math.pow(indelp[1]-0.5, 2)
                if b_prs >= b_amb:
                    # gap or amino-acid state plausible
                    pchars += [ self.indelchars[0] ]
                    pprobs += [ indelp[0] ]
            retstr += ''.join(pchars) + ' '
            retstr += ','.join([ '{:.4f}'.format(round(x,4)) for x in pprobs ])
            retstr += '\n'
            #retstr += str(self.indelprobs[i]) + ' ' + str(self.aaprobs[i]) + '\n'
        return retstr

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

def _reconstruct_ancseqs(aln_fasta_string, tree,
                         raxml_alnfname, raxml_trefname, raxml_prefix,
                         model, threads, verbose):

    # write alignment file
    with open(raxml_alnfname, 'w') as handle:
        handle.write(aln_fasta_string)
    # write tree file
    tree.write(path=raxml_trefname, schema='newick')

    # set up raxml run
    if verbose:
        raxml_redirect = ''
    else:
        raxml_redirect = '> /dev/null'

    raxml_exe  = 'raxml-ng'
    raxml_io   = '--msa {} --tree {} --prefix {}'.format(raxml_alnfname,
                                                         raxml_trefname,
                                                         raxml_prefix)
    raxml_thrd = '--force --threads {}'.format(threads)
    raxml_opts = '--ancestral --blopt nr_safe --lh-epsilon 0.01 --model {}'.format(model)
    # execute raxml asr
    run_raxml(raxml_exe, raxml_opts, raxml_thrd, raxml_io, raxml_redirect,
              verbose)

    return raxml_prefix + '.raxml.ancestralProbs'

# END HELPER FUNCTION DEFINITIONS
################################################################################

################################################################################
# BEG MODULE INTERFACE

def asr(alnfname, treefname, model='LG+G8+F', consensus=80,
        threads=0, verbose=True):
    """
    Reconstructs ancestral sequences.

    Infers maximum-likelihood ancestral protein sequences, assuming the
    aligned sequences, phylogenetic tree and the input model. Nodes with
    < consensus support are excluded.

    Parameters:
      alnfname   (string): aligned protein sequences in FASTA format
      treefname  (string): phylogenetic tree in Newick format
      model      (string): evolutionary model to use
      consensus (integer): collapse nodes with < this support
      threads   (integer): Number of threads to run for inference
      verbose   (boolean): Show function execution to stdout
    Returns:
      (tree, ancestral_sequences) pair
    """

    # set threads to max available if not set manually
    if threads <= 0:
        threads = multiprocessing.cpu_count()

    # read tree object
    if verbose:
        sys.stdout.write('BEG Tree Organization\n')
    tree = dendropy.Tree.get(path=treefname, schema='newick',
                             extract_comment_metadata=False,
                             case_sensitive_taxon_labels=True,
                             preserve_underscores=True)

    # assign ancestral node labels to tree
    nodecount = len(tree.internal_nodes(exclude_seed_node=True))
    total_zeros = math.floor(math.log10(nodecount))
    anc_node_prefix = 'n'
    anc_node_number = 1
    anc_sequences   = collections.OrderedDict()
    for node in tree.preorder_internal_node_iter(exclude_seed_node=True):
        supp_val = float(node.label)
        # convert support values to percent for consensus comparison
        if supp_val < 1.1:
            supp_val *= 100
        if supp_val > 100.0:
            supp_val = 100.0
        # exclude nodes with support < consensus
        if supp_val >= float(consensus):
            zeros  = total_zeros - math.floor(math.log10(anc_node_number))
            nlabel = anc_node_prefix + '0'*zeros + str(anc_node_number)
            node.label = nlabel
            anc_sequences[nlabel] = None
            anc_node_number += 1
        else:
            node.label = ''
    if verbose:
        sys.stdout.write('found {} supported ancestral nodes at {}% consensus\n'.format(anc_node_number-1, consensus))
        sys.stdout.write('END Tree Organization\n')

    # parse alignment file
    if verbose:
        sys.stdout.write('BEG Alignment Conversion\n')
    aln = proteinRavel.align.Alignment(alnfname)
    protein_aln_fasta = aln.FASTA()
    binary_aln_fasta  = aln.binary()
    if verbose:
        sys.stdout.write('END Alignment Conversion\n')

    # reconstruct ancestral sequences
    if verbose:
        sys.stdout.write('BEG Protein Sequence Reconstruction\n')

    # set up temporary directory for raxml output
    with tempfile.TemporaryDirectory(prefix=dnameprefix) as tempdir:
        # write alignment and tree files
        raxml_alnfname = tempdir + os.path.sep + fnameprefix + '.fasta'
        raxml_trefname = tempdir + os.path.sep + fnameprefix + '.tre'
        raxml_prefix   = tempdir + os.path.sep + fnameprefix

        # run raxml-ng to reconstruct ancestral sequences
        raxml_outfname = _reconstruct_ancseqs(protein_aln_fasta, tree,
                                              raxml_alnfname, raxml_trefname,
                                              raxml_prefix, model, threads,
                                              verbose)

        # parse raxml asr output
        for nodeid in anc_sequences.keys():
            anc_sequences[nodeid] = AncestralSequence(raxml_outfname, nodeid)

    if verbose:
        sys.stdout.write('END Protein Sequence Reconstruction\n')

    # reconstruct ancestral indels
    if verbose:
        sys.stdout.write('BEG Indel Reconstruction\n')

    # set up temporary directory for raxml output
    with tempfile.TemporaryDirectory(prefix=dnameprefix) as tempdir:
        # write alignment and tree files
        raxml_alnfname = tempdir + os.path.sep + fnameprefix + '.fasta'
        raxml_trefname = tempdir + os.path.sep + fnameprefix + '.tre'
        raxml_prefix   = tempdir + os.path.sep + fnameprefix

        # run raxml-ng to reconstruct ancestral sequences
        binmodel = 'BIN+G8+F'
        raxml_outfname = _reconstruct_ancseqs(binary_aln_fasta, tree,
                                              raxml_alnfname, raxml_trefname,
                                              raxml_prefix, binmodel, threads,
                                              verbose)

        # parse raxml asr output
        for nodeid in anc_sequences.keys():
            anc_sequences[nodeid].add_indels(raxml_outfname)

    if verbose:
        sys.stdout.write('END Indel Reconstruction\n')

    return (tree, anc_sequences)

# END MODULE INTERFACE
################################################################################

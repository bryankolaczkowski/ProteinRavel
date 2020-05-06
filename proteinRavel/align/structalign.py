################################################################################
# structure-guided multiple sequence alignment with sequence trimming and
# filtering options
#
# (C) 2020 Bryan Kolaczkowski, University of Florida, Gainesville, FL USA
# Released under GNU General Public License (GPL)
# bryank@ufl.edu
################################################################################

import sys
import os
import re
import collections
import tempfile
import multiprocessing
import modeller
import modeller.salign

################################################################################
# BEG VARIABLE DEFINITIONS

fnameprefix = 'structalign_'  # prefix to use for temporary files

# END VARIABLE DEFINITIONS
################################################################################

################################################################################
# BEG CLASS DEFINITIONS

class AlignmentError(Exception):
    """Exceptions raised when alignment has an issue."""
    pass


class Sequence(object):
    """Protein sequence."""

    idprefix  = '>'  # comes before FASTA and PIR identifiers
    piridsep  = ';'  # between pirid and sequence id in PIR format
    pirlinsep = ':'  # between elements in the PIR line
    pirseqend = '*'  # after sequence in PIR format
    preschar  = '.'  # an amino-acid is present
    abschar   = '-'  # no amino-acid is present (ie, a gap)

    def __init__(self, id, seq, pirid=None, pirline=None):
        """Initialize a new Sequence object."""
        self.identifier  = id
        self.sequence    = seq
        if pirid:
            self.pirid   = pirid
        else:
            self.pirid   = 'P1'
        if pirline:
            self.pirline = pirline
        else:
            self.pirline = 'sequence:::::::::'
        return

    def __str__(self):
        """String representation in PIR format."""
        return self.PIR()

    def __repr__(self):
        return self.__str__()

    def __len__(self):
        """Length of the amino-acid sequence, including gaps."""
        return len(self.sequence)

    def FASTA(self):
        """Sequence as a FASTA-formatted string."""
        return self.idprefix + self.identifier + '\n' + \
               self.sequence +                   '\n'

    def PIR(self):
        """Sequence as a PIR-formatted string, just as modeller likes it."""
        return self.idprefix   + self.pirid + self.piridsep + \
               self.identifier +                       '\n' + \
               self.pirline    +                       '\n' + \
               self.sequence   + self.pirseqend +      '\n'

    def binary(self):
        """
        Sequence as a binary presence-absence string.

        1 = amino acid present
        0 = no amino acid present
        """
        pas     = self.pres_abs()
        bin_seq = pas.translate(pas.maketrans(self.preschar + self.abschar,
                                              '10'))
        return self.idprefix + self.identifier + '\n' + \
               bin_seq +                         '\n'

    def pres_abs(self):
        """Returns a presence-absence representation of this sequence."""
        return re.sub(r'[a-zA-Z]', self.preschar, self.sequence)


class Alignment(Sequence):
    """Alignment of Sequences."""

    def __init__(self, filename, includeprefix=None, excludeprefix=None):
        """Initialize a new Alignment object from a FASTA or PIR file."""
        self.sequencelist = []
        self.length       = 0
        self.initFromFile(filename, includeprefix, excludeprefix)
        return

    def initFromFile(self, filename, includeprefix=None, excludeprefix=None):
        """Initialize (or re-initialize) this alignment from a file."""
        current_id_list = [ x.identifier for x in self.sequencelist ]
        self.length     = 0
        with open(filename, 'r') as handle:
            line = handle.readline()
            while line:
                if line[0] == self.idprefix:
                    id  = line[1:].strip()
                    pid = ''
                    # check for include prefix, jump out early?
                    if includeprefix:
                        if id.find(includeprefix) != 0:
                            line = handle.readline()
                            continue
                        else:
                            id = id.replace(includeprefix, '', 1)
                    # check for exclude prefix, jump out early?
                    if excludeprefix:
                        if id.find(excludeprefix) == 0:
                            line = handle.readline()
                            continue
                    # separate PIR id if found
                    if self.piridsep in id:
                        pid = id.split(self.piridsep)[0]
                        id  = id.split(self.piridsep)[1]
                    pil = ''
                    se  = ''
                    line = handle.readline()
                    # grab PIR line, if found
                    if self.pirlinsep in line:
                        pil  = line.strip()
                        line = handle.readline()
                    while line and line[0] != self.idprefix:
                        se += line.strip()
                        line = handle.readline()
                    # remove PIR sequence end, if found
                    if se[-1] == self.pirseqend:
                        se = se[:-1]
                    # set length, if needed
                    if self.length == 0:
                        self.length = len(se)
                    # check length
                    if len(se) != self.length:
                        raise AlignmentError('sequence {} has length {}; alignment has length {}'.format(newseq.identifier, len(newseq), self.length))
                    # create new sequences
                    if not current_id_list:
                        newseq = Sequence(id, se, pid, pil)
                        self.sequencelist.append(newseq)
                    # adding to an existing alignment; only change the sequence
                    else:
                        try:
                            idx = current_id_list.index(id)
                            self.sequencelist[idx].sequence = se
                        except:
                            raise AlignmentError('trying to add a new sequence {} to an existing alignment'.format(id))
                # skip over blank lines and other stuff we don't want
                else:
                    line = handle.readline()
        return

    def __str__(self):
        """String representation in PIR format."""
        return self.PIR()

    def __repr__(self):
        return self.__str__()

    def __len__(self):
        """Number of characters in this alignment."""
        return self.length

    def size(self):
        """Number of Sequences in this alignment."""
        return len(self.sequencelist)

    def FASTA(self):
        """Alignment as FASTA-formatted string."""
        aln = ''
        for seq in self.sequencelist:
            aln += seq.FASTA()
        return aln

    def PIR(self):
        """Alignment as PIR-formatted string, just as modeller likes it."""
        aln = ''
        for seq in self.sequencelist:
            aln += seq.PIR()
        return aln

    def binary(self):
        """
        Alignment as a binary presence-absence string.

        In FASTA format
        1 = amino acid present
        0 = no amino acid present
        """
        aln = ''
        for seq in self.sequencelist:
            aln += seq.binary()
        return aln

    def pres_abs(self):
        """Returns a presence-absence representation of this alignment"""
        pas = self.abschar*self.length
        for s_pas in [ x.pres_abs() for x in self.sequencelist ]:
            for i in range(len(s_pas)):
                if s_pas[i] == self.preschar:
                    pas = pas[:i] + self.preschar + pas[i+1:]
        return pas

    def trim(self, other):
        """
        Trims begining and end of self and other, removing characters that
        are all gaps in self
        """
        pa_seq = self.pres_abs()
        bindex = 0
        while pa_seq[bindex] == self.abschar:
            bindex += 1
        eindex = len(pa_seq) - 1
        while pa_seq[eindex] == self.abschar:
            eindex -= 1
        eindex += 1
        pa_seq = pa_seq[bindex:eindex]
        self.length  = len(pa_seq)
        other.length = len(pa_seq)
        for seq in self.sequencelist:
            seq.sequence = seq.sequence[bindex:eindex]
        for seq in other.sequencelist:
            seq.sequence = seq.sequence[bindex:eindex]
        return

    def filter(self, other, maxindel, verbose=False):
        """
        Filters the other alignment, removing any sequences that have insertions
        or deletions > maxindel.
        """

        pa_seq = self.pres_abs()
        seq_ids_removed = set([])
        seqs_removed = set([])

        # get sequences with large inserts
        gappattern = re.compile('\%s{%d,}' % (self.abschar, maxindel+1))
        for match in gappattern.finditer(pa_seq):
            beg,end = match.span()
            end += 1
            # for each match, identify sequences with large inserts
            for i in range(len(other.sequencelist)):
                insertlen = len(other.sequencelist[i].sequence[beg:end].replace(other.abschar,''))
                if insertlen > maxindel:
                    seqs_removed.add(i)
                    seq_ids_removed.add(other.sequencelist[i].identifier)

        # get sequences with large deletions
        prespattern = re.compile('\%s{%d,}' % (self.preschar, maxindel+1))
        for match in prespattern.finditer(pa_seq):
            beg,end = match.span()
            end += 1
            # for each match, identify sequences with large deletions
            for i in range(len(other.sequencelist)):
                for gapmatch in gappattern.finditer(other.sequencelist[i].sequence[beg:end]):
                    del_len = (match.end()-match.start()) + 1
                    if del_len > maxindel:
                        seqs_removed.add(i)
                        seq_ids_removed.add(other.sequencelist[i].identifier)

        # actually remove the sequences
        for i in sorted(list(seqs_removed), reverse=True):
            del other.sequencelist[i]
        return seq_ids_removed

    def removeColumns(self, columns):
        """Removes a list of column indices from this alignment."""
        self.length = self.length - len(columns)
        for col in columns:
            for seq in self.sequencelist:
                seq.sequence = seq.sequence[:col] + seq.sequence[col+1:]
        return

    def removeIdenticalSequences(self):
        """Removes identical sequences from this alignment."""
        seqs_found = set([])
        indices_to_remove = []
        for i in range(len(self.sequencelist)):
            if self.sequencelist[i].sequence in seqs_found:
                indices_to_remove.append(i)
            else:
                seqs_found.add(self.sequencelist[i].sequence)
        for i in sorted(indices_to_remove, reverse=True):
            del self.sequencelist[i]
        return

# END CLASS DEFINITIONS
################################################################################

################################################################################
# BEG HELPER FUNCTION DEFINITIONS

def _remove_gap_only(seqaln, straln):
    """Removes columns tat are gap-only in both alignments."""
    cols_to_remove = []
    pa_seq = seqaln.pres_abs()
    pa_str = straln.pres_abs()
    for i in range(len(pa_seq)):
        if pa_seq[i] == seqaln.abschar and pa_str[i] == straln.abschar:
            cols_to_remove.append(i)
    seqaln.removeColumns(cols_to_remove)
    straln.removeColumns(cols_to_remove)
    return


def _align_structures(structures, verbose):
    """Aligns structures using iterative structural alignment."""

    # set up modeller environment
    if verbose:
        modeller.log.verbose()
    else:
        modeller.log.none()
    env = modeller.environ()
    aln = modeller.alignment(env)

    # read structures into modeller environment
    for (id,structure) in structures.items():
        mdl = modeller.model(env, file=structure)
        aln.append_model(mdl, align_codes=id, atom_files=structure)

    # align structures using iterative structural alignment
    modeller.salign.iterative_structural_align(aln)

    # convert modeller alignment to Alignment object
    mod_aln_f = tempfile.NamedTemporaryFile(mode='w',
                                            prefix=fnameprefix,
                                            suffix='.ali',
                                            delete=False)
    mod_aln_fname = mod_aln_f.name
    mod_aln_f.close()
    aln.write(mod_aln_fname, alignment_format='PIR')
    alnobj = Alignment(mod_aln_fname)
    os.remove(mod_aln_fname)
    return alnobj


def _align_sequences(sequencefname, struct_aln, threads, verbose):
    """
    Aligns sequences guided by aligned structures.

    SIDE EFFECT: struct_aln is aligned to the sequence alignment that is
                 returned
    """

    # write aligned structures as FASTA seed to temporary file
    struct_aln_seedf = tempfile.NamedTemporaryFile(mode='w',
                                                   prefix=fnameprefix,
                                                   suffix='.fasta',
                                                   delete=False)
    struct_aln_seedf_name = struct_aln_seedf.name
    struct_aln_seedf.write('{}\n'.format(struct_aln.FASTA()))
    struct_aln_seedf.close()

    # set up mafft run
    if verbose:
        mafft_redirect = ''
    else:
        mafft_redirect = '2> /dev/null'
    mafft_exe   = 'mafft'
    mafft_opts  = '--localpair --maxiterate 1000 --ep 0.123 --anysymbol'
    mafft_thrds = '--thread {} --threadtb {} --threadit {}'.format(threads,
                                                                   threads,
                                                                   threads)
    mafft_seed  = '--seed {}'.format(struct_aln_seedf_name)

    mafft_aln_outf = tempfile.NamedTemporaryFile(mode='w',
                                                 prefix=fnameprefix,
                                                 suffix='.fasta',
                                                 delete=False)
    mafft_aln_outf_name = mafft_aln_outf.name
    mafft_aln_outf.close()

    # execute mafft run
    mafft_cmd = '{} {} {} {} {} {} > {}'.format(mafft_exe,
                                                mafft_opts,
                                                mafft_seed,
                                                mafft_thrds,
                                                sequencefname,
                                                mafft_redirect,
                                                mafft_aln_outf_name)
    if verbose:
        sys.stdout.write('{}\n'.format(mafft_cmd))

    os.system(mafft_cmd)

    # parse mafft alignment into seed- and non-seed sequences
    mafft_seed_ident = '_seed_'  # prefix mafft uses to identify seed sequences
    struct_aln.initFromFile(mafft_aln_outf_name, includeprefix=mafft_seed_ident)
    seq_aln = Alignment(mafft_aln_outf_name, excludeprefix=mafft_seed_ident)

    # clean up all temporary files we made
    os.remove(struct_aln_seedf_name)
    os.remove(mafft_aln_outf_name)

    return seq_aln

# END HELPER FUNCTION DEFINITIONS
################################################################################

################################################################################
# BEG MODULE INTERFACE

def structalign(sequence_fname, structures_fnames,
                trim=True, filter=True, maxindel=8,
                threads=0, verbose=True):
    """
    Aligns input structures and sequences.

    Aligns input structures based on iterative structural superposition. That
    structural alignment is then used to guide the alignment of the input
    sequences. Optionally, the beginning and ending of sequences that don't
    overlap the structural alignment can be trimmed, and sequences that contain
    large internal insertions or deletions - relative to the structural
    alignment - can be removed.

    Parameters:
      sequence_fname        (string): unaligned FASTA filename
      structures_fnames (dictionary): ID:protein structure filename
      trim                 (boolean): Remove beginning and ending parts of
                                    sequences not aligning to the structure
                                    alignment
      filter               (boolean): Remove sequences with large indels
      maxindel             (integer): Largest indel to accept
      threads              (integer): Number of threads to run for alignment
      verbose              (boolean): Show function execution to stdout
    Returns:
      (aligned_sequences, aligned_structures) pair
    """

    # set threads to max available if not set manually
    if threads <= 0:
        threads = multiprocessing.cpu_count()

    # align structures
    if verbose:
        sys.stdout.write('BEG Structure Alignment\n')
    aln_struct = _align_structures(structures_fnames, verbose)
    if verbose:
        #sys.stdout.write(str(aln_struct))
        sys.stdout.write('END Structure Alignment\n')

    # align sequences
    if verbose:
        sys.stdout.write('BEG Sequence Alignment\n')
    aln_seq = _align_sequences(sequence_fname, aln_struct, threads, verbose)
    if verbose:
        #sys.stdout.write(str(aln_struct))
        #sys.stdout.write(str(aln_seq))
        sys.stdout.write('END Sequence Alignment\n')

    # trim sequences
    if trim:
        if verbose:
            sys.stdout.write('BEG Trim Alignment\n')
            sys.stdout.write('raw alignment is {} characters long\n'.format(aln_struct.length))

        aln_struct.trim(aln_seq)

        if verbose:
            sys.stdout.write('trimmed alignment is {} characters long\n'.format(aln_struct.length))
            #sys.stdout.write(str(aln_struct))
            #sys.stdout.write(str(aln_seq))
            sys.stdout.write('END Trim Alignment\n')

    # filter sequences
    if filter:
        if verbose:
            sys.stdout.write('BEG Filter Alignment\n')
            sys.stdout.write('sequence alignment has {} sequences\n'.format(aln_seq.size()))
        seqs_removed = aln_struct.filter(aln_seq, maxindel)
        _remove_gap_only(aln_seq, aln_struct)
        if verbose:
            sys.stdout.write('removed sequences: {}\n'.format(seqs_removed))
            sys.stdout.write('filtered sequence alignment has {} sequences\n'.format(aln_seq.size()))
            #sys.stdout.write(str(aln_struct))
            #sys.stdout.write(str(aln_seq))
            sys.stdout.write('END Filter Alignment\n')

    # remove identical sequences from the sequence alignment
    if verbose:
        sys.stdout.write('BEG Remove identical sequences\n')
    orig_aln_size = aln_seq.size()
    aln_seq.removeIdenticalSequences()
    if verbose:
        sys.stdout.write('removed {} identical sequences\n'.format(orig_aln_size-aln_seq.size()))
        sys.stdout.write('sequence alignment now has {} sequences\n'.format(aln_seq.size()))
        sys.stdout.write('END Remove identical sequences\n')

    return (aln_seq, aln_struct)

# END MODULE INTERFACE
################################################################################

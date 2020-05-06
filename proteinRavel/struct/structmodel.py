################################################################################
# homology modeling of protein structures, and lots of it
#
# (C) 2020 Bryan Kolaczkowski, University of Florida, Gainesville, FL USA
# Released under GNU General Public License (GPL)
# bryank@ufl.edu
################################################################################

import sys
import os
import glob
import tempfile
import multiprocessing
import functools
import modeller
import modeller.automodel

import proteinRavel.align

################################################################################
# BEG VARIABLE DEFINITIONS

dnameprefix = 'structmodel_'  # prefix to use for temporary directories
fnameprefix = 'structmodel'   # prefix to use for temporary files

# END VARIABLE DEFINITIONS
################################################################################

################################################################################
# BEG CLASS DEFINITIONS

class StructmodelError(Exception):
    """Exceptions raised when structural modeling has an issue."""
    pass

# END CLASS DEFINITIONS
################################################################################

################################################################################
# BEG HELPER FUNCTION DEFINITIONS

def _build_models(structfname, basedir, nmodels, refstructure, verbose,
                  seq_rep_list):
    """
    Builds replicate structural models of a list of protein sequences.

    seq_rep_list is a list of (sequence,replicates) pairs, giving each
    sequence object to be modeled and the number of replicates needed for
    that sequence object

    SIDE EFFECT: models are placed in basedir/sequence_id directory
    """

    # set up path links, assuming current working directory
    workingdir  = os.getcwd()
    structfname = os.path.normpath(os.path.join(workingdir, structfname))
    basedir     = os.path.normpath(os.path.join(workingdir, basedir))

    for seq,reps in seq_rep_list:
        # check this sequence's existing structures; bail out if done
        mindex = 1
        outdir = basedir + os.path.sep + seq.identifier
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        else:
            existing_fnames = [ x.split(os.path.sep)[-1] for x in \
                                glob.glob(outdir + os.path.sep + 'rep*.pdb') ]
            existing_reps = [ int(x.split('rep')[1].split('.pdb')[0]) for \
                              x in existing_fnames]
            if existing_reps:
                existing_reps.sort(reverse=True)
                last_rep = existing_reps[0]
                if last_rep < reps:
                    mindex = existing_reps[0] + 1
                else:
                    continue

        # set up temporary directory for modeller execution
        with tempfile.TemporaryDirectory(prefix=dnameprefix) as tempdir:
            os.chdir(tempdir)

            # set up modeller environment
            if verbose:
                modeller.log.verbose()
            else:
                modeller.log.none()
            env = modeller.environ()
            env.io.atom_files_directory = [workingdir]

            # set up complete alignment
            aln = modeller.alignment(env)
            aln.append(file=structfname, remove_gaps=False)
            knowns = [s.code for s in aln]
            aln.append_sequence(seq.sequence)
            aln[-1].code = seq.identifier

            # write alignment - modeller doesn't like alignment in memory
            full_aln_fname = 'structaligntemp.ali'
            aln.write(full_aln_fname, alignment_format='PIR')

            # set up model assessments
            ASSESS_METHODS = [modeller.automodel.assess.DOPE,
                              modeller.automodel.assess.DOPEHR]
            ASSESS_NAMES   = ["DOPE score", "DOPE-HR score"]

            a = modeller.automodel.dope_loopmodel(env, alnfile=full_aln_fname,
                                                  knowns=knowns,
                                                  sequence=seq.identifier,
                                                  assess_methods=ASSESS_METHODS)
            a.starting_model = 1        # index of the first model
            a.ending_model   = nmodels  # index of the last model
            # adjust optimization parameters
            a.library_schedule = modeller.automodel.autosched.slow
            a.md_level         = modeller.automodel.refine.slow
            a.make()  # do homology modeling

            # evaluate structural models
            ok_models = [ x for x in a.outputs if x["failure"] is None ]
            score_results = []

            for data in ok_models:
                fname  = data["name"]
                myscrs = []
                for score_name in ASSESS_NAMES:
                    myscrs.append(data[score_name])
                ave_score = sum(myscrs) / len(myscrs)
                score_results.append((ave_score, fname, myscrs))

            score_results.sort()
            best_models = score_results[:reps]
            rest_models = score_results[reps:]

            # map to reference structure
            refseq = aln[0]
            if refstructure:
                refseq = aln[refstructure]

            refcode  = refseq.code
            refpdbf  = refseq.atom_file
            refrange = refseq.range
            refmdl   = modeller.model(env, file=refpdbf, model_segment=refrange)
            refpos   = modeller.selection(refmdl).only_atom_types('CA')

            # get best models
            final_files = []
            for (score,infname,scores) in best_models:
                outfname = outdir + os.path.sep + 'rep{}.pdb'.format(mindex)
                final_files.append(outfname)

                # build alignment
                myaln = modeller.alignment(env)
                myaln.append(file=structfname, align_codes=(refcode),
                             remove_gaps=False)
                myaln.append_sequence(seq.sequence)
                myaln[-1].code      = seq.identifier
                myaln[-1].atom_file = infname

                # read pdb file
                mymodel = modeller.model(env, file=infname)
                # translate to reference coordinates
                r = refpos.superpose(mymodel, myaln)
                # write translated pdb file
                mymodel.write(file=outfname)

                mindex += 1

            os.chdir(workingdir)

    return

# END HELPER FUNCTION DEFINITIONS
################################################################################

################################################################################
# BEG MODULE INTERFACE

def structmodel(ext_seqfname, anc_seqfname, struct_fname, smodel_dirname,
                replicates=5, nmodels=20, probcutoff=0.2, refstructure=None,
                threads=0, verbose=True):
    """
    Infers replicate structural models of extant and ancestral sequences.

    For each sequence in ext_seqfname and anc_seqfname, infers replicate
    structural models by homology modeling, based on structures in
    struct_fname. Note that ext_seqfname, anc_seqfname and struct_fname
    must already be aligned to one another.

    Parameters:
      ext_seqfname   (string): extant protein sequences in aligned FASTA format
      anc_seqfname   (string): ancestral sequences in ravel's .asr.seq format
      struct_fname   (string): structural templates in aligned PIR format
      smodel_dirname (string): directory to place inferred models
      replicates    (integer): number of final replicate models per sequence
      nmodels       (integer): number of initial models to build per sequence
      probcutoff      (float): don't include ancestral state with low probability
      refstructure   (string): map structure models to this structure
      threads       (integer): Number of threads to run for inference
      verbose       (boolean): Show function execution to stdout
    Returns:
      None
      SIDE EFFECT: structural models are in smodel_dirname, which must exist
    """

    # set threads to max available if not set manually
    if threads <= 0:
        threads = multiprocessing.cpu_count()

    # build models for extant sequences
    if verbose:
        sys.stdout.write('BEG Extant Protein Modeling\n')

    # set up base directory
    basedir = smodel_dirname + os.path.sep + 'extant'
    if not os.path.isdir(basedir):
        os.makedirs(basedir)

    # read alignment
    aln = proteinRavel.align.Alignment(ext_seqfname)

    # build replicate models for each extant sequence
    seq_rep_list = []
    for sequence in aln.sequencelist:
        seq_rep_list.append([[sequence,replicates]])

    with multiprocessing.Pool(threads) as threadpool:
        targetfn = functools.partial(_build_models, struct_fname, basedir,
                                                    nmodels, refstructure,
                                                    verbose)
        threadpool.map(targetfn, seq_rep_list)

    if verbose:
        sys.stdout.write('END Extant Protein Modeling\n')

    # build models for ancestral sequences
    if verbose:
        sys.stdout.write('BEG Ancestral Protein Modeling\n')

    # set up base directory
    basedir = smodel_dirname + os.path.sep + 'ancestral'
    if not os.path.isdir(basedir):
        os.makedirs(basedir)

    if verbose:
        sys.stdout.write('END Ancestral Protein Modeling\n')

    return

# END MODULE INTERFACE
################################################################################

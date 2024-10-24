"""
This file contains functions that assist in loading and pre-processin g the
pdb files.
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import Align
from urllib.request import urlopen
import os
import subprocess
import re
import warnings
import glob

# Internal
from . import load

def generate_reference(ref_pps, uniprot_ids):
    """
    Generate list of reference sequence objects. Either from uniprot_ids
    entries, or from the reference structure.
    """
    # Debug input
    assert len(ref_pps) == len(uniprot_ids),\
        "More/less chains than uniprot ids"    
    
    # Loop over uniprot chain ids
    reference_seqs = []
    for ind, up_id in enumerate(uniprot_ids):
        # Download and read fasta seq
        if up_id:
            download_fasta(up_id)
            ref_seq = SeqIO.read('data/seq/' + up_id + '.fa',
                                 format='fasta')
        # Or generate seq from ref_pps
        else:
            # Get pps in chain
            pp_seqs = [pp.get_sequence() for pp in ref_pps[ind]]
            seq_join = pp_seqs[0]
            # Stitch pps and generate seq
            for seq in pp_seqs[1:]:
                seq_join += seq
            # Generate reference sequence
            ref_seq = SeqIO.SeqRecord(seq_join,
                                      id='query',
                                      description='')
            # Fix sequences of unknown residues
            ref_seq = fix_unknown_sequence(ref_seq)
        reference_seqs.append(ref_seq)
    return reference_seqs

def fix_unknown_sequence(sequence):
    """
    If a sequence is made of unknwon residues, labeled X (6J5I, u), replace by
    a repetitive sequence of DTS residues, so that alignment is possible. This
    only occurs if in build_peptides we have aa_only=0.
    """

    if sequence.seq.count('X') == len(sequence):
        print(sequence)
        warnings.warn('Aligning pps made of unknown res.', stacklevel=2)
        s = 'DTS' * (len(sequence)//3) + 'DTS'[:len(sequence)%3]
        fixed_sequence = SeqRecord(id='query', description='', seq=Seq(s))
        return sequence
    else:
        return sequence


def download_fasta(up_id):
    """
    Download fasta sequence if it is not present
    """
    filename = 'data/seq/' + up_id + '.fa'
    if not os.path.exists(filename):
        target_url = 'http://www.uniprot.org/uniprot/' + up_id + '.fasta'
        url = urlopen(target_url)
        infile = open(filename, 'w')
        infile.write(url.read())
        infile.close()

        
def align_pair(structure1, structure2, all_atoms = False):
    """
    Receives a pair of structures and performs pairwise sequence 
    alignemnt for each of the chains. Returns a list with 2 dicts
    of labels.
    """
    
    n = len(structure1)
    temp1 = dict()
    temp2 = dict()
    
    for i in range(n):
        uni_ids = [None]
        my_ref_seqs = seq.generate_reference([structure1[i]], uni_ids)

        if all_atoms:
            aligned = seq.aligned_dict_all_atoms
        else:
            aligned = seq.aligned_dict

        label = aligned([structure1[i]], 
                        seq.align([structure1[i]], my_ref_seqs))
        temp1.update(label)

        label = aligned([structure2[i]], 
                        seq.align([structure2[i]], my_ref_seqs))
        temp2.update(label)
    
    labels = [temp1, temp2]
    
    return labels


def align_multiple(structures, all_atoms = False, ref_struc = 0):
    """
    Receives a list of structures and performs pairwise sequence 
    alignemnt for each of the chains using one of the structures
    as reference. Returns a lists with n = len(sturctures) dicts 
    of labels.
    """
    n = len(structures)
    
    labels = []
    for i in range(n):
        temp = align_pair(structures[ref_struc], 
                          structures[i], 
                          all_atoms)
        
        labels.append(temp[1])
    
    return labels


def align_monomers(structures, all_atoms = False, ref_struc = 0, ref_chain = 0):
    """
    Receives a list of structures from a homopolymer and performs 
    pairwise sequence alignemnt for each of the monomers using 
    one monomer as reference. Returns a list with n*m dicts of
    labels, with n = len(structures), m = len(structures[0]).
    """
    n = len(structures)
    m = len(structures[0])
    
    labels = []
    for i in range(n):
        for j in range(m):
            temp = align_pair([structures[ref_struc][ref_chain]], 
                              [structures[i][j]], 
                              all_atoms)
            
            labels.append(temp[1])
    
    return labels

def align_homopolymers(structures, all_atoms = False, ref_struc = 0, ref_chain = 0):
    """
    Receives a list of structures from a homopolymer and performs 
    pairwise sequence alignemnt for each of the structures using 
    one monomer as reference. Returns a lists with 
    n = len(sturctures) dicts of labels.
    """
    n = len(structures)
    m = len(structures[0])
    
    labels = []
    for i in range(n):
        temp = dict()
        for j in range(m):
            label0, label1 = align_pair([structures[ref_struc][ref_chain]], 
                                        [structures[i][j]], 
                                        all_atoms)
            
            temp.update(label1)
        labels.append(temp)
    
    return labels


def align(pps, ref_seqs):
    """
    Alignment the polypeptides in each chain of each protein to the reference
    sequences using local alignment. It returns the indices of the
    ref-seq and pps where the alignment of a segment starts and stops.
    """
    # Load sequence aligner
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score = -.5
    aligner.extend_gap_score = -.1
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0
    
    # For each chain&pep, start/stop indices of aligned segments from ref/pep
    start_stop = [] #(((ref_st1,ref_sp1),...,(N)), ((pep_st1,pep_sp1),...,(N))
   
    # Loop over proteins, chains and peptides
    for (chain, ref_seq) in zip(pps, ref_seqs):
        start_stop.append([])
        for pep in chain:
            # Align peptide ro reference sequence
            alignments = aligner.align(ref_seq.seq, pep.get_sequence())
            start_stop[-1].append(alignments[0].aligned)
            
    # Remove temp files
    for f in glob.glob("tmp/alignment_*"):
        os.remove(f)
    
    return start_stop


def align_score(pps, ref_seqs):
    """
    Alignment the polypeptides in each chain of each protein to the reference
    sequences using local alignment. It returns the indices of the
    ref-seq and pps where the alignment of a segment starts and stops.
    """
    # Load sequence aligner
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score = -.5
    aligner.extend_gap_score = -.1
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0
    score = 0

    # For each chain&pep, start/stop indices of aligned segments from ref/pep
    start_stop = [] #(((ref_st1,ref_sp1),...,(N)), ((pep_st1,pep_sp1),...,(N))
   
    # Loop over proteins, chains and peptides
    for (chain, ref_seq) in zip(pps, ref_seqs):
        start_stop.append([])
        for pep in chain:
            # Align peptide ro reference sequence
            alignments = aligner.align(ref_seq.seq, pep.get_sequence())
            score += alignments.score
            start_stop[-1].append(alignments[0].aligned)
            
    # Remove temp files
    for f in glob.glob("tmp/alignment_*"):
        os.remove(f)
    
    return start_stop, score



def aligned_dict(pps, start_stop):
    """
    A dictionary with keys the labels of the aligned residues in the pdb
    language, and values the labels in a self-made language that uses the
    reference sequence to identify each residue. This second language is
    common to all structures
    """
    
    aligned = {}
    
    # Loop over chains, peptides, residues and aligned segments
    for c_i, chain in enumerate(pps):
        for p_i, pep in enumerate(chain):
            for r_i, res in enumerate(pep):
                for s_i in range(len(start_stop[c_i][p_i][1])):
                    seg_ref_st = start_stop[c_i][p_i][0][s_i][0]
                    seg_ref_sp = start_stop[c_i][p_i][0][s_i][1]
                    seg_pep_st = start_stop[c_i][p_i][1][s_i][0]
                    seg_pep_sp = start_stop[c_i][p_i][1][s_i][1]
                    if r_i >= seg_pep_st and r_i < seg_pep_sp:
                        n_num = r_i + seg_ref_st - seg_pep_st
                        # Test duplicates
                        #assert res.full_id not in aligned.keys(),\
                        #'Duplicated key in aligned dict'
                        #assert (c_i, n_num) not in aligned.values(),\
                        #'Duplicated value in aligned dict'

                        # Store entry if residue passes test
                        if load.test_residue(res):
                            aligned[res.full_id] = (c_i, n_num)

    return aligned

def aligned_dict_all_atoms(pps, start_stop):
    """
    A dictionary with keys the labels of the aligned residues and atoms
    in the pdb language, and values the labels in a self-made language 
    that uses the reference sequence to identify each residue. This 
    second language is common to all structures
    """
    
    aligned = {}
    
    # Loop over chains, peptides, residues and aligned segments
    for c_i, chain in enumerate(pps):
        for p_i, pep in enumerate(chain):
            for r_i, res in enumerate(pep):
                for s_i in range(len(start_stop[c_i][p_i][1])):
                    seg_ref_st = start_stop[c_i][p_i][0][s_i][0]
                    seg_ref_sp = start_stop[c_i][p_i][0][s_i][1]
                    seg_pep_st = start_stop[c_i][p_i][1][s_i][0]
                    seg_pep_sp = start_stop[c_i][p_i][1][s_i][1]
                    if r_i >= seg_pep_st and r_i < seg_pep_sp:
                        n_num = r_i + seg_ref_st - seg_pep_st
                        
                        # Store entry if residue passes test
                        if load.test_residue(res):
                            for a_i, at in enumerate(res):
                                aligned[at.full_id] = (c_i, n_num, at.id)

    return aligned

def dict_from_label(label):
    """
    A dictionary with keys the labels of the aligned residues and atoms
    in the pdb language, and values the labels in a self-made language 
    that uses the reference sequence to identify each residue
    """
    
    output = {}
    keys = list(label.keys())
    values = list(label.values())
    
    for i in range(len(keys)):
        output[keys[i]] = (keys[i][3][1], keys[i][4][0], values[i][0])

    return output


def common(rel_dict, def_dict):
    """
    Compute the set of common residues by finding the intersect of two
    dictionaries. We use the self-made language (i.e., dict values)
    """
    # Create sets with residues labels
    rel_set = set(rel_dict.values())
    def_set = set(def_dict.values())

    # Calculate the intersect set
    common_res = rel_set.intersection(def_set)

    return common_res


def common_multiple(vec_dict):
    """
    Compute the set of common residues by finding the intersect of two 
    or more dictionaries. We use the self-made language (i.e., dict values)
    """
    # Create sets with residues labels
    vec_set = []
    for i in range(len(vec_dict)):
        vec_set.append(set(vec_dict[i].values()))
        
    common_res = vec_set[0].intersection(vec_set[1])
    for i in range(1, len(vec_dict)-1):
        temp_set = vec_set[i].intersection(vec_set[i+1])
        common_res = common_res.intersection(temp_set)
    
    return common_res

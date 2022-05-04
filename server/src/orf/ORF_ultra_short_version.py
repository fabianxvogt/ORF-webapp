from re import finditer, compile
from orf.fasta_utils import fasta_to_dict, write_fasta_record
import os
dir_path = os.path.dirname(os.path.realpath(__file__))


IN_DNA, OUT_DNA, OUT_AA, IUPAC_FNAME, LENGTH_FNAME = dir_path + "/input/in.fasta", dir_path + "/output/orfs_dna", dir_path + "/output/orfs_aa", dir_path + "/data/iupac.txt", dir_path + "/data/orf_length.txt" # FILENAMES
IUPAC = dict(x.split(" ") for x in open(IUPAC_FNAME).read().splitlines()) # Prepares IUPAC dictionary
def iupac_transl(orf): return ''.join([IUPAC[orf[i:i+3]] for i in range(0, len(orf)-2, 3)]) # Translates an ORF using UIPAC dict
def orf_regex(min, max): return compile(r"(?=(ATG([ACTG]{3}(?<!ATG|TAA|TAG|TGA)){"+str(min-1)+r","+str(max-1)+r"}(TAA|TAG|TGA)))") # ORF RegEx pattern with min and max length
def revcomp(dna): return dna[::-1].translate(str.maketrans("ATGC","TACG")) # Generates the reverse complement of a DNA string
def write_orf(dna_f, aa_f, orf, start, end): return [write_fasta_record(f, "ORF_"+str(start)+"_"+str(end), d) for f,d in ((dna_f,orf),(aa_f, iupac_transl(orf)))] # Saves an ORF
def print_output(out_files): [[print(l) for l in (f, open(f).read())] for f in [f for gen in out_files.values() for f in gen]] # Print the output files
def findORFs_from_file(infile, outfile_dna, outfile_aa, minlen, maxlen): # FindORFs Function
    inp_genoms = fasta_to_dict(infile) # Open input file
    out_files = {} # Variable to store output file names
    regex_pattern = orf_regex(minlen, maxlen) # Build the RegEx pattern
    for gen_id, dna in inp_genoms.items(): # Loop over all genoms 
        rev = revcomp(dna) # Get the reverse complement
        matches_dna, matches_rev = [tuple(finditer(regex_pattern, g)) for g in (dna, rev)] # Find all matches 
        out_files[gen_id] = [f + "_" + gen_id.strip()[0:10] + ".fasta" for f in (outfile_dna,outfile_aa)]
        dna_f, aa_f = [open(out_files[gen_id][i], "w") for i in (0,1)] # Open output files
        [[write_orf(dna_f, aa_f, g, s, e) for g,(s,e) in [(gen[i[0]:i[1]], ((len(gen)-i[0], len(gen)-i[1]+1) if r else (i[0]+1, i[1]))) for i in [m.regs[1] for m in ms]]] for ms, gen, r in [(matches_dna, dna, False),(matches_rev, rev, True)]] # Write ORFs -- I've taken it a bit too far on this line ... :-)
        
        [f.close() for f in (dna_f, aa_f)] # Close files
    return out_files # Success


def findORFs(dna, min, max):
    output_orfs = list()

    #inp_genoms = fasta_to_dict(infile) # Open input file
    #out_files = {} # Variable to store output file names
    regex_pattern = orf_regex(min, max) # Build the RegEx pattern
    rev = revcomp(dna) # Get the reverse complement
    matches_dna, matches_rev = [tuple(finditer(regex_pattern, g)) for g in (dna, rev)] # Find all matches 
    [[output_orfs.append({'START': s, 'END': e, 'ORF': g, 'IUPAC': iupac_transl(g)}) for g,(s,e) in [(gen[i[0]:i[1]], ((len(gen)-i[0], len(gen)-i[1]+1) if r else (i[0]+1, i[1]))) for i in [m.regs[1] for m in ms]]] for ms, gen, r in [(matches_dna, dna, False),(matches_rev, rev, True)]] # Write ORFs -- I've taken it a bit too far on this line ... :-)

    return sorted(output_orfs, key=lambda d: d['START'])  # Success
#findORFs(IN_DNA, OUT_DNA, OUT_AA, min, max)
#[(print("Testing ORF length "+str(l)+"-"+str(h)), print_output(findORFs(IN_DNA, OUT_DNA, OUT_AA, l, h))) for l,h in ((4,10),(1,10))] # Test with lengths 1-10 and 4-10 
from io import TextIOWrapper

# Reads FASTA file from path into dictionary
# Params:
# - fasta_file_path: FASTA file path to read from
def fasta_to_dict(fasta_file_path: str):
    d = {}
    k = ""
    v = []
    with open(fasta_file_path) as f:
        for l in f.read().splitlines():
            if l.startswith(">"):
                # Merge previous data string
                if (k != ""):
                    d[k] = "".join(v)
                # Set new key
                k = l[1:]
                # Init dict
                d[k] = ""
            else: 
                v += [l]
        # Merge last data string
        if (k != ""):
            d[k] = "".join(v)
    return d

# Write a record to a file in FASTA format
# Params:
# - FASTA_file: FASTA file to write to
# - descr: FASTA description line
# - data: FASTA data line 
def write_fasta_record(fasta_file: TextIOWrapper, descr: str, data: str): 
    [fasta_file.write(x) for x in (">" + descr + "\n", data + "\n")]
    return fasta_file
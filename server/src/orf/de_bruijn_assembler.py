import sys

# Read-Class: For FASTA file reads
# DBGnode-Class: Single node in De Bruijn graph
# DBGraph-Class: De Bruijn graph implementation

class Read:

    lines = []
    name = ""
    bases = ""

    def __init__(self, lines):
        if (lines[0][0] == '>'):
            self.name = lines[0][1:]
            lines = lines[1:]
        self.lines = lines
        self.bases = ''.join(lines)

    def get_kmers(self, kmersize):
        result = {}
        seq = ''.join(self.lines)
        for i in range(0,len(seq)-kmersize+1):
            if seq[i:i+kmersize] in result:
                result[seq[i:i+kmersize]] += 1
            else:
                result[seq[i:i+kmersize]] = 1
        return result

    def __str__(self):
        if (len(self.bases) > 20):
            bases = self.bases[0:20] + "..."
        else:
            bases = self.bases
        return self.name + ": " + bases

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        if (type(self) != type(other)):
            return False

        if ((self.name == other.name) and
           (self.bases == other.bases)):
           return True
        return False

class DBGnode:
    seq = ""
    alphabet = ['A', 'G', 'T', 'C']
    edges_to = {}
    edges_from = {}

    def __init__(self, seq):
        self.edges_from = {}
        self.edges_to = {}
        self.seq = seq

    def add_edge_to(self, eto):
        if eto.seq in self.edges_to:
            self.edges_to[eto.seq] += 1
        else:
            self.edges_to[eto.seq] = 1

    def add_edge_from(self, efrom):
        if efrom.seq in self.edges_from:
            self.edges_from[efrom.seq] += 1
        else:
            self.edges_from[efrom.seq] = 1

    def get_potential_from(self):
        result = []
        for char in self.alphabet:
            result.append(char + self.seq[:len(self.seq)-1])
        return result

    def get_potential_to(self):
        result = []
        for char in self.alphabet:
            result.append(self.seq[1:] + char)
        return result

    def get_edge_to_weight(self, other):
        if other.seq not in self.edges_to:
            return 0
        return self.edges_to[other.seq]

    def get_edge_from_weight(self, other):
        if other.seq not in self.edges_from:
            return 0
        return self.edges_from[other.seq]

    def count_edges(self):
        return len(self.edges_from) + len(self.edges_to)

    def __str__(self):
        return self.seq

    def __eq__(self, other):
        if (type(self) != type(other)):
            return False

        if (self.seq == other.seq):
           return True
        return False

class DBGraph:
    nodes = {}
    kmer_len = -1
    def __init__(self):
        self.nodes = {}
        return

    def add_kmers(self, kmers):
        if (self.kmer_len == -1):
            self.kmer_len = len(list(kmers.items())[0][0])

        for i in range(0, len(kmers)):
            seq = list(kmers)[i]
            if (len(seq) != self.kmer_len):
                raise ValueError('K-mer lengths not equal') 
            
            # check if node exists already
            if not seq in self.nodes:
                # add the node
                self.nodes[seq] = DBGnode(seq) 

            # create edges
            for from_seq in self.nodes[seq].get_potential_from():
                if from_seq in self.nodes:
                    self.nodes[from_seq].add_edge_to(self.nodes[seq])
                    self.nodes[seq].add_edge_from(self.nodes[from_seq])
            for to_seq in self.nodes[seq].get_potential_to():
                if to_seq in self.nodes:
                    self.nodes[to_seq].add_edge_from(self.nodes[seq])
                    self.nodes[seq].add_edge_to(self.nodes[to_seq])

    def count_edges(self):
        sum = 0
        for node in self.nodes.values():
            sum += node.count_edges()
        return sum/2

    def count_nodes(self):
        length = len(self.nodes)
        return length

    def __str__(self):
        pass

def read_fasta(readfile):
    reads = []
    start = 0
    end = 0
    
    with open(readfile) as f:
        fasta_lines = f.read().splitlines()
        for i in range(0, len(fasta_lines)):
            if fasta_lines[i][0] == '>':
                if i == 0:
                    continue
                else:
                    end = i
                    reads.append(Read(fasta_lines[start: end]))
                    start = i
        reads.append(Read(fasta_lines[start:]))

    return reads

def build_graph(filename, kmersize):
    reads = read_fasta(filename)
    graph = DBGraph()

    for read in reads:
        graph.add_kmers(read.get_kmers(kmersize))
    
    return graph
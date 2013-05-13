import random
from Bio import SeqIO
from Bio import Seq
from Bio.Alphabet import generic_dna

dnaindex = {'A':0, 'C':1, 'G':2, 'T':3, 'a':0, 'c':1, 'g':2, 't':3}
nucleotide = {0:'A', 1:'C', 2:'G', 3:'T'}
class SequenceStats:
    def __init__(self, A,G,C,T,X,avglen,count):
        self.A = A
        self.G = G
        self.C = C
        self.T = T
        self.X = X
        self.avglen = avglen
        self.count = count
        self.Accumulative = [A, A+C, A+C+G, 1.0]


def _Nucleotide(N):
    if N==0:
        return 'A'
    elif N==1:
        return 'C'
    elif N==2:
        return 'G'
    elif N==3:
        return 'T'
    else:
        return 'N'

def _index(N):
    global dnaindex

    if N in dnaindex:
        return dnaindex[N]
    else:
        return N


def ReverseComplement(seq_str):
    seq = Seq.Seq(seq_str, generic_dna)
    return str(seq.reverse_complement())


def _BgEColi(N):

    #regulonDB: these numbers come from RegulonDB.  more similar to intergenic    
    bg_R = [ 0.29088, 0.20805, 0.20465, 0.29641]

    #myCalculations:
    bg_G = [0.2461870712927091, 0.25423203133840194, 0.2536649657572998, 0.24591593161158917]

    bg_I = [0.2820738035859304, 0.21793195601965176, 0.22167713954948365, 0.27831710084493416]

    bg_C = [0.2411798593237505, 0.24509166256283507, 0.27318020947390914, 0.24054826863950526]

    return bg_I[_index(N)]




def createARandomSeq(stats):
    seq = ''

    for l in xrange(random.randint(stats.avglen * 7 / 10 , stats.avglen * 13 / 10)):
        r = random.random()
        i = 0
        while stats.Accumulative[i]<=r:
            i += 1
        n = _Nucleotide(i)
        seq += n
    return seq



def creatRandomSeq(outputfilename, SEQ_COUNT, SEQ_LENGTH, dist):
    file = open(outputfilename, 'w')
    random.seed()

    for i in [1,2]:
        dist[i] = (dist[i-1] + _BgEColi(i))
    

    for s in xrange(SEQ_COUNT):
        file.write(">randomseq" + str(s) + "\n")

        s = ''
        for l in xrange(random.randint(SEQ_LENGTH * 7 / 10 , SEQ_LENGTH * 13 / 10)):
            r = random.random()
            i = 0
            while dist[i]<=r:
                i += 1
            n = _Nucleotide(i)
            s += n

        file.write(s + '\n')
    file.close()




def calculateDist(fastafile):
    input_file = SeqIO.parse(fastafile, "fasta")

    cnt = avglen = A = C = G = T = X = 0
        
    while True:
        try:
            entry = input_file.next()
            
        except Exception as ex:
            print ex
            break

        s = entry.seq.tostring()
        for N in s:
            if N=='A' or N=='a':
                A += 1
            elif N=='C' or N=='c':
                C += 1
            elif N=='G' or N=='g':
                G += 1
            elif N=='T' or N=='t':
                T += 1
            else:
                X += 1
        
        avglen += len(s)
        cnt += 1

    input_file.close()
    N = 1.0 *( A + C + G + T + X)
    return SequenceStats(A/N, C/N, G/N, T/N, X/N, avglen/cnt,cnt)
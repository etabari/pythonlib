#!/usr/bin/python
import sys, os
from Bio import SeqIO

from dna import ReverseComplement
### [ Gene, ...]
###
###

def parseValue(x):
    if x == '' or x == 'NULL':
        return None
    try:
        y = x
        y = float(x)
        y = int(x)
    finally:
        return y



class Gene:
    def __init__(self):
        self.posStart = -1
        self.posEnd = -1
        self.strand = None
        self.aaLength = -1 
        self.gi = -1
        self.name = None
        self.bnumber = -1
        self.code = None
        self.cog = None
        self.product = None
        self.chr = None

    def __str__(self):
        return str(self.gi)+'('+self.strand+')'

    def __cmp__(self, other):
        if self.strand == other.strand:
            return self.posStart - other.posStart
        elif self.strand =='+':
            return sys.maxint
        else:
            return -sys.maxint

    def __eq__(self, other):
        return self.gi == other.gi

    def __hash__(self):
        return hash(self.gi)

     
    def intergenic(self, other):
        if self.strand == other.strand and self.chr == other.chr:
            return self.chr.sequence[self.posEnd:other.posStart-1]
        else:
            raise Exception('Erorr')
    
    def intergenicDistance(self, other):
        if self.strand == other.strand and self.chr == other.chr:
            intdis = other.posStart - self.posEnd
            return intdis
        else:
            raise Exception('Erorr')


    ### if gene is not found (!) it will return (None,None) 
    ### if the gene is the first gene ===> (-1, None)
    ### if it is the last gene ===> (genecount, None) 
    def previousRegardlessOfStrand(self):
        prevGene = None
        prevIndex = None
        for i in xrange(len(self.chr.genes_regardless_of_strand)):
            if self.chr.genes_regardless_of_strand[i].gi == self.gi:
                #print i
                selfIndex = i
                if self.strand == '+':
                    prevIndex = i - 1
                else:
                    prevIndex = i + 1
                
                if prevIndex>=0 and prevIndex<len(self.chr.genes_regardless_of_strand):
                    prevGene = self.chr.genes_regardless_of_strand[prevIndex]

                break
        #print prevIndex
        return (prevIndex, prevGene)
    
    def previous(self):
        prevGene = None
        prevIndex = None
        for i in xrange(len(self.chr.genes)):
            if self.chr.genes[i].gi == self.gi:
                #print i
                selfIndex = i
                if self.strand == '+':
                    prevIndex = i - 1
                else:
                    prevIndex = i + 1
                
                if prevIndex>=0 and prevIndex<len(self.chr.genes): 
                    prevGene = self.chr.genes[prevIndex]
                    if prevGene.strand != self.strand:
                        prevGene = None
                break
        #print prevIndex
        return (prevIndex, prevGene)
    

    def upstreamCoordinate(self, maxLength):

        (prevIndex, prevGene) = self.previous()
        
        if not prevIndex:
            return (-1,-1)
        
        if self.strand == '+':
            if prevGene:
                st = max( 0 , self.posStart - min(maxLength, self.posStart - prevGene.posEnd) )
            else:
                st = max(0, self.posStart - maxLength)
            return  (st,self.posStart-1)
        else:
            if prevGene:
                en =  min( len(self.chr.sequence) , self.posEnd + min(maxLength, prevGene.posStart - self.posEnd) )
            else:
                en =  min( len(self.chr.sequence) , self.posEnd + maxLength )
            #print self, prevGene

            return (self.posEnd,en-1)


    def upstream(self, maxLength):

        (start, end) = self.upstreamCoordinate(maxLength)
        
        if start<0:
            return ''
        
        if self.strand == '+':
            return  self.chr.sequence[start:end]
        else:
            return ReverseComplement(self.chr.sequence[start:end])
            
            




class Chromosome:
    def __init__(self):
        self.seqFile = None
        self.sequence = None
        self.length = 0
        self.name = None
        self.genes = []
        self.genes_regardless_of_strand = []

        
    def importData(self, folder, name):
        self.name = name

        self._importPtt(os.path.join(folder, name +'.ptt'))

        self.seqFile = SeqIO.parse(os.path.join(folder, name + '.fna'), "fasta")
        self.sequence = str(self.seqFile.next().seq)
        self.length = len(self.sequence)

        

    def _importPtt(self, filename):

        #print "Importing Ptt ..."
        
        colCount = 9 # there should be 9 columns in the file

        genes = []
        
        with open(filename) as input_file:
            line = input_file.readline()
            line_count = 1
            insert_count = 0
            while line:

                line = line.strip()
                values = line.split("\t")

                if line=='' or line[0] == '#' or not '..' in values[0] or len(values)<4:
                    line = input_file.readline()
                    continue
            
                while len(values)<9:
                    values += ['']

                g = Gene()
                #start and end
                se = values[0].split('..')
                
                g.posStart = parseValue(se[0])
                g.posEnd = parseValue(se[1])

                g.strand = parseValue(values[1])
                g.aaLength = parseValue(values[2])
                g.gi = parseValue(values[3])
                g.name = parseValue(values[4])
                g.bnumber = parseValue(values[5])
                g.code = parseValue(values[6])
                g.cog = parseValue(values[7])
                g.product = parseValue(values[8])
                g.chr = self
           
                #posstart, posend, strand,aa_length, gi, name, bnumber, code, cog, product
                genes += [g]

                line = input_file.readline()

        self.genes = sorted(genes)
        self.genes_regardless_of_strand = sorted(genes, lambda g1, g2: cmp(g1.posStart,g2.posStart)) 
        


class Genome:
    def __init__(self):
        self.path = None
        self.name = None
        self.geneCount = 0
        self.size = 0
        self.chromosomes = []
        self.operons = {}
    
    #def _readOperonFile(self, fileName):
    #    f = open(fileName)
    #    first = True
    #    for l in f:
    #        if first:
    #            First = False
    #            continue
    #        col = l.split(',')
    #        self.operons[(int(col[0]), int(col[1]))] = col[3]
    #    f.close()

    def readFolder(self, folder):
        self.path = folder
        for fileName in os.listdir(folder):
            (fileBaseName, fileExtension) = os.path.splitext(fileName)
            if fileExtension == '.fna':
                chr = Chromosome()
                chr.importData(folder, fileBaseName)
                self.chromosomes += [chr]
                self.geneCount += len(chr.genes)
                self.size += len(chr.sequence)
            if fileExtension == '.aa':
                self.name = fileBaseName
            #if fileExtension =='.operons':
            #    self._readOperonFile(os.path.join(folder, fileName))


    def getGenePairs(self):
        gene_pairs = []
        for chr in self.chromosomes:
            for i in xrange(len(chr.genes)-1):
                if chr.genes[i].strand == chr.genes[i+1].strand:
                    gene_pair = (chr.genes[i], chr.genes[i+1])
                    gene_pairs += [gene_pair]

        return gene_pairs

     
    def getGeneFromGi(self, gi):
        for chr in self.chromosomes:
            for gene in chr.genes:
                if gene.gi == gi:
                    return gene
        return None

    def getGeneFrom(self, attr, value):
        v = value.lower()
        for chr in self.chromosomes:
            for gene in chr.genes:
                if getattr(gene,attr).lower() == v:
                    return gene
        return None

    def getGeneFromGname(self, gname):
        return self.getGeneFrom('name', gname)


    def distanceIgnoringDirectionality(self, gi1, gi2):
        for chr in self.chromosomes:
            for i in xrange(len(chr.genes_regardless_of_strand)):
                if chr.genes_regardless_of_strand[i].gi == gi1:
                    
                    for j in xrange(len(chr.genes_regardless_of_strand)):
                        if chr.genes_regardless_of_strand[j].gi == gi2:
                            if gi1 == 16077929:
                                for k in xrange(i-1,j+1):
                                    print chr.genes_regardless_of_strand[k]

                            return abs(j-i)

        return self.geneCount        
                          
                               
    def distance(self, gi1, gi2):
        g1chr = None
        g1gene = None
        g1index = None
        for chr in self.chromosomes:
            for i in xrange(len(chr.genes)):
                if chr.genes[i].gi == gi1:
                    g1chr = chr
                    g1gene = chr.genes[i]
                    g1index = i
                    break

        for i in xrange(len(g1chr.genes)):
            if g1chr.genes[i].gi == gi2:
                if g1gene.strand == chr.genes[i].strand:
                    return abs(i-g1index)
                else:
                    return len(g1chr.genes)

        return self.geneCount        
            


def ReadOperonPairs(file_name):
    operons = {}
    operon_file = open(file_name)
    for line in operon_file:
         values = line.strip().split(',')
         try:
            gi1 = int(values[0])
            gi2 = int(values[1])
            opName = values[3]
            operons[(gi1, gi2)] = opName
         except:
            pass
    return operons



if __name__ == '__main__':

    #External initialization
    #import os,sys
    #os.chdir('c:/Users/Ehsan/Work/Lab/operon/python')
    #sys.path += ['./']
    #import genome
    #gn = genome.Genome()
    #


    print 'TESTING'

    
    gn = Genome()
    gn.readFolder('../genomes/Acetobacter_pasteurianus_IFO_3283_01_uid59279')


    print len(gn.chromosomes)
    gp = gn.getGenePairs()

 
    for (g1,g2) in gp:
        if g1.strand == '-':
            print g1, g2
            break

    raw_input()

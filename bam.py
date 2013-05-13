from itertools import imap
import pysam
import dna



####reads_GC



def GenomicReads(bamfile, bam_chrname, genome, extend):
	for chr in genome.chromosomes:
		for gene in chr.genes:
			start = max(0, gene.posStart - extend)
			end = min(chr.length, gene.posEnd + extend)
			is_reverse = gene.strand=='-'
			for read in bamfile.fetch(bam_chrname, start, end):
				if read.pos >= start and read.pos + read.rlen < end and read.is_reverse == is_reverse:
					read.gene = gene
					yield read


def GenomicReadsWithGene(bamfile, bam_chrname, genome, extend):
	for chr in genome.chromosomes:
		for gene in chr.genes:
			start = max(0, gene.posStart - extend)
			end = min(chr.length, gene.posEnd + extend)
			is_reverse = gene.strand=='-'
			for read in bamfile.fetch(bam_chrname, start, end):
				if read.pos >= start and read.pos + read.rlen < end and read.is_reverse == is_reverse:
					yield (read, gene)

#for read in GenomicReads(bamfile, 'gi|49175990|ref|NC_000913.2|', g, 20):
#	break





def _applyWindow(ll, w):
	rr = []
	for i in xrange(len(ll)):
		s = max(0, i - w/2)
		e = min(i + w/2, len(ll))
		rr += [1.0 * sum(ll[s:e]) / (e-s)]


	return rr
	

def _normalize(ll):
	if ll and len(ll)>0:
		min_v = min(ll)
		size_v = max(ll) - min_v * 1.0
		if size_v > 0:
			return [(x-min_v)/size_v for x in ll]

	return [0]*len(ll)

#
# BOWTIE reverese complements the reads to match the leading strand.
# So, if the read is on reverse strand, start of read is read.pos + read.rlen
#
# from bowtie refrence:
#    POS: 0-based offset into the forward reference strand where leftmost character of the alignment occurs
#
#
def countExprReadStart(sam_file, chromosome_name, Start, End, Strand, windowSize,  normalize):
	result = [0]*(End-Start)
	is_reverse = Strand=='-'
	for read in sam_file.fetch(reference=chromosome_name, start=Start, end=End):
		if  is_reverse:
			read_pos = read.pos + read.rlen
		else:
			read_pos = read.pos

		if read_pos >= Start and read_pos < End and is_reverse==read.is_reverse:
			result[read_pos-Start] += 1
	#
	if windowSize>1:
		result = _applyWindow(result, windowSize)
	if normalize:
		result = _normalize(result)
	return result

def countExprReadPileUp(sam_file, chromosome_name, Start, End, Strand, windowSize,  normalize):
	result = [0]*(End-Start)
	is_reverse = Strand=='-'
	for pileupcolumn in sam_file.pileup(reference=chromosome_name , start=Start, end=End):
			if pileupcolumn.pos >= Start and pileupcolumn.pos < End and is_reverse==read.is_reverse:
				result[pileupcolumn.pos-Start] = pileupcolumn.n

	
	#return result
	if windowSize>1:
		result = _applyWindow(result, windowSize)

	if normalize :
		result = _normalize(result)
	return result


####################


def countExprKmer(sam_file, chromosome_name, start, end, kmerBias, windowSize,  normalize):
	result = [0]*(end-start)

	kmer_size = len(kmerBias.keys()[0])

	pileups = get_pileups(sam_file, chromosome_name, start, end)
	for i in xrange(end-start):
		if i in pileups:
			for pu in pileups[i]:
				kmer = pu.alignment.seq[0:kmer_size]
				if kmer in kmerBias:
					kmer_weight = kmerBias[kmer]
					result[i] += kmer_weight

	#return result
	if windowSize>1:
		result = _applyWindow(result, windowSize)
	
	if normalize:
		result = _normalize(result)
	return result



def countExprKmerReadStart(sam_file, chromosome_name, Start, End, Strand, kmerBias, windowSize,  normalize):
	
	result = [0]*(End-Start)
	kmer_size = len(kmerBias.keys()[0])
	is_reverse = Strand=='-'

	for read in sam_file.fetch(reference=chromosome_name, start=Start, end=End):
		if is_reverse:
			read_pos = read.pos + read.rlen
			kmer = read.seq[-kmer_size:]
			kmer = dna.ReverseComplement(kmer)
		else:
			read_pos = read.pos
			kmer = read.seq[:kmer_size]

		if read_pos >= Start and read_pos < End and is_reverse==read.is_reverse and kmer in kmerBias:
			kmer_weight = kmerBias[kmer]
			result[read_pos-Start] += kmer_weight
	#
	if windowSize>1:
		result = _applyWindow(result, windowSize)
	if normalize:
		result = _normalize(result)
	return result



########
"""
from itertools import imap
import pysam
bamfile = pysam.Samfile('Cluster_6_trimmedAdapter.sorted.bam')
chrname = 'gi|49175990|ref|NC_000913.2|'
col = None
for pileupcolumn in bamfile.pileup(chrname , 1000, 1200):
			if pileupcolumn.pos >= 1020 and pileupcolumn.pos < 1200:
				col = pileupcolumn
				break
pileup = None
for pu in col.pileups:
	pileup = pu
	break


col = get_pileups(bamfile, chrname, 1200)

pileups = get_pileups(bamfile, chrname, 1000, 1200)


"""

def get_pileups(bam, chrom, pos):
     for pileupcol in bam.pileup(chrom, pos, pos+100):
         if pileupcol.pos == pos:
             return pileupcol.pileups

def get_pileups(bam, chrom, start, end):
    result = {}
    for pileupcol in bam.pileup(chrom, start, end):
        if pileupcol.pos >= start and pileupcol.pos < end:
            result[pileupcol.pos-start] = pileupcol.pileups
    return result





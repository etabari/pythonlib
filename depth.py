from itertools import imap

class Peak:
	def __init__(self,start,end,size, change):
		self.start = start
		self.end = end
		self.size = size
		self.change = change
	#
	#
	def __cmp__(self, other):
		return self.start.__cmp__(other.start)
	#
	#
	def __str__(self):
		change = '*'
		if self.change>0:
			change = '+'
		elif self.change<0:
			change = '-'
		return '['+str(self.start) +','+str(self.end)+'] @' + change +str(self.size)


class GeneExpression:
	def __init__(self):
		self.name = None
		self.data = None
	#
	def setFromDepthLine(self, line):
		data = line.strip().split('\t')
		self.name = data[0]
		self._setGeneInfo()
		#	
		self.data = []
		for d in data[1:]:
			self.data += [float(d)]
		self.max_val = max(self.data) - min(self.data)
		self.min_val = min(self.data)
	#
	#
	def setFromNameData(self, gene_str, expr):
		self.name = gene_str
		self._setGeneInfo()
		self.data = expr
		self.max_val = max(self.data) - min(self.data)
		self.min_val = min(self.data)
	#
	#
	def __eq__(self, other):
		return self.gi == other.gi
	#
	def __hash__(self):
		return hash(self.gi)
#
	def __cmp__(self, other):
		return cmp(self.name, other.name)
	#
	#
	def _setGeneInfo(self):
		gene = self.name.split(',')
		self.gName = gene[0][1:]
		self.strand = gene[2]
		self.start = int(gene[3])
		self.end = int(gene[4])
		self.length = self.end - self.start 
		self.gi = int(gene[5])
	#
	#
	def peakRegion(self, start=0, end=-1):
		if end == -1:
			end =len(self.data)
		max_val = max(self.data[start:end])
		max_pos = -1
		for i in xrange(start,end):
			if self.data[i] == max_val:
				max_pos = i
				break
		peak_s = max_pos - 1
		while peak_s>=start and self.data[peak_s]<=self.data[peak_s+1]:
			peak_s -= 1
		#
		peak_e = max_pos + 1
		while peak_e<end and self.data[peak_e]<=self.data[peak_e-1]:
			peak_e += 1
		return (peak_s + 1, peak_e)
	#
	#
	def peaks(self, start=0, end=-1):
		if end == -1:
			end =len(self.data)
		s,e = self.peakRegion(start,end)
		#
		peakMax = max(self.data[s:e])
		peakSv = self.data[s]
		peakEv = self.data[e-1]
		peak = Peak(s, e, peakMax-min(peakSv,peakEv),  peakEv - peakSv)
		#
		print '(',start,':',end,')--->','[',s,':',e,']', e-s, '@', peakSv, peakMax, peakEv, 'C', peakEv - peakSv
	 	if start == s and e == end:
	 		return [peak]
	 	else:
	 		right = left = []
			if s > start:	 		
	 			right = self.peaks(start,s)
	 		if e<end:
	 			left = self.peaks(e,end)
	 		return right + left + [peak]
	#
	#
	def cliffRegion(self, start=0, end=-1, ratio=.8):
		if end == -1:
			end =len(self.data)
		max_val = max(self.data[start:end])
		max_pos = -1
		for i in xrange(start,end):
			if self.data[i] == max_val:
				max_pos = i
				break
		peak_s = max_pos - 1
		while peak_s>=start and self.data[peak_s] >  ratio * self.max_val + self.min_val:
			peak_s -= 1
		#
		peak_e = max_pos + 1
		while peak_e<end and self.data[peak_e] >  ratio * self.max_val + self.min_val:
			peak_e += 1
		return (peak_s + 1, peak_e)
	#
	#
	#
	def cliffs(self, start=0, end=-1, ratio=.8):
		try:
			if end == -1:
				end =len(self.data)
			s,e = self.cliffRegion(start, end, ratio)
			#
			cliffMax = max(self.data[s:e])
			#print '(',start,':',end,')--->','[',s,':',e,']',  '@', self.data[s], cliffMax, self.data[e-1], '|', .5*max(self.data), .8*max(self.data)
			result = []
			if cliffMax > 0 and cliffMax > ratio * self.max_val + self.min_val:
				result = [Peak(s, e, cliffMax, 0)] 
			#
			#
		 	if (start == s and e == end) or cliffMax==0:
		 		return result
		 	else:
		 		right = left = []
				if s > start:	 		
		 			right = self.cliffs(start, s, ratio)
		 		if e<end:
		 			left = self.cliffs(e, end, ratio)
		 		return right + left + result
		except:
		 	print 'ERROR ', self.name

	#
	#
	#
	def plateauRegion(self, start=0, end=-1, ratio=.1):
		if end == -1:
			end =len(self.data)
		#
		min_val = min(self.data[start:end])
		min_pos = -1
		for i in xrange(start,end):
			if self.data[i] == min_val:
				min_pos = i
				break
		peak_s = min_pos - 1
		while peak_s>=start and self.data[peak_s] < ratio * self.max_val + self.min_val:
			peak_s -= 1
		#
		peak_e = min_pos + 1
		while peak_e<end and self.data[peak_e] < ratio * self.max_val + self.min_val:
			peak_e += 1
		return (peak_s + 1, peak_e)
	#
	#
	def plateaus(self, start=0, end=-1, ratio=.1):
		if end == -1:
			end =len(self.data)
		s,e = self.plateauRegion(start, end, ratio)
		#
		plateauMax = max(self.data[s:e])
		result = []
		if plateauMax > 0 and plateauMax < ratio * self.max_val + self.min_val:
			result = [Peak(s, e, plateauMax, 0)] 
		#
		#
	 	if (start == s and e == end) or plateauMax==0:
	 		return result
	 	else:
	 		right = left = []
			if s > start:	 		
	 			right = self.plateaus(start, s, ratio)
	 		if e<end:
	 			left = self.plateaus(e, end, ratio)
	 		return right + left + result

	#
	#



def readFile(file_name):
	genes = {}
	with open(file_name) as f:
		line = f.readline()
		while line:
			ge = GeneExpression()
			ge.setFromDepthLine(line)
			genes[ge.name] = ge
			line = f.readline()
	return genes


def writeFile(file_name, genes):
	with open(file_name, 'w') as f:
		for gene_name in genes:

			f.write(gene_name +'\t')
			for x in genes[gene_name].data:
				f.write('%.2f\t' % x)
			f.write('\n')


######################################
#
#
#

def gcCount(seq):
	return (seq.count('C') + seq.count('G') + seq.count('S') + seq.count('c') + seq.count('g') + seq.count('s') ) 



def gcPercent(big_seq, pos, window):
	start = pos - window/2
	end = pos + window/2
	if start<0:
		start = 0
	if end>=len(big_seq):
		end = len(big_seq)-1

	return 100.0*gcCount(big_seq[start:end])/(end-start)
	


def countGC(big_seq, w, start, end):
	rr = []
	for i in xrange(end-start):
		s = start + i - w/2
		e = start + i + w/2
		if s<0:
			s=0
		if e>=len(big_seq):
			e = len(big_seq)
		rr += [100.0*gcCount(big_seq[s:e])/(e-s)]
	return rr

def geneToStr(gene):
	return '>{0},1,{1},{2},{3},{4}'.format(gene.name, gene.strand, gene.posStart, gene.posEnd, gene.gi)



def pearsonr(x, y):
	n = min(len(x), len(y))

	sum_x = float(sum(x[:n]))
	sum_y = float(sum(y[:n]))
	
	sum_x_sq = sum(map(lambda ex: ex * ex, x[:n]))
	sum_y_sq = sum(map(lambda ey: ey * ey, y[:n]))
	psum = sum(imap(lambda ex, ey: ex * ey, x[:n], y[:n]))

	num = psum - (sum_x * sum_y/n)
	den = pow((sum_x_sq - pow(sum_x, 2) / n) * (sum_y_sq - pow(sum_y, 2) / n), 0.5)

	if den == 0: 
		return 0
	return num / den

##############
#
# genes = readFile('ecoli-57779.genes.cluster6.w10.depth')
# pnhD = genes['>phnD,1,-,4321359,4322375,16131931']
# genes = readFile('ecoli-57779.genes.cluster4.w10.depth')
# gnError = genes['>wcaJ,1,-,2118184,2119578,16129987']
# pnhD.peakRegion()
# ps = sorted(pnhD.peaks())
# for p in ps:
# 	print p

# genes = readFile('ecoli-57779.genes.cluster6.w10.depth')
# pnhD = genes['>phnD,1,-,4321359,4322375,16131931']
# for p in sorted(pnhD.cliffs()):
# 	print p
# for p in sorted(pnhD.plateaus()):
# 	print p

def makeRcommands(phnD):
	print 'data <- c(', phnD.data[0],
	for x in pnhD.data[1:]:
		print ',',x,
	print ')'
	#
	ps = sorted(pnhD.peakCount())
	#
	color='red'
	print 'plot(data, col="white")'
	for p in ps:
		if color == 'red':
			color = 'blue'
		else:
			color = 'red'
		print 'points( c(', p.start,':',p.end, '),data[',p.start,':',p.end,'], col="'+color+'")'


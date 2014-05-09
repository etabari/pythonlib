

class ProgressBar:
	def __init__(self, maxProgress, onProgressPercent=10, printSymbol='.'):
		self.max = maxProgress
		self.current = 0
		self.ratio = 100 / onProgressPercent
		self.symbol = printSymbol
		self.last_printed_precent = 0
		print '('+str(maxProgress)+')',


	def step(self):
		if self.current<self.max:
			self.current += 1

	def stepPrint(self):

		if self.current<self.max:
			self.current += 1

		new =  self.current * self.ratio / self.max 

		#print self.current
		if new > self.last_printed_precent:
			self.last_printed_precent = new
			print self.symbol,






if __name__=='__main__':
	print 'testing progress bar '
	p = ProgressBar(500)
	for x in xrange(500):
		p.stepPrint()
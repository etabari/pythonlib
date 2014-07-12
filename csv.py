#!/bin/usr/python


def parseValue(strV):
	if strV == '' or strV == 'NULL':
		return None
	if not strV.isalpha():
		if strV.isdigit():
			return int(strV)
		elif '..' in strV:
			values = strV.split('..')
			result = []
			for value in values:
				result += [parseValue(value)]
			return result
		elif '.' in strV:
			return float(strV)	
	return strV




class CSVRow(object):
	def __init__(self):
		self._values = []
		self._names = []

	def __getattr__(self, name):
		#print 'geting',name
		if not name in self._names:
			msg = '{0}:{1}'
			raise AttributeError(msg.format(type(self).__name__, name))
		else:
			index = self._names.index(name)
			#print 'get',name,index
			return self._values[index]

	def __setattr__(self, name, value):
		#print 'setting', name, value
		if name in ['_names', '_values']:
			return  object.__setattr__(self, name, value)
		elif not name in self._names:
			self._values += [value]
			self._names += [name]
		else:
			i = self._names.index(name)
			#print 'set',name,i
			self._values[i] = value
		return value

	def __delattr__(self, name):
		if not name in self._names:
			msg = '{0}:{1}'
			raise AttributeError(msg.format(type(self).__name__, name))
		else:
			index = self._names.index(name)
			#print 'del',name,index
			del self._values[index]
			del self._names[index]

	def __str__(self):
		result = ''
		for value in self._values:
			if isinstance(value, list):
				result += '\t' + '..'.join([str(x) for x in value])
			else:
				result += '\t'+str(value)
		return result[1:]


class CSV:
	def __init__(self, sep='\t'):
		self._colList = []
		self.file = None
		self.separator = sep
		

	def __iter__(self):
		return self

	def open(self, file_name, has_header=True):
		self.file = open(file_name)
		self.current_line = 0
		self.has_header = has_header
		if has_header:
			self.current_line = 1
			self._parse_header()


	def _parse_header(self):
		header = self.file.readline().strip()
		self._colList = header.split(self.separator)


	def colName(self, i):
		if self.has_header:
			return self._colList[i]
		else: 
			return 'V'+str(i+1) # just as in R, Yoha! 



	def next(self):
		line = self.file.readline()
		#print 'reading ', line
		if not line:
			self.file.close() 
			raise StopIteration

		cols = line.strip().split(self.separator)
		result = CSVRow()
		
		if self.has_header and len(cols) != len(self._colList):
			raise Exception("Incompatible line (" + str(self.current_line) + "): " + line)
		
		for i in xrange(len(cols)):
			value = parseValue(cols[i])

			setattr(result, self.colName(i), value)

		self.current_line += 1

		return result



	def create(self, file_name, colList, has_header=True):
		self.close()
		self.file = open(file_name, 'w')
		self._colList = colList
		self.current_line = 0
		if has_header:
			self.current_line = 1
			for col in colList:
				self.file.write(col+self.separator)
			self.file.write('\n')		

	def writeln(self, row):
		if len(row._names) != len(self._colList):
			raise Exception("Incompatible row (" + str(self.current_line) + "): " + str(row))
		self.file.write(str(row)+'\n')
		self.current_line += 1

	def close(self):
		if self.file:
			self.file.close()

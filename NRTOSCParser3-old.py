import sc, time, struct, math, os, string
from scosc import OSC
from scosc import tools
from scosc import osctools

class NRTOSCParser3:
	def __init__(self, sr=44100, anchor='/Users/kfl/'):
		self.lib = { 'bndl': self.padBytes('#bundle'), 'zero': self.padBytes(0), 'anchor': anchor }
	
	def absToOSCTimestamp(self, abs):
		return struct.pack('!LL', math.floor(abs), long(float(abs - long(abs)) * 4294967296L))
	
	def padBytes(self, val):
		if type(val) == type('str'):
			pad = (math.ceil(((len(val) + 1) / 4.0)) * 4)
			ba = bytearray(val)
			while len(ba) < pad:
				ba.append(0)
			return ba
		elif type(val) == type(1):
			return struct.pack('!i', val)
		elif type(val) == type(0.2):
			return struct.pack('!f', val)

	def processAndWriteFile(self, score, output):
		header = bytearray("\x00\x00\x00$#bundle\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x10/g_new\x00\x00,i\x00\x00\x00\x00\x00\x01\x00\x00\x00$#bundle\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x10/g_new\x00\x00,i\x00\x00\x00\x00\x00\x01")
		total_length = 0
		with open(output, 'wb') as f:
			f.write(header)
			total_length += len(header)
			max_timestamp = 1.0
			for line in score:
				bundle = bytearray("\x00\x00\x00\x00#bundle\x00")
				timestamp = self.absToOSCTimestamp(line[0])
				max_timestamp = max(max_timestamp, (line[0] + 0.01))
				bundle.extend(timestamp)
# 				print line[1]
				bundle.extend("\x00\x00\x00\x00")
				commands = line[1]
				bundle.extend(self.padBytes(commands[0]))
				types = ","
				for item in commands[1:]:
					if type(item) == type('str'):
						types += 's'
					elif type(item) == type(1):
						types += 'i'
					elif type(item) == type(0.2):
						types += 'f'
				types = self.padBytes(types)
# 				print "BUNDLE: " + bundle
				bundle.extend(types)
				for item in commands[1:]:
						translated = self.padBytes(item)
						#print "        ", translated
						bundle.extend(translated)
# 				print "***LENGTH: ", len(bundle)
				total_length += len(bundle)
# 				for c in bundle: print("::: ", c, chr(c))
				bundle[3] = len(bundle) - 4
				bundle[23] = len(bundle) - 24
				f.write(bundle)
			footer = bytearray("\x00\x00\x00\x1c#bundle\x00")
			
			footer.extend(self.absToOSCTimestamp(max_timestamp))
			
			footer += "\x00\x00\x00\x00\x00\x00\x00\x08\x00\x00\x00\x00,\x00\x00\x00"
			total_length += len(footer)
# 			print "footer:"
# 			for c in footer: print("::: ", c, chr(c))
			
			f.write(footer)
#			print "total_length: " + `total_length`

	def createNRTScore(self, sfpath = '/Users/kfl/76444.wav',
						index = 0,
						pBuf = 10,
						aBuf = 11,
# 						currBus = 10, # internally?
						tratio = 1.0, # arg?
						srate = 44100,
						duration = 1.0,
						params = None):
		
		cwd = os.path.abspath(sfpath)
		fname = string.split(os.path.basename(sfpath), '.')
		cwd = os.path.dirname(cwd)
# 		print 'CWD: ',  cwd
		
		## create a subdirectory 'md' if it doesn't already exist
		self.make_a_dir(cwd, 'md')
		self.make_a_dir(os.path.join(cwd, 'md'), `index`)
		os.chdir(cwd)
		self.make_a_dir(cwd, 'osc')
		self.make_a_dir(os.path.join(cwd, 'osc'), `index`)
		os.chdir(cwd)
		
		mdpath = os.path.join(cwd , 'md', `index`, (fname[0] + '.md.' + fname[1]))
# 		print 'mdpath: ', mdpath
		oscpath = os.path.join(cwd , 'osc', `index`, (fname[0] + '_sigmaAnalyzer.osc'))
# 		print 'oscpath: ', oscpath
		
		## the two alloc calls
		oscList = [[0.0, ['/b_allocReadChannel', pBuf, sfpath, 0, -1, '[0]']]]
		oscList += [[0.01, ['/b_alloc', aBuf, int(math.ceil((duration/0.04) / tratio)), 25]]]
		
		sdef = 'sigmaSynthAnalyzer'
				
		row = ['srcbufNum', pBuf, 'savebufNum', aBuf, 'outbus', 0, 'dur', duration, 'transp', tratio] + params
				
			
		oscList += [[0.02, (['/s_new', sdef, -1, 0, 0] + row)]]
		
		oscList += [[((duration / tratio) + 0.03), ['/b_write', aBuf, mdpath, 'wav', 'float32']]]
		## don't free any buffers (yet)
		oscList += [[((duration / tratio) + 0.04), ['/c_set', 0, 0]]]
	
#  		print '\nTHE LIST:: ' + `oscList` + '\n'
	
		self.processAndWriteFile(oscList, oscpath)
		
# 		print self.lib['anchor']
		os.chdir(self.lib['anchor']) # always return pwd to the anchor dir
		
		return oscpath, mdpath

	def make_a_dir(self, parent, dirname):
		try:
			os.chdir(parent)
		except OSError:
# 			print 'OSERROR: ', parent, 'does not exist!'
			return None
		try:
			os.mkdir(dirname)
		except OSError:
# 			print 'OSERROR: ', dirname, 'already exists!'
			pass
		
import sc, random, contextlib, wave, os, math
import shlex, subprocess, signal
import NRTOSCParser3

import numpy as np
import scipy.signal


# generator class for weighted random numbers
#
# Pass in one or the other:
# - weights: custom weights array
# - size: size of "standard" weights array that algo should make on its own
#
# call next to actually make the random selection
#
class RandomGenerator_8Bit(object):

	def __init__(self, initval=-1):
		if initval >= 0:
			self.val = initval
		else:
			self.val = random.randint(0,128)
		
	def next(self, scale=1.0):
		self.val = random.randint(0,128)

	def __call__(self): return self.next()


# helper function
def midi2hz(m): return pow(2.0, (m/12.0))

# slot assignments for sigmaSynth
ALPHA = 0
C_DELAY = 1
BETA = 2
D_MULT = 3
GAMMA = 4
MS_BINS = 5

class GenomicExplorer:

	def __init__(self, anchor, sfilenames, size=20, kdepth=10): #, start_state=[1.0, 0.0, 1.0, 1.0, 1.0, 0.0]
				
		self.anchor = anchor
		self.sfpaths = [(anchor + '/snd/' + sfilename) for sfilename in sfilenames]
		self.filenames = sfilenames
		self.sfinfos = []
		
		for path in self.sfpaths:
			with contextlib.closing(wave.open(path,'r')) as f: 
				frames = f.getnframes()
				rate = f.getframerate()
				dur = frames/float(rate)
				self.sfinfos += [{'rate':rate,'dur':dur}]
		
		self.mutation_prob = 0.05
		#self.xover_prob = 0.05
		self.depth = kdepth
 		# 'alpha', 'c_delay', 'beta', 'd_mult', 'gamma', 'ms_bins'
		self.parser = NRTOSCParser3.NRTOSCParser3(anchor=self.anchor)
		self.rawtable, self.rawmaps, self.dists = dict(), dict(), dict()

		print self.sfpaths
		print self.sfinfos

		self.init_population(size=size)

	
	def init_population(self, size):
		self.population = []
		for n in range(size):
  			#start = [random.randrange(500, 1000)*0.001, random.randrange(0,50)*0.001, random.randrange(500, 1000)*0.001, random.randrange(100,1000)*0.01, random.randrange(500, 1000)*0.001, random.randrange(0,5000)*0.01]
 			self.population += [Genome()] #random seed
# 			self.population += [Genome(starter)]
		self.population[0] = Genome(values=[0,0,0,0,0,0])
		self.analyze_individual(0)
		self.activate_raw_data(0)
		self.compare_all_individuals(aflag=True)
			
	def mutate_pop(self):
		for indiv in range(1, len(self.population)):
			if random.random() < self.mutation_prob:
				print "indiv: ", indiv
				self.population[ indiv ].mutate()
				self.do_update_cascade(indiv)
				
	def do_update_cascade(self, index, clearedits=False):
		if clearedits is True:
			self.population[ index ].edits = 0
		else:
			self.population[ index ].edits += 1
		self.analyze_individual( index )
		self.activate_raw_data( index )
		self.compare_individual_chi_squared( index )
# 		self.compare_individual( index )
	
	def mate(self, a, b, kill_index):
		
# 		cut = random.randint(0,5)
		offspring = None		
		
		if random.random() < 0.5:
			offspring = self.population[a].values[:]
			
		else:
			offspring = self.population[b].values[:]
		
		# basic gene selection from 2 parents
		for i in range(6):
			if random.random() < 0.5:
				offspring[i] = self.population[a].values[i]
			else:
				offspring[i] = self.population[b].values[i]
		
		self.population[kill_index] = Genome(offspring)
		
		self.do_update_cascade(kill_index, True)
	
	def sort_by_distances(self, depth):
		sorted_dists = [[k, self.dists[k], self.population[k].age, self.population[k].edits] for k in sorted(self.dists.keys())]
		sorted_dists = sorted(sorted_dists[1:], key = lambda row: row[1]) # + (maxedits - row[3])))
		return sorted_dists[:depth], sorted_dists[(-1*depth):]
		
	def reproduce(self, depth=25):
		
		kills, duplicates = self.sort_by_distances(depth)
		print 'depth: ', depth
		# depth # of times: choose 2 random parents to mate and overwrite replacement in unfit individual's slot
		for n in range(depth):
			print 'num. duplicates: ', len(duplicates)
			aidx = duplicates[ random.randint(0, depth-1) ][0]
			bidx = duplicates[ random.randint(0, depth-1) ][0]
			
			kidx = kills[ random.randint(0, depth-1) ][0]
			
			self.mate(aidx, bidx, kidx)
	
	def age_pop(self):
		for i in range(len(self.population)): self.population[i].age += 1
	
	def iterate(self, iters=1):
		sc.quit()
		for iter in range(iters):
			self.age_pop()
			self.mutate_pop()
# 			self.crossover()
			if (iter%20)==0:
				print self.population[0].age
				self.reproduce(self.depth)
	
	def print_all_individuals(self):
		print '== pop ==========================='
		for g in self.population: print g
	
	def start_sc(self):
		try:
			sc.start(verbose=1, spew=1, startscsynth=1)
		except OSError: # in case we've already started the synth
			print 'QUIT!'
			sc.quit()
		print 'sfpath: ', self.sfpath
		for i, sfpath in enumerate(self.sfpaths):
			bnum = sc.loadSnd(sfpath, wait=False)
			print 'bnum: ', bnum
			self.infos[i]['bnum'] = bnum
		return 1
	
	# |outbus=20, srcbufNum, start=0.0, dur=1.0, transp=1.0, c_delay=0.0, c_decay=0.0, d_mult=1.0, d_amp=0.7, ms_bins=0, alpha=1, beta=1, gamma=1|
	def play_genome(self, index):
		
		vals = self.population[index].realvalues
		if vals[C_DELAY] < 1.0:
			cdelay = 0.0
		else:
			cdelay = vals[C_DELAY]
		decay = 0.9

		tr = self.population[index].tratio
		
		if index == 0:
			slot = 0
		else:
			slot = 1

		print '===================\n', self.infos[slot]['dur']
		
		sc.Synth('sigmaSynth', 
			args=[
				'srcbufNum', self.infos[slot]['bnum'], 
				'start', 0,
				'dur', self.infos[slot]['dur']*1000,
				'transp', tr,
				'c_delay', cdelay,
				'c_decay', decay,
				'd_mult', vals[D_MULT],
				'ms_bins', vals[MS_BINS],
				'alpha', vals[ALPHA],
				'beta', vals[BETA],
				'gamma', vals[GAMMA]])
	
	def analyze_individual(self, index):

#		oscpath = os.path.join(self.anchor, 'snd', 'osc', `index`, (os.path.splitext(self.filename)[0] + '_sigmaAnalyzer.osc'))
# 		mdpath = os.path.join(self.anchor, 'snd', 'md', `index`, self.filename)
		
		vals = self.population[index].realvalues
		if vals[C_DELAY] < 1.0:
			cdelay = 0.0
		else:
			cdelay = vals[C_DELAY]
		decay = 0.9
		
		tr = self.population[index].tratio

		if index == 0:
			slot = 0
		else:
			slot = 1

		oscpath, mdpath = self.parser.createNRTScore(self.sfpaths[slot],
							index=index, 
							tratio=tr,
							srate=self.sfinfos[slot]['rate'],
							duration=self.sfinfos[slot]['dur'],
							params=[
								'c_delay', cdelay,
								'c_decay', decay,
								'd_mult', vals[D_MULT],
								'ms_bins', vals[MS_BINS],
								'alpha', vals[ALPHA],
								'beta', vals[BETA],
								'gamma', vals[GAMMA]])

		cmd = 'scsynth -N ' + oscpath + ' _ _ 44100 WAVE float32 -o 1'
# 		print cmd
		args = shlex.split(cmd)
		p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE) #, shell=True, close_fds=True)
		
# 		print 'PID: ', p.pid
		rc = p.wait()
				
# 		print 'RC: ', rc
		if rc == 1:
			num_frames = int(math.ceil(self.sfinfos[slot]['dur'] / 0.04 / tr))
# 			print 'num frames: ', num_frames
			self.rawtable[index] = (mdpath, num_frames)
		
# 		print self.rawtable
	
	def render_individual(self, index):
	
		vals = self.population[index].realvalues
		if vals[C_DELAY] < 1.0:
			cdelay = 0.0
		else:
			cdelay = vals[C_DELAY]
		decay = 0.9
		
		tr = self.population[index].tratio

		if index == 0:
			slot = 0
		else:
			slot = 1

		oscpath, mdpath = self.parser.createNRTScore(self.sfpaths[slot], 
							index=index, 
							tratio=tr,
							srate=self.sfinfos[slot]['rate'],
							duration=self.sfinfos[slot]['dur'],
							params=[
								'c_delay', cdelay,
								'c_decay', decay,
								'd_mult', vals[D_MULT],
								'ms_bins', vals[MS_BINS],
								'alpha', vals[ALPHA],
								'beta', vals[BETA],
								'gamma', vals[GAMMA]])

		cmd = 'scsynth -N ' + oscpath + ' _ ' + os.path.join(self.anchor, 'snd', 'out', (str(index) + '.aiff')) + ' 44100 AIFF int16 -o 1'
		args = shlex.split(cmd)
		p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE) #, shell=True, close_fds=True)
		rc = p.wait()
		if rc == 1:
		 	print 'SUCCESS: ', os.path.join(self.anchor, 'snd', 'out', (str(index) + '.aiff'))
		 	rc = 0
		else:
			return None
		cmd = 'sox -b 16 ' + os.path.join(self.anchor, 'snd', 'out', (str(index) + '.aiff')) + ' ' + os.path.join(self.anchor, 'snd', 'out', (str(index) + '.wav')) # + '; rm ' + os.path.join(self.anchor, 'snd', 'out', (str(index) + '.aiff'))
		print cmd
		args = shlex.split(cmd)
		p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE) #, shell=True, close_fds=True)
		rc = p.wait()
		print rc
		if rc == 1: ' DOUBLE SUCCESS!!'
	
	def activate_raw_data(self, index):
		
		mdpath = self.rawtable[index][0]
		num_frames = self.rawtable[index][1]
		self.rawmaps[index] = np.memmap(mdpath, dtype=np.float32, mode='r', offset=272, shape=(num_frames, 25))
	
	"""
	COMPARE_ALL_INDIVIDUALS:
		... to individual in slot 0!
	"""
	def compare_all_individuals(self, aflag=False):
		for i in range(1, len(self.population)):
 			if aflag:
 				self.analyze_individual(i)
				self.activate_raw_data(i)
 			self.compare_individual_chi_squared(i)
#			self.compare_individual(i)
		print self.dists
		return self.dists
	"""
	COMPARE_INDIVIDUAL:
		... to individual in the slot that is stipulated by the arg zeroindex!
		-	by convention, we should usually put what we are comparing to in slot 0
	"""
	def compare_individual(self, index, zeroindex=0):
		i_length = self.rawmaps[index].shape[0]
		zr0_length = self.rawmaps[zeroindex].shape[0]
		print i_length, ' | ', zr0_length

		# i1_length = self.rawmaps[index-1].shape[0] ## <--- NEIGHBOR comparison
		# print i_length, ' | ', i1_length, ' | ', zr0_length

		# based on length comparison, resample the mutated individuals so that they are same length as the zeroth individual (that does not mutate)
		# if indiv. is longer, resample indiv., take abs. diff., sum, div. by length
		if zr0_length < i_length:
			zero_dist = float(np.sum(np.abs(scipy.signal.signaltools.resample(self.rawmaps[index], zr0_length, window='hanning') - self.rawmaps[0]))) / float(zr0_length)
# 			print self.dists[index]
		# if zeroth indiv. is longer, resample zeroth indiv., take abs. diff., sum, div. by length, then do same comparison with "neighbor"
		elif i_length < zr0_length:
			zero_dist = float(np.sum(np.abs(self.rawmaps[index] - scipy.signal.signaltools.resample(self.rawmaps[0], float(i_length), window='hanning')))) / float(i_length)
		else:
		# otherwise, take abs. diff., sum, div. by length, then do same comparison with "neighbor"
# 			print 'ZERO'
			zero_dist = float(np.sum(np.abs(self.rawmaps[index][:,1:] - self.rawmaps[0][:,1:]))) / float(zr0_length)
			
			### CHECK THIS DISTANCE CALCULATION!!!!!!
			power_dist = float(np.sqrt(np.sum(np.abs(self.rawmaps[index][:,0] - self.rawmaps[0][:,0])))) / float(zr0_length)
			print (zero_dist, (power_dist * 10.0))
			zero_dist += (power_dist * 10.0)
		
# 		if i1_length < i_length:
# 			neighbor_dist = float(np.sum(np.abs(scipy.signal.signaltools.resample(self.rawmaps[index-1], i_length, window='hanning') - self.rawmaps[index]))) / float(i_length)
# 		elif i_length < i1_length:
# 			neighbor_dist = float(np.sum(np.abs(self.rawmaps[index-1] - scipy.signal.signaltools.resample(self.rawmaps[index], i1_length, window='hanning')))) / float(i1_length)
# 		else:
# 			print 'ZERO-NEIGHBOR'
# 			neighbor_dist = float(np.sum(np.abs(self.rawmaps[index-1] - scipy.signal.signaltools.resample(self.rawmaps[index], i1_length, window='hanning')))) / float(i1_length)
		
# 		self.dists[index] = zero_dist + neighbor_dist
		self.dists[index] = zero_dist
	
	
	def compare_individual_chi_squared(self, index):
		i_length = self.rawmaps[index].shape[0]
		i1_length = self.rawmaps[index-1].shape[0]
		zr0_length = self.rawmaps[0].shape[0]
# 		print i_length, '|', zr0_length
		# based on length comparison, resample the mutated individuals so that they are same length as the zeroth individual (that does not mutate)
		# if indiv. is longer, resample indiv., take abs. diff., sum, div. by length
		if zr0_length < i_length:
			zero_dist = scipy.stats.mstats.chisquare(scipy.signal.signaltools.resample(self.rawmaps[index], zr0_length, window='hanning'), self.rawmaps[0])
# 			print self.dists[index]
		# if zeroth indiv. is longer, resample zeroth indiv., take abs. diff., sum, div. by length, then do same comparison with "neighbor"
		elif i_length < zr0_length:
			zero_dist = scipy.stats.mstats.chisquare(self.rawmaps[index], scipy.signal.signaltools.resample(self.rawmaps[0], i_length, window='hanning'))
		else:
		# otherwise, take abs. diff., sum, div. by length, then do same comparison with "neighbor"
			 print 'CHI-ZERO'
			 zero_dist = scipy.stats.mstats.chisquare(self.rawmaps[index], self.rawmaps[0])
		
		if i1_length < i_length:
			neighbor_dist = scipy.stats.mstats.chisquare(scipy.signal.signaltools.resample(self.rawmaps[index-1], i_length, window='hanning') - self.rawmaps[index])
		elif i_length < i1_length:
			neighbor_dist = scipy.stats.mstats.chisquare(self.rawmaps[index-1], scipy.signal.signaltools.resample(self.rawmaps[index], i1_length, window='hanning'))
		else:
			print 'CHI-NEIGHBOR'
			neighbor_dist = scipy.stats.mstats.chisquare(self.rawmaps[index-1], scipy.signal.signaltools.resample(self.rawmaps[index], i1_length, window='hanning'))
		
		nsum = np.sum(np.abs(neighbor_dist[0].data[:24]))
		zsum = np.sum(np.abs(zero_dist[0].data[:24]))
		nasum = neighbor_dist[0].data[24]
		zasum = zero_dist[0].data[24]
		
		self.dists[index] = nsum + zsum - (24.0 * nasum) - (24.0 * zasum)


class Genome:

	def __init__(self, values=None, slotranges=[[1.0,0.5],[0.0,0.05],[1.0, 0.5],[1.0,10.],[1.0,0.5],[0.0,50.]]):
		
		# """
		# 'alpha',   'c_delay', 'beta',   'd_mult', 'gamma',  'ms_bins'
		# [[1.0,0.5],[0.0,0.05],[1.0,0.5],[1.0,10.],[1.0,0.5],[0.0,50.]]
		# """
		
		self.tratio	= 1.0		# CHECK THIS... WHY IS IT HERE/in Hertz!!! ???
		
		self.boundaries = slotranges
		self.generators = [RandomGenerator_8Bit(-1) for n in range(6)] ### CONSTANT WARNING
		#StaticGenerator_8Bit(VAL) ???
		
		if values is None:
			print 'values is None, generators are seeded randomly!'
			self.values = [gen.val for gen in self.generators]
		else:
			self.values = values
		self.bitlength = len(self.values) * 8
		self.binarystring = vals_to_binarystring(self.values)
# 		print self.values
# 		print type(self.values[0])
		self.realvalues = [lininterp(val,self.boundaries[i]) for i,val in enumerate(self.values)]
		
		self.age = 0
		self.edits = 0
	
	def __repr__(self):
		print tuple(self.values)
		print ((self.age, self.edits) + tuple(self.values) + tuple(self.binarystring))
		return "%9i/%9i || %.6f|%.6f|%.6f|%.6f|%.6f|%.6f" % ((self.age, self.edits) + tuple(self.realvalues)) # + tuple(self.binarystring)
	
	def mutate(self):
		pos = random.randint(0,(self.bitlength-1))
		# flip bit
		print 'bit flipped to: ', abs(1 - int(self.binarystring[pos],2))
		self.binarystring = substitute_char_in_string(self.binarystring, pos, abs(1 - int(self.binarystring[pos],2)))
		# recalc binary string
		self.values = binarystring_to_vals(self.binarystring)
		print "values: ", self.values
		self.realvalues = [lininterp(val,self.boundaries[i]) for i,val in enumerate(self.values)]
	
# 	def xover_sub(self, pos, incomingSeq, headortail=0):
# 		if headortail == 0:
# 			print '<<>> ', self.binarystring
# 			print '<<>> ', pos
# 			print '<<>> ', incomingSeq
# 			self.binarystring = incomingSeq[:pos] + self.binarystring[pos:]
# 		else:
# 			print '<<>> ', self.binarystring
# 			print '<<>> ', pos
# 			print '<<>> ', incomingSeq
# 			self.binarystring = self.binarystring[:pos] + incomingSeq[:(len(self.binarystring)-pos)]
# 		# recalc binary string
# 		print '==== ', self.binarystring
# 		self.values = binarystring_to_vals(self.binarystring)
# 		print "values: ", self.values
# 		self.realvalues = [lininterp(val,self.boundaries[i]) for i,val in enumerate(self.values)]


def lininterp(val,bounds=[0.,1.]):
	return (((val/128.0)*(bounds[1]-bounds[0]))+bounds[0])

def substitute_char_in_string(s, p, c):
	l = list(s)
	l[p] = str(c)
	return "".join(l)

# def substitute_string_head(s, p, snew):
# 	s1 = snew[:]
# 	print '++++ ', s1
# 	s2 = s[p:]
# 	print '++++ ', s2
# 	return (s1 + s2)[:len(s)]
# 
# def substitute_string_tail(s, p, snew):
# 	s1 = s[:p]
# 	print '==== ', s1
# 	print len(s)
# 	print p
# 	s2 = snew[:(len(s)-p)]
# 	print '==== ', s2
# 	return (s1 + s2)[:len(s)]

def vals_to_binarystring(vals = [0, 0, 0, 0, 0]):
	return ''.join((("{0:08b}".format(val)) for val in vals))

# never a '0bXXX' string!
def binarystring_to_vals(binstring):
	mystring = binstring[:]
	length = len(mystring) / 8 # ignore the last digits if it doesn't chunk into 8-item substrings
	res = []
# 	print mystring[(n*8):((n+1)*8)]
	return [int(mystring[(n*8):((n+1)*8)], 2) for n in range(length)]
	

# if __name__=='__main__':
# 	genex = GenomicExplorer('/Users/kfl/dev/python/sc-0.3.1/genomic', 'test.wav')
# 	genex.analyze_genome(1)
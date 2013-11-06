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
class WeightedRandomGenerator(object):
	"""
	mu - should be 0!
	stdev - set very lo (at most 1/4 of the desired range)

	"""
	def __init__(self, mn=0.0, stdev=1.0, lo=0.0, hi=1.0, initval=0):
		self.mean = mn
		self.stdev = stdev
		self.lo = lo
		self.hi = hi
		self.val = initval
		
	def next(self, scale=1.0):
		rnd = random.gauss(self.mean, self.stdev)
		self.val = min(max((self.val+rnd), self.lo), self.hi)
		return self.val

	def __call__(self): return self.next()


# helper function
def midi2hz(m): return pow(2.0, (m/12.0))

ALPHA = 0
C_DELAY = 1
BETA = 2
D_MULT = 3
GAMMA = 4
MS_BINS = 5

class GenomicExplorer:

	def __init__(self, anchor, sfilename, size=10, start_state=[1.0, 0.0, 1.0, 1.0, 1.0, 0.0]):
				
		self.anchor = anchor
		self.sfpath = anchor + '/snd/' + sfilename
		self.filename = sfilename
		
		with contextlib.closing(wave.open(self.sfpath,'r')) as f: 
			self.frames = f.getnframes()
			self.rate = f.getframerate()
			self.dur = self.frames/float(self.rate)
		
		self.mutation_prob = 0.01
		self.depth = 25
 		# 'alpha', 'c_delay', 'beta', 'd_mult', 'gamma', 'ms_bins'
		self.parser = NRTOSCParser3.NRTOSCParser3(anchor=self.anchor)
		self.rawtable, self.rawmaps, self.dists = dict(), dict(), dict()

		self.init_population(size=size, starter=start_state)

	
	def init_population(self, size, starter=None):
		self.population = []
		for n in range(size):
  			start = [random.randrange(500, 1000)*0.001, random.randrange(0,50)*0.001, random.randrange(500, 1000)*0.001, random.randrange(100,1000)*0.01, random.randrange(500, 1000)*0.001, random.randrange(0,5000)*0.01]
 			self.population += [Genome(start)]
# 			self.population += [Genome(starter)]
		self.population[0] = Genome(starter)
		self.analyze_individual(0)
		self.activate_raw_data(0)
		self.compare_all_individuals(aflag=True)
			
	def mutate(self):
		
		if random.random() < self.mutation_prob:
			indiv = random.randint(1, len(self.population)-1)
			self.population[ indiv ].mutate()
			self.population[ indiv ].edits += 1
			self.analyze_individual( indiv )
			self.activate_raw_data( indiv )
# 			self.compare_individual_chi_squared( indiv )
			self.compare_individual( indiv )
	
	def mate(self, a, b, kill_index):
		
		cut = random.randint(0,5)
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
		
		self.analyze_individual(kill_index)
		self.activate_raw_data(kill_index)
# 		self.compare_individual_chi_squared( indiv )
		self.compare_individual(kill_index)
	
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
			self.mutate()
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
		# in case we've already started the synth
		except OSError:
			print 'QUIT!'
			sc.quit()
		print 'sfpath: ', self.sfpath
 		self.bnum = sc.loadSnd(self.sfpath, wait=False)
 		print 'bnum: ', self.bnum

# |outbus=20, srcbufNum, start=0.0, dur=1.0, transp=1.0, c_delay=0.0, c_decay=0.0, d_mult=1.0, d_amp=0.7, ms_bins=0, alpha=1, beta=1, gamma=1|

	def play_genome(self, index):
		
		vals = self.population[index].values
		if vals[C_DELAY] < 1.0:
			cdelay = 0.0
		else:
			cdelay = vals[C_DELAY]
		decay = 0.9

		tr = self.population[index].tratio
		
		print '===================\n', self.dur
		
		sc.Synth('sigmaSynth', 
			args=[
				'srcbufNum', self.bnum, 
				'start', 0,
				'dur', self.dur*1000,
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
		
		vals = self.population[index].values
		if vals[C_DELAY] < 1.0:
			cdelay = 0.0
		else:
			cdelay = vals[C_DELAY]
		decay = 0.9
		
		tr = self.population[index].tratio

		oscpath, mdpath = self.parser.createNRTScore(self.sfpath, 
							index=index, 
							tratio=tr,
							srate=self.rate,
							duration=self.dur,
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
			num_frames = int(math.ceil(self.dur / 0.04 / tr))
# 			print 'num frames: ', num_frames
			self.rawtable[index] = (mdpath, num_frames)
		
# 		print self.rawtable
	
	def render_individual(self, index):
	
		vals = self.population[index].values
		if vals[C_DELAY] < 1.0:
			cdelay = 0.0
		else:
			cdelay = vals[C_DELAY]
		decay = 0.9
		
		tr = self.population[index].tratio

		oscpath, mdpath = self.parser.createNRTScore(self.sfpath, 
							index=index, 
							tratio=tr,
							srate=self.rate,
							duration=self.dur,
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
			cmd = 'sox -b 16 ' + os.path.join(self.anchor, 'snd', 'out', (str(index) + '.aiff')) + ' ' + os.path.join(self.anchor, 'snd', 'out', (str(index) + '.wav')) + '; rm ' + os.path.join(self.anchor, 'snd', 'out', (str(index) + '.aiff'))
			args = shlex.split(cmd)
			p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE) #, shell=True, close_fds=True)
			rc = p.wait()
			if rc == 1: ' DOUBLE SUCCESS!!'
	
	def activate_raw_data(self, index):
		
		mdpath = self.rawtable[index][0]
		num_frames = self.rawtable[index][1]
		self.rawmaps[index] = np.memmap(mdpath, dtype=np.float32, mode='r', offset=272, shape=(num_frames, 25))
	
	def compare_all_individuals(self, aflag=False):
		for i in range(1, len(self.population)):
 			if aflag:
 				self.analyze_individual(i)
				self.activate_raw_data(i)
# 			self.compare_individual_chi_squared(i)
			self.compare_individual(i)
		return self.dists
	
	def compare_individual(self, index):
		i_length = self.rawmaps[index].shape[0]
		i1_length = self.rawmaps[index-1].shape[0]
		zr0_length = self.rawmaps[0].shape[0]
		print i_length, ' | ', i1_length, ' | ', zr0_length
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
			power_dist = float(np.sum(self.rawmaps[index][:,0] - self.rawmaps[0][:,0])) / float(zr0_length)
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

	def __init__(self, starter):

		print 'starter: ', starter
		
		self.values = [starter[ALPHA], starter[C_DELAY], starter[BETA], starter[D_MULT], starter[GAMMA], starter[MS_BINS]]
		self.tratio	= 1.0		# in Hertz!!!
		
		self.generators = [
			WeightedRandomGenerator(0.0, 	0.1, 	0.5,	1.0,	self.values[ALPHA]),
			WeightedRandomGenerator(0.0, 	0.0025,	0.0, 	0.05,	self.values[C_DELAY]),
			WeightedRandomGenerator(0.0, 	0.05, 	0.5,	1.0,	self.values[BETA]),
			WeightedRandomGenerator(0.0, 	0.5, 	1.0,	10.,	self.values[D_MULT]),
			WeightedRandomGenerator(0.0, 	0.05, 	0.5,	1.0,	self.values[GAMMA]),
			WeightedRandomGenerator(0.0, 	2.5, 	0.0,	50.,	self.values[MS_BINS])
		]
		self.age = 0
		self.edits = 0
	
	def __repr__(self):
		return "%9i/%9i || %.6f|%.6f|%.6f|%.6f|%.6f|%.6f" % (self.age, self.edits, self.values[ALPHA], self.values[C_DELAY], self.values[BETA], self.values[D_MULT], self.values[GAMMA], self.values[MS_BINS])
	
	def mutate(self):
		choice = random.randint(0,5)

		self.values[choice] = self.generators[choice].next()
	

# if __name__=='__main__':
# 	genex = GenomicExplorer('/Users/kfl/dev/python/sc-0.3.1/genomic', 'test.wav')
# 	genex.analyze_genome(1)
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

T12 = 0
T3 = 1
T1 = 2
ALPHA = 3
C_DELAY = 4
BETA = 5
D_MULT = 6
GAMMA = 7
MS_BINS = 8

class GenomicExplorer:

	def __init__(self, anchor, sfilename, size=10, start_state=[0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0]):
				
		self.anchor = anchor
		self.sfpath = anchor + '/snd/' + sfilename
		self.filename = sfilename
		
		with contextlib.closing(wave.open(self.sfpath,'r')) as f: 
			self.frames = f.getnframes()
			self.rate = f.getframerate()
			self.dur = self.frames/float(self.rate)
		
		self.mutation_prob = 0.05
		self.xover_prob = 0.05
		self.replication_prob = 0.05
# 		self.map=['transp', 'alpha', 'c_delay', 'beta', 'd_mult', 'gamma', 'ms_bins']
# 		self.limits = {
# 			'transp'	:	(-12.,	12., 	12.,	0.),
# 			'alpha'		:	(0.,	1.,		0.5,	0.5),
# 			'c_delay'	:	(0.0, 	0.05,	0.025,	0.025),
# 			'beta'		:	(0.,	1.,		0.5,	0.5),
# 			'd_mult'	:	(1.,	10.,	4.5,	5.5),
# 			'gamma'		:	(0.,	1.,		0.5,	0.5),
# 			'ms_bins'	:	(0.,	50.,	25.,	25.)
# 		}
		self.parser = NRTOSCParser3.NRTOSCParser3(anchor=self.anchor)
		self.rawtable, self.rawmaps, self.dists = dict(), dict(), dict()

		self.init_population(size=size, starter=start_state)

	
	def init_population(self, size, starter=None):
		self.population = []
		for n in range(size):
#  			start = [random.randrange(-1,1), random.randrange(-3,3), random.randrange(-7500,7500)*0.0001, random.random(), random.randrange(0,50)*0.001, random.random(), random.randrange(100,1000)*0.01, random.random(), random.randrange(0,5000)*0.01]
#  			self.population += [Genome(start)]
			self.population += [Genome(starter)]
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
			self.compare_individual( indiv )
	
	def crossover(self, a, b):
		
		if random.random() < self.xover_prob:
# 			a = random.randint(1, len(self.population)-1)
# 			b = random.randint(1, len(self.population)-1)
			cut = random.randint(0,8)
			
	# 			print a, '|', b, '|', cut
			self.population[a].edits += 1
			self.population[b].edits += 1
			
			temp = self.population[a].values
			self.population[a].values[:cut] = self.population[b].values[:cut]
			self.population[b].values[cut:] = temp[cut:]
	
			self.analyze_individual(a)
			self.analyze_individual(b)
			self.activate_raw_data(a)
			self.activate_raw_data(b)
			self.compare_individual(a)
			self.compare_individual(b)

	
	def replicate(self, depth=25):
		sorted_dists = [[k, self.dists[k], self.population[k].age, self.population[k].edits] for k in sorted(self.dists.keys())]
		#maxedits = int(np.max(np.array(sorted_dists[1:])[:,3]))
		#print 'max edits: ', maxedits
		sorted_dists = sorted(sorted_dists[1:], key = lambda row: row[1]) # + (maxedits - row[3])))
# 		print ''
# 		print '==================='
# 		print sorted_dists[1:]
# 		print ''
		kills = sorted_dists[:depth]
		duplicates = sorted_dists[(-1*depth):]
		
# 		print kills
# 		print duplicates
		
		# either replicate or xover!!!
		for n in range(depth):
			if random.random() < self.replication_prob:
# 				print 'REP!'
				seed = self.population[duplicates[n][0]]
				
				self.population[kills[n][0]] = Genome(seed.values)
				self.dists[ kills[n][0] ] =  self.dists[ duplicates[n][0] ]
				# perform cross-over with newly born indiv.
				
				bnum = random.randint(0, depth-1)
				if random.random < 0.5:
					self.crossover(kills[n][0], duplicates[bnum][0])
				else:
					self.crossover(duplicates[bnum][0], kills[n][0])
			else:
				bnum = random.randint(0, depth-1)
				self.crossover(duplicates[n][0], duplicates[bnum][0])
	
	def age_pop(self):
		for i in range(len(self.population)): self.population[i].age += 1
	
	def iterate(self, iters=1):
		sc.quit()
		for iter in range(iters):
			self.age_pop()
			self.mutate()
# 			self.crossover()
			if (iter%20)==0:
				self.replicate()
	
	def print_all_individuals(self):
		print '== pop ==========================='
		for g in self.population: print g
	
	def start_sc(self):
		sc.start()
 		self.bnum = sc.loadSnd(self.sfpath, wait=False)

# |outbus=20, srcbufNum, start=0.0, dur=1.0, transp=1.0, c_delay=0.0, c_decay=0.0, d_mult=1.0, d_amp=0.7, ms_bins=0, alpha=1, beta=1, gamma=1|

	def play_genome(self, index):
		
		vals = self.population[index].values
		if vals[C_DELAY] < 1.0:
			cdelay = 0.0
		else:
			cdelay = vals[C_DELAY]
		decay = 0.9

		tr = self.population[index].tratio
		
		sc.Synth('sigmaSynth', 
			arglist=[
				'srcbufNum', self.bnum, 
				'start', 0,
				'dur', self.dur,
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
	
	def activate_raw_data(self, index):
		
		mdpath = self.rawtable[index][0]
		num_frames = self.rawtable[index][1]
		self.rawmaps[index] = np.memmap(mdpath, dtype=np.float32, mode='r', offset=272, shape=(num_frames, 25))
	
	def compare_all_individuals(self, aflag=False):
		for i in range(1, len(self.population)):
 			if aflag:
 				self.analyze_individual(i)
				self.activate_raw_data(i)
			self.compare_individual(i)
		return self.dists
	
	def compare_individual(self, index):
		i_length = self.rawmaps[index].shape[0]
		i1_length = self.rawmaps[index-1].shape[0]
		zr0_length = self.rawmaps[0].shape[0]
# 		print i_length, '|', zr0_length
		# based on length comparison, resample the mutated individuals so that they are same length as the zeroth individual (that does not mutate)
		# if indiv. is longer, resample indiv., take abs. diff., sum, div. by length
		if zr0_length < i_length:
			zero_dist = float(np.sum(np.abs(scipy.signal.signaltools.resample(self.rawmaps[index], zr0_length, window='hanning') - self.rawmaps[0]))) / float(zr0_length)
# 			print self.dists[index]
		# if zeroth indiv. is longer, resample zeroth indiv., take abs. diff., sum, div. by length, then do same comparison with "neighbor"
		elif i_length < zr0_length:
			zero_dist = float(np.sum(np.abs(self.rawmaps[index] - scipy.signal.signaltools.resample(self.rawmaps[0], i_length, window='hanning')))) / float(i_length)
		else:
		# otherwise, take abs. diff., sum, div. by length, then do same comparison with "neighbor"
			 zero_dist = float(np.sum(np.abs(self.rawmaps[index] - self.rawmaps[0]))) / float(zr0_length)
		
		if i1_length < i_length:
			neighbor_dist = float(np.sum(np.abs(scipy.signal.signaltools.resample(self.rawmaps[index-1], i_length, window='hanning') - self.rawmaps[index]))) / float(i_length)
		elif i_length < i1_length:
			neighbor_dist = float(np.sum(np.abs(self.rawmaps[index-1] - scipy.signal.signaltools.resample(self.rawmaps[index], i1_length, window='hanning')))) / float(i1_length)
		else:
			neighbor_dist = float(np.sum(np.abs(self.rawmaps[index-1] - scipy.signal.signaltools.resample(self.rawmaps[index], i1_length, window='hanning')))) / float(i1_length)
		
		self.dists[index] = zero_dist + neighbor_dist


# 	def compare_individual_to_mean(self, index):
# 		i_length = self.rawmaps[index].shape[0]
# 		i1_length = self.rawmaps[index-1].shape[0]
# 		zr0_length = self.rawmaps[0].shape[0]
# 
# 		if zr0_length < i_length:
# 			zero_dist = float(np.sum(np.abs(scipy.signal.signaltools.resample(self.rawmaps[index], zr0_length, window='hanning') - self.rawmaps[0]))) / float(zr0_length)
# # 			print self.dists[index]
# 		# if zeroth indiv. is longer, resample zeroth indiv., take abs. diff., sum, div. by length, then do same comparison with "neighbor"
# 		elif i_length < zr0_length:
# 			zero_dist = float(np.sum(np.abs(self.rawmaps[index] - scipy.signal.signaltools.resample(self.rawmaps[0], i_length)))) / float(i_length)
# 		else:
# 		# otherwise, take abs. diff., sum, div. by length, then do same comparison with "neighbor"
# 			 zero_dist = float(np.sum(np.abs(self.rawmaps[index] - self.rawmaps[0]))) / float(zr0_length)
# 
# 		# compute population mean
		

class Genome:

	def __init__(self, starter):

		self.values = [(round(starter[0]) * 12.0), round(starter[T3]), starter[T1], starter[ALPHA], starter[C_DELAY], starter[BETA], starter[D_MULT], starter[GAMMA], starter[MS_BINS]]
		self.transp = self.values[0] + self.values[1] + self.values[2]
		self.tratio	= midi2hz(self.transp)		# in Hertz!!!
		
		self.generators = [
			WeightedRandomGenerator(0.0, 	0.25, 	-1.0,	1.0,	self.values[T12]),
			WeightedRandomGenerator(0.0, 	1.0, 	-3.0,	3.0,	self.values[T3]),
			WeightedRandomGenerator(0.0, 	0.15, 	-0.75,	0.75,	self.values[T1]),
			WeightedRandomGenerator(0.0, 	0.1, 	0.,		1.0,	self.values[ALPHA]),
			WeightedRandomGenerator(0.0, 	0.0025,	0.0, 	0.05,	self.values[C_DELAY]),
			WeightedRandomGenerator(0.0, 	0.05, 	0.0,	1.0,	self.values[BETA]),
			WeightedRandomGenerator(0.0, 	0.5, 	1.0,	10.,	self.values[D_MULT]),
			WeightedRandomGenerator(0.0, 	0.05, 	0.0,	1.0,	self.values[GAMMA]),
			WeightedRandomGenerator(0.0, 	2.5, 	0.0,	50.,	self.values[MS_BINS])
		]
		self.age = 0
		self.edits = 0
	
	def __repr__(self):
		return "%9i/%9i || %.3f|%.3f|%.3f|%.3f|%.3f|%.3f|%.3f" % (self.age, self.edits, self.transp, self.values[ALPHA], self.values[C_DELAY], self.values[BETA], self.values[D_MULT], self.values[GAMMA], self.values[MS_BINS])

	def mutate(self):
		choice = random.randint(0,8)

		self.values[choice] = self.generators[choice].next()

		if choice < 3:
			self.transp = round(self.values[0])*12.0 + round(self.values[1]) + self.values[2]
			self.tratio = midi2hz(self.transp)
	

if __name__=='__main__':
	genex = GenomicExplorer('/Users/kfl/dev/python/sc-0.3.1/genomic', 'test.wav')
	genex.analyze_genome(1)
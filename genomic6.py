import sc, random, contextlib, wave, os, math
import shlex, subprocess, signal
import NRTOSCParser3

import numpy as np
import scipy.signal
import matplotlib.pyplot as plt

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
			self.val = random.randint(0,256)
		
	def next(self, scale=1.0):
		self.val = random.randint(0,256)

	def __call__(self): return self.next()

# slot assignments for sigmaSynth
ALPHA = 0
C_DELAY = 1
C_DECAY = 2
BETA = 3
D_MULT = 4
GAMMA = 5
MS_BINS = 6
DELTA = 7

class GenomicExplorer:

	def __init__(self, anchor, sfilenames, subdir=None, out_dir='out', size=50, margin=10, report_interval=20, mut_prob=0.01, stop_slope=0.1, stop_maxgens=5, rev_flag=False): #, start_state=[1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0]
				
		self.anchor = anchor
		if subdir is not None:
			self.sfpaths = [os.path.join(anchor, 'snd', subdir, sfilename) for sfilename in sfilenames]
		else:
			self.sfpaths = [os.path.join(anchor, 'snd', sfilename) for sfilename in sfilenames]
		self.out_dir = out_dir
		self.filenames = sfilenames
		self.sfinfos = []
		
		for path in self.sfpaths:
			with contextlib.closing(wave.open(path,'r')) as f: 
				frames = f.getnframes()
				rate = f.getframerate()
				dur = frames/float(rate)
				self.sfinfos += [{'rate':rate,'dur':dur}]
		
		self.mutation_prob = mut_prob
		self.depth = margin
		self.rev_flag = rev_flag
 		# 'alpha', 'c_delay', 'c_decay', 'beta', 'd_mult', 'gamma', 'ms_bins'
		self.parser = NRTOSCParser3.NRTOSCParser3(anchor=self.anchor)
		self.rawtable, self.rawmaps, self.dists, self.pool_means, self.pool_stdevs = dict(), dict(), dict(), dict(), dict()

		self.reporting_interval = report_interval
		
		self.stopping_slope = stop_slope
		self.running_avg_mean_stdevs = dict()
		self.stopping_crit_min_gens = stop_maxgens
		self.init_population(size=size)

	
	def init_population(self, size):
		self.population = []
		for n in range(size):
 			self.population += [Genome()] #random seed
		self.population[0] = Genome(values=[0,0,0,0,0,0,0,0])
		self.analyze_individual(0)
		self.activate_raw_data(0)
#		print self.population[0]
		self.compare_all_individuals(aflag=True)
			
	def mutate_pop(self):
		for indiv in range(1, len(self.population)):
			if random.random() < self.mutation_prob:
				print "indiv: ", indiv
				self.population[ indiv ].mutate()
				self.do_update_cascade(indiv)
				
	# This is performed once per mutation
	# - only affects the individual being mutated
	def do_update_cascade(self, index, clearedits=False):
		if clearedits is True:
			self.population[ index ].edits = 0
		else:
			self.population[ index ].edits += 1
		self.analyze_individual( index )
		self.activate_raw_data( index )
# 		self.compare_individual_chi_squared( index )
		self.compare_individual( index )
	
	def mate(self, a, b, kill_index):
		
		offspring = None		
		
		if random.random() < 0.5:
			offspring = self.population[a].values[:]
			
		else:
			offspring = self.population[b].values[:]
		
		# basic random gene selection from 2 parents
		for i in range(7):
			if random.random() < 0.5:
				offspring[i] = self.population[a].values[i]
			else:
				offspring[i] = self.population[b].values[i]
		# replace the killed individual with our new individual
		self.population[kill_index] = Genome(offspring)
		self.do_update_cascade(kill_index, True)
	
	def sort_by_distances(self, depth, rev=False):
		sorted_dists = [[k, self.dists[k], self.population[k].age, self.population[k].edits] for k in sorted(self.dists.keys())]
		sorted_dists = sorted(sorted_dists[1:], key = lambda row: row[1], reverse=rev) # + (maxedits - row[3])))
		print 'sorted dists: '
		print sorted_dists
		return sorted_dists[:depth], sorted_dists[(-1*depth):]
		
	def reproduce(self, depth=25):
		
		#kills, duplicates = self.sort_by_distances(depth)
		duplicates, kills  = self.sort_by_distances(depth, self.rev_flag)
		for dup in duplicates:
			self.render_individual(dup[0], 'gen'+str(self.population[0].age))
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
			if (iter%self.reporting_interval)==0:
				print self.population[0].age
				self.reproduce(self.depth)
				self.collect_population_data()
				res = self.check_for_stopping_conditions()
				if (res == 1) and (self.population[0].age > self.reporting_interval):
					return
	
	def print_all_individuals(self):
		print '== pop ==========================='
		for g in self.population: print g
	
	def start_sc(self):
		try:
			sc.start(verbose=1, spew=1, startscsynth=1)
		except OSError: # in case we've already started the synth
			print 'QUIT!'
			sc.quit()
		print 'sfpaths: ', self.sfpaths
		for i, sfpath in enumerate(self.sfpaths):
			bnum = sc.loadSnd(os.path.basename(sfpath), wait=False)
			print 'bnum: ', bnum
			self.sfinfos[i]['bnum'] = bnum
		return 1
	

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

		print '===================\n', self.sfinfos[slot]['dur']
		
		# |outbus=20, srcbufNum, start=0.0, dur=1.0, transp=1.0, c_delay=0.0, c_decay=0.0, d_mult=1.0, d_amp=0.7, ms_bins=0, alpha=1, beta=1, gamma=1|
		sc.Synth('sigmaSynth', 
		args=[
			'srcbufNum', self.sfinfos[slot]['bnum'], 
			'start', 0,
			'dur', self.sfinfos[slot]['dur']*1000,
			'transp', tr,
			'c_delay', cdelay,
			'c_decay', decay,
			'd_mult', vals[D_MULT],
			'ms_bins', vals[MS_BINS],
			'alpha', vals[ALPHA],
			'beta', vals[BETA],
			'gamma', vals[GAMMA],
			'delta', vals[DELTA]])
	
	def analyze_individual(self, index):

		print "%%%%%%%%%%%%%%%%%%"
		print "INDEX: ", index
		print len(self.sfpaths)

		if index == 0:
			oscpath = os.path.join(self.anchor, 'snd', 'osc', `index`, (os.path.splitext(self.filenames[0])[0] + '_sigmaAnalyzer2.osc'))
# 			mdpath = os.path.join(self.anchor, 'snd', 'md', `index`, self.filenames[0])
			mdpath = os.path.join(self.anchor, 'snd', 'md', `index`, (os.path.splitext(self.filenames[0])[0] + '.md.wav'))
		else:
			oscpath = os.path.join(self.anchor, 'snd', 'osc', `index`, (os.path.splitext(self.filenames[1])[0] + '_sigmaAnalyzer2.osc'))
 			mdpath = os.path.join(self.anchor, 'snd', 'md', `index`, self.filenames[0])
			mdpath = os.path.join(self.anchor, 'snd', 'md', `index`, (os.path.splitext(self.filenames[1])[0] + '.md.wav'))
		
		print "-----------------------------"
		print oscpath
		print mdpath
		
		vals = self.population[index].realvalues
		if vals[C_DELAY] < 0.01:
			cdelay = 0.0
		else:
			cdelay = vals[C_DELAY]
# 		decay = 0.9
		
		tr = self.population[index].tratio

		if index == 0:
			slot = 0
		else:
			slot = 1
		
		print (self.sfpaths[slot], index, tr, self.sfinfos[slot]['rate'], self.sfinfos[slot]['dur'])
		print ''
		print ['c_delay', cdelay, 'c_decay', vals[C_DECAY], 'd_mult', vals[D_MULT], 'ms_bins', vals[MS_BINS], 'alpha', vals[ALPHA], 'beta', vals[BETA], 'gamma', vals[GAMMA], 'delta', vals[DELTA]]
		print ''
		oscpath, mdpath = self.parser.createNRTScore(self.sfpaths[slot],
							index=index, 
							tratio=tr,
							srate=self.sfinfos[slot]['rate'],
							duration=self.sfinfos[slot]['dur'],
							params=[
								'c_delay', cdelay,
								'c_decay', vals[C_DECAY],
								'd_mult', vals[D_MULT],
								'ms_bins', vals[MS_BINS],
								'alpha', vals[ALPHA],
								'beta', vals[BETA],
								'gamma', vals[GAMMA],
								'delta', vals[DELTA]])

		cmd = 'scsynth -N ' + oscpath + ' _ _ 44100 AIFF int16 -o 1'
		print cmd
		args = shlex.split(cmd)
		p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE) #, shell=True, close_fds=True)
		
		print 'PID: ', p.pid
		rc = p.wait()
				
		print 'RC: ', rc
		if rc == 1:
			num_frames = int(math.ceil(self.sfinfos[slot]['dur'] / 0.04 / tr))
# 			print 'num frames: ', num_frames
			self.rawtable[index] = (mdpath, num_frames)
		
# 		print self.rawtable
	
	def render_individual(self, index, generation_subdir='gen0'):
	
		vals = self.population[index].realvalues
		if vals[C_DELAY] < 0.01:
			cdelay = 0.0
		else:
			cdelay = vals[C_DELAY]
# 		decay = 0.9
		
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
								'c_decay', vals[C_DECAY],
								'd_mult', vals[D_MULT],
								'ms_bins', vals[MS_BINS],
								'alpha', vals[ALPHA],
								'beta', vals[BETA],
								'gamma', vals[GAMMA],
								'delta', vals[DELTA]])

		if os.path.exists(os.path.join(self.anchor, 'snd', self.out_dir, str(generation_subdir))) is False:
			os.mkdir(os.path.join(self.anchor, 'snd', self.out_dir, str(generation_subdir)))

		cmd = 'scsynth -N ' + oscpath + ' _ ' + os.path.join(self.anchor, 'snd', self.out_dir, generation_subdir, (str(index) + '.aiff')) + ' 44100 AIFF int16 -o 1'
		args = shlex.split(cmd)
		p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE) #, shell=True, close_fds=True)
		rc = p.wait()
		if rc == 1:
		 	print 'SUCCESS: ', os.path.join(self.anchor, 'snd', self.out_dir, (str(index) + '.aiff'))
		 	rc = 0
		else:
			return None
		# cannot get this to work:
		# cmd = 'sox -b 16 ' + os.path.join(self.anchor, 'snd', self.out_dir, str(generation_subdir), (str(index) + '.aiff')) + ' ' + os.path.join(self.anchor, 'snd', self.out_dir, str(generation_subdir), (str(index) + '.wav')) + '; rm ' + os.path.join(self.anchor, 'snd', self.out_dir, str(generation_subdir), (str(index) + '.aiff'))
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
		print self.population
		for i in range(1, len(self.population)):
 			if aflag:
 				self.analyze_individual(i)
				self.activate_raw_data(i)
#  			self.compare_individual_chi_squared(i)
			self.compare_individual(i)
		print self.dists
		return self.dists
	"""
	COMPARE_INDIVIDUAL:
		... to individual in the slot that is stipulated by the arg zeroindex!
		-	by convention, we should usually put what we are comparing to in slot 0
	"""
	def compare_individual_resample(self, index, zeroindex=0):
		i_length = self.rawmaps[index].shape[0]
		zr0_length = self.rawmaps[zeroindex].shape[0]
		print i_length, ' | ', zr0_length

		# i1_length = self.rawmaps[index-1].shape[0] ## <--- NEIGHBOR comparison
		# print i_length, ' | ', i1_length, ' | ', zr0_length

		# based on length comparison, resample the mutated individuals so that they are same length as the zeroth individual (that does not mutate)
		# if indiv. is longer, resample indiv., take abs. diff., sum, div. by length
		if zr0_length < i_length:
			mfccs_dist_to_zero = float(np.sum(np.abs(scipy.signal.signaltools.resample(self.rawmaps[index][:,1:14], zr0_length, window='hanning') - self.rawmaps[0][:,1:14]))) / float(zr0_length)
			total_dist = (float(np.sqrt(np.sum(np.abs(self.rawmaps[index][:zr0_length,0] - self.rawmaps[0][:zr0_length,0])))) / float(zr0_length)) + mfccs_dist_to_zero
		# if zeroth indiv. is longer, resample zeroth indiv., take abs. diff., sum, div. by length, then do same comparison with "neighbor"
		elif i_length < zr0_length:
			mfccs_dist_to_zero = float(np.sum(np.abs(self.rawmaps[index][:,1:14] - scipy.signal.signaltools.resample(self.rawmaps[0][:,1:14], float(i_length), window='hanning')))) / float(i_length)
			total_dist = (float(np.sqrt(np.sum(np.abs(self.rawmaps[index][:i_length,0] - self.rawmaps[0][:i_length,0])))) / float(i_length)) + mfccs_dist_to_zero
		else:
		# otherwise, take abs. diff., sum, div. by length, then do amp 
			mfccs_dist_to_zero = float(np.sum(np.abs(self.rawmaps[index][:,1:14] - self.rawmaps[0][:,1:14]))) / float(zr0_length)
			total_dist = float(np.sqrt(np.sum(np.abs(self.rawmaps[index][:,0] - self.rawmaps[0][:,0])))) / float(zr0_length) + mfccs_dist_to_zero
		
		self.dists[index] = total_dist

	def compare_individual(self, index, zeroindex=0):
		i_length = self.rawmaps[index].shape[0]
		zr0_length = self.rawmaps[zeroindex].shape[0]
		print i_length, ' | ', zr0_length

		min_length = min(i_length, zr0_length)
		print i_length, ' | ', zr0_length, ' | ', min_length

		# based on length comparison, resample the mutated individuals so that they are same length as the zeroth individual (that does not mutate)
		# if indiv. is longer, resample indiv., take abs. diff., sum, div. by length
		mfccs_dist_to_zero = float(np.sum(np.abs( self.rawmaps[index][:zr0_length,1:14] - self.rawmaps[0][:min_length,1:14]))) / float(min_length)
		total_dist = (float(np.sqrt(np.sum(np.abs(self.rawmaps[index][:min_length,0] - self.rawmaps[0][:min_length,0])))) / float(min_length)) + mfccs_dist_to_zero
		
		self.dists[index] = (total_dist / float(min_length))
	
	def compare_individual_chi_squared(self, index):
		i_length = self.rawmaps[index].shape[0]
		i1_length = self.rawmaps[index-1].shape[0]
		zr0_length = self.rawmaps[0].shape[0]
# 		print i_length, '|', zr0_length
		# based on length comparison, resample the mutated individuals so that they are same length as the zeroth individual (that does not mutate)
		# if indiv. is longer, resample indiv., take abs. diff., sum, div. by length
		if zr0_length < i_length:
			mfccs_dist_to_zero = scipy.stats.mstats.chisquare(scipy.signal.signaltools.resample(self.rawmaps[index], zr0_length, window='hanning'), self.rawmaps[0])
# 			print self.dists[index]
		# if zeroth indiv. is longer, resample zeroth indiv., take abs. diff., sum, div. by length, then do same comparison with "neighbor"
		elif i_length < zr0_length:
			mfccs_dist_to_zero = scipy.stats.mstats.chisquare(self.rawmaps[index], scipy.signal.signaltools.resample(self.rawmaps[0], i_length, window='hanning'))
		else:
		# otherwise, take abs. diff., sum, div. by length, then do same comparison with "neighbor"
			 print 'CHI-ZERO'
			 mfccs_dist_to_zero = scipy.stats.mstats.chisquare(self.rawmaps[index], self.rawmaps[0])
		
		if i1_length < i_length:
			neighbor_dist = scipy.stats.mstats.chisquare(scipy.signal.signaltools.resample(self.rawmaps[index-1], i_length, window='hanning') - self.rawmaps[index])
		elif i_length < i1_length:
			neighbor_dist = scipy.stats.mstats.chisquare(self.rawmaps[index-1], scipy.signal.signaltools.resample(self.rawmaps[index], i1_length, window='hanning'))
		else:
			print 'CHI-NEIGHBOR'
			neighbor_dist = scipy.stats.mstats.chisquare(self.rawmaps[index-1], scipy.signal.signaltools.resample(self.rawmaps[index], i1_length, window='hanning'))
		
		nsum = np.sum(np.abs(neighbor_dist[0].data[:24]))
		zsum = np.sum(np.abs(mfccs_dist_to_zero[0].data[:24]))
		nasum = neighbor_dist[0].data[24]
		zasum = mfccs_dist_to_zero[0].data[24]
		
		self.dists[index] = nsum + zsum - (24.0 * nasum) - (24.0 * zasum)
	
	def collect_population_data(self):
		diffs = [self.dists[k] for k in self.dists.keys()]
		print 'diffs: ', diffs
		age0 = self.population[0].age
		age1 = self.population[1].age
		print 'ages: ', age0, '|', age1
		self.pool_means[age1] = np.mean(diffs)
		self.pool_stdevs[age1] = np.std(diffs)

	def collect_population_data_resample(self):
		zero_data = np.array(self.rawmaps[0][:,1:14])
		zr0_length = zero_data.shape[0]
		diffs = []
		for indiv in range(1, len(self.population)):
			data = np.array(self.rawmaps[indiv][:,1:14])
			i_length = data.shape[0]
			if (i_length > zr0_length):
				diffs += [np.sum(np.abs(scipy.signal.signaltools.resample(data, zr0_length, window='hanning') - zero_data))]
			elif (i_length < zr0_length):
				diffs += [np.sum(np.abs(data - scipy.signal.signaltools.resample(zero_data, i_length, window='hanning')))]
			else:
				diffs += [np.sum(np.abs(data - zero_data))]
			
		diffs = np.array(diffs)
		age0 = self.population[0].age
		age1 = self.population[1].age
		print 'ages: ', age0, '|', age1
		self.pool_means[age1] = np.mean(diffs)
		self.pool_stdevs[age1] = np.std(diffs)

	def check_for_stopping_conditions(self):
		"""
		0 = continue
		1 = stop
		"""
		age0 = self.population[0].age
		stdevs_skeys = sorted(self.pool_stdevs.keys())
		self.stdevs_ordered = [self.pool_stdevs[key] for key in stdevs_skeys]
		lastNstdevs = self.stdevs_ordered[(-1*self.stopping_crit_min_gens):]
		self.running_avg_mean_stdevs[age0] = mean(lastNstdevs)
		print ">>>>>>>>>>>>>>>>> STOP??? ::: ", (abs(max(lastNstdevs) - min(lastNstdevs)) /  self.stopping_crit_min_gens)
		print ">>>>>>>>>>>>>>>>> MEAN    ::: ", mean(lastNstdevs)
		if (len(lastNstdevs) < self.stopping_crit_min_gens) or ((abs(max(lastNstdevs) - min(lastNstdevs)) /  self.stopping_crit_min_gens) > self.stopping_slope):
			print " continue ..."
			return 0
		else:
			print "**STOP**"
			return 1 # signal stop
	
	def population_realvalues_as_array(self):
		realvals = []
		for indiv in self.population:
			realvals += indiv.realvalues
		return np.array(realvals).reshape((-1,7))
	
	def population_8bitvalues_as_array(self):
		vals = []
		for indiv in self.population:
			vals += indiv.values
		return np.array(vals).reshape((-1,7))
		
	


class Genome:

	def __init__(self, values=None, slotranges=[[1.0,0.5],[0.0,0.05],[0.0,1.0],[1.0, 0.5],[1.0,10.],[1.0,0.5],[0.0,50.],[1.0,0.5]]):
		
		# """
		# 'alpha',   'c_delay', 'decay'    'beta',   'd_mult', 'gamma',  'ms_bins', 'delta'
		# [[1.0,0.5],[0.0,0.05],[0.0,1.0],[1.0,0.5],[1.0,10.],[1.0,0.5],[0.0,50.]],[1.0,0.5]
		# """
		
		self.tratio	= 1.0		# CHECK THIS... WHY IS IT HERE/in Hertz!!! ???
		
		self.boundaries = slotranges
		self.generators = [RandomGenerator_8Bit(-1) for n in range(8)] ### CONSTANT WARNING
		#StaticGenerator_8Bit(VAL) ???
		
		if values is None:
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
		#print '1. ', tuple(self.values)
		#print '2. ', ((self.age, self.edits) + tuple(self.values) + tuple(self.binarystring))
		return "%9i/%9i || %.6f|%.6f|%.6f|%.6f|%.6f|%.6f|%.6f" % ((self.age, self.edits) + tuple(self.realvalues)) # + tuple(self.binarystring)
	
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

def mean(arr):
	return sum([(float(val)/len(arr)) for val in arr])

def lininterp(val,bounds=[0.,1.]):
	return (((val/128.0)*(bounds[1]-bounds[0]))+bounds[0])

def substitute_char_in_string(s, p, c):
	l = list(s)
	l[p] = str(c)
	return "".join(l)

# conversion function
def midi2hz(m): return pow(2.0, (m/12.0))
	
def vals_to_binarystring(vals = [0, 0, 0, 0, 0]):
	return ''.join((("{0:08b}".format(val)) for val in vals))

# never a '0bXXX' string!
def binarystring_to_vals(binstring):
	mystring = binstring[:]
	length = len(mystring) / 8 # ignore the last digits if it doesn't chunk into 8-item substrings
	res = []
# 	print mystring[(n*8):((n+1)*8)]
	return [int(mystring[(n*8):((n+1)*8)], 2) for n in range(length)]

def plot_one_generational_means_stdevs(pool_means1, pool_stdevs1, poolsize1):
	
	fig, ax1 = plt.subplots(nrows=1, ncols=1)
	
	timepoints1 = sorted(pool_means1.keys())
	means1 = np.array([pool_means1[tp] for tp in timepoints1])
	stdevs1 = np.array([pool_stdevs1[tp] for tp in timepoints1])
	lower_stdevs1 = np.where(np.subtract(means1, (stdevs1/2))>0, stdevs1/2, means1) # 0
	ax1.errorbar(timepoints1, means1, yerr=[lower_stdevs1, stdevs1/2], fmt='o')
		
	ax1.set_xlabel('Number of generations')
	ax1.set_ylabel('Fitness score/dissimilarity')
	plt.show()


def plot_generational_means_stdevs(pool_means1, pool_stdevs1, poolsize1, pool_means2, pool_stdevs2, poolsize2):
	
	fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharey=True)
	
	timepoints1 = sorted(pool_means1.keys())
	means1 = np.array([pool_means1[tp] for tp in timepoints1])
	stdevs1 = np.array([pool_stdevs1[tp] for tp in timepoints1])
	lower_stdevs1 = np.where(np.subtract(means1, stdevs1)>0, stdevs1, means1) # 0
	ax1.errorbar(timepoints1, means1, yerr=[lower_stdevs1/2, stdevs1/2], fmt='o')
	#
	timepoints2 = sorted(pool_means2.keys())
	means2 = np.array([pool_means2[tp] for tp in timepoints2])
	stdevs2 = np.array([pool_stdevs2[tp] for tp in timepoints2])
	lower_stdevs2 = np.where(np.subtract(means2, stdevs2)>0, stdevs2, means2) # 0
	ax2.errorbar(timepoints2, means2, yerr=[lower_stdevs2/2, stdevs2/2], fmt='o')
	
	ax1.set_xlabel('Number of generations')
	ax2.set_xlabel('Number of generations')
	ax1.set_ylabel('Fitness score/dissimilarity')
	plt.show()
	

# if __name__=='__main__':
# 	genex = GenomicExplorer('/Users/kfl/dev/python/sc-0.3.1/genomic', 'test.wav')
# 	genex.analyze_genome(1)
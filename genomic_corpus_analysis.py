import os, glob

import corpusdb
from genomic5 import *


#	create corpus
anchorpath = os.path.expanduser('~/dev/git/Genomic')
genomic_corpus = corpusdb.CorpusDB(anchorpath)

sndpath = os.path.join(anchorpath, 'snd')
sfid_filename_map = dict()

# iterate on all files in snd dir
allfiles = glob.glob( os.path.join(sndpath, 'combos/*.wav') )

for a, sndA in enumerate(allfiles):
	# for A in each of 10 sounds:
	
	# add untreated sound to corpus
	parent_node = genomic_corpus.add_sound_file(filename=os.path.basename(sndA), subdir='combos',)
	parent_sfid = parent_node.sfid
	
	sfid_filename_map[parent_sfid] = os.path.basename(sndA)
	sfid_filename_map[os.path.basename(sndA)] = parent_sfid
	
	#analyze and add file-length unit
	genomic_corpus.analyze_sound_file(os.path.basename(sndA), parent_sfid)
	
	sfdur = genomic_corpus.sftree.nodes[parent_sfid].duration
	
	# sfid and tag *should* match for parents...!
	parent_sfid = genomic_corpus.add_sound_file_unit(parent_sfid, onset=0, dur=min(0.5, sfdur), tag=sfid_filename_map[os.path.basename(sndA)])
	
	genomic_corpus.segment_units(parent_sfid)
	genomic_corpus._deactivate_raw_sound_file(parent_sfid)
	
# we now have 10 parent nodes; the rest will be child nodes
for a, sndA in enumerate(allfiles):

	# for B in each of the other 9 sounds, do 3 times:
	for sndB in sorted(list(set(allfiles).difference(set([allfiles[a]])))):
		print '    ', sndB
		# create gexex with A,B
		genex = GenomicExplorer('/Users/kfl/dev/git/Genomic', [os.path.basename(sndA), os.path.basename(sndB)], subdir='combos',size=50, margin=10, report_interval=5, mut_prob=0.02, stop_slope=0.001)
		
		for gens in range(100):
			# iterate until stop
			genex.iterate(20)
		
			rand5 = range(len(genex.population))
			random.shuffle(rand5)
		
			# sample 5 indivs, add each to corpus
			for slot in rand5[:5]:
				curragent = genex.population[slot]
				print 'curragent: ', curragent
			
				# We now have:
				#	list of 2 snd files
				#	list of 7 params
			
				genomic_params = curragent.realvalues
				print '>>>> ', genomic_params
				gparams = ['alpha', genomic_params[0],
					'c_delay', genomic_params[1],
					'c_decay', genomic_params[2],
					'beta', genomic_params[3],
					'd_mult', (genomic_params[4] + random.random()),
					'gamma', genomic_params[5],
					'ms_bins', genomic_params[6]]
			
				print 'src. file id: ', sfid_filename_map[ os.path.basename(sndB) ]
				child_node = genomic_corpus.add_sound_file(filename=None, 
															sfid=None, 
															srcFileID=sfid_filename_map[ os.path.basename(sndB) ],
															subdir='combos',
															synthdef=['sigmaMod'], 
															params=[(['inbus', 1000, 'outbus', 1000] + gparams)])
				child_sfid = child_node.sfid
				print 'child_sfid: ', child_sfid
				#analyze and add file-length unit
				genomic_corpus.analyze_sound_file(os.path.basename(sndB), child_sfid)

				sfdur = genomic_corpus.sftree.nodes[child_sfid].duration

				child_sfid = genomic_corpus.add_sound_file_unit(child_sfid, onset=0, dur=min(0.5, sfdur), tag=sfid_filename_map[os.path.basename(sndB)])

				genomic_corpus.segment_units(child_sfid)
				genomic_corpus._deactivate_raw_sound_file(child_sfid)
			
#	save corpus
corpus_array = genomic_corpus.convert_corpus_to_array(map_flag=True)
print corpus_array.shape
genomic_corpus.export_corpus_to_json('violas.json', 'violas.corpusdata')
#	make music







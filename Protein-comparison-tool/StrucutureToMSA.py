#run fatcat with structures in a folder
#condense alignment into MSA
#it requires Fatcat and clustalo to work
#fatcat is the java based version included with the protein comparison tool
#from the PDB website.
#the phylogeny is derived from a distance matrix

import glob
import subprocess
import itertools
import shlex
from Bio import AlignIO , SeqIO
from Bio.PDB import *

import pickle
import sys
import re
import numpy as np
import pandas as pd
import math
#import pylab
from itertools import combinations, permutations
import multiprocessing as MP
import os
import random
import hashlib
import tempfile as tmp



random.seed(0)

overwrite = False
TMoverwrite = False

#find input pdbs
grabPDBs = True

#steps in the process of alignment and outputing an MSA and tree
align = True
#compute all pairwise alignments with fatcat
Fatcatalign = True
#realign flexed structures with TM align
TMalign_from_Fatcat = True
gesamt_from_Fatcat = False
#use fatcat output to generate pairwise fastas
parse = False
#use TMalign results on flexed fatcat structures to generate pairwise fastas
maketmfastas = False
#use clustalo to create MSA of the strucutres
merge = False

#output a distance matrices
make_distmat = True
make_TMdistmat = True
SDM = False
show_distmat = False
splitFatcat = False

#where is fatcat?
FatcatPath =  './runFATCAT.sh'
TMalignPath = './tmalign'
Dalipath = './DaliLite3.3/DaliLite '
clustaloPath = 'clustalo'
fastmepath = 'fastme '
gesamtpath = '/scratch/dmoi/software/ccp4/ccp4-7.1/bin/gesamt'

#split directory
splitpath = 'split/'

def save_obj(obj, name ):
	with open( name + '.pkl', 'wb') as f:
		pickle.dump(obj, f)

def load_obj(name ):
	with open( name + '.pkl', 'rb') as f:
		return pickle.loads(f.read())



def runFastme( fastmepath , clusterfile ):
		args =  fastmepath +  ' -i ' + clusterfile + ' -o ' + clusterfile+'_tree.txt'
		print( args)
		p = subprocess.call(shlex.split(args) , stdout=subprocess.PIPE )
		return p,[clusterfile+'_tree.txt' ]

def distmat_to_txt( pdblist , distmat, filedir , name):
	#write out distmat in phylip compatible format
	print(distmat)
	outstr = str(len(pdblist)) + '\n'
	for i,pdb in enumerate(pdblist):
		if '/' in pdb:
			pdb = pdb.split('/')[-1]

		namestr = pdb.replace('.','').replace('_','')[0:20]
		outstr += namestr+ ' ' +''.join([ "{:.2f}".format(n) + ' ' for n in list(distmat[i,:])])  + '\n'
	print( outstr)
	handle = open(filedir + name + 'fastmemat.txt' , 'w')
	handle.write(outstr)
	handle.close()
	return filedir + name + 'fastmemat.txt'

	
def check_pair(pairvec):
	occurences = {}
	if pairvec[0] == pairvec[1]:
		return False
	lists = [[],[]]
	for i,fasta in enumerate(pairvec):
		with open(fasta, 'r') as lines:
			for line in lines:
				if '>' in line :
					ID = line.split('.')[0]
					lists[i].append(ID)
	commonIDs = list(set(lists[0]) & set(lists[1]))
	if len(commonIDs) > 0 :
		return True
	else :
		return False

def check_identical(fasta):
	seq= ['','']
	i = 0
	try:
		if fasta:
			record_iterator = iter(AlignIO.read(fasta, "fasta"))
			first = next(record_iterator)
			second = next(record_iterator)
		else:
			return None
	except:
		print('read error')
		print( fasta )

	if str(first.seq).replace('\n','') not in str(second.seq).replace('\n',''):
		if '-' in str(first.seq)+str(second.seq):

			return fasta
	else:
		print( fasta )
		return None


def check_ids(fasta):
	record_iterator = iter(AlignIO.read(fasta, "fasta"))
	return [ rec.id for rec in record_iterator]



def runFatcat(pargs):

	FatcatPath , struct1, struct2, outfile , outPDB = pargs

	'''if os.path.isfile(outfile):
					with open( outfile , 'r')as preout:
						if len(preout.read())>0:
							print(outfile + ' already exists	')
							return outfile'''
	args =  FatcatPath + ' -file1 ' + struct1 + ' -file2 ' + struct2 + '  -printFatCat -flexible True -outputPDB  -outFile '+ outPDB
	print( args)
	args = shlex.split(args)
	fatcatfile = open(outfile, 'w')
	p = subprocess.run(args  , stdout= fatcatfile )
	try:
		p.wait()
	except:
		pass
	fatcatfile.close()
	return outfile , [ struct1 , struct2 ]

def runTMalign(pargs):
	FatcatPath , struct1, struct2, outfile , outPDB = pargs
	args =  TMalignPath + ' ' + struct1 + ' ' + struct2
	print( args)
	args = shlex.split(args)
	TMfile = open(outfile, 'w')
	p = subprocess.call(args , stdout= TMfile )

	return outfile


def splitFile(Fatcatstruct, outdir , pdbs):
	"""
	TITLE  jFatCat_flexible 1.0 domain3_4hj1.pdb vs. domain3_4adi.pdb
	EXPDTA    NMR, 2 STRUCTURES
	MODEL      1
	"""

	coords = {}
	header = ''
	i =0
	start = False
	with open( Fatcatstruct , 'r' ) as PDBfile:
		for line in PDBfile:
			if start == False:
				if 'EXPDTA' in line:
					header += line.replace('2','1')
				else:
					header += line
			if 'TITLE' in line:
				pdb2 = line.replace('TITLE  jFatCat_flexible 1.0' , '' ).split('vs.')[0].strip()
				pdb1 = line.replace('TITLE  jFatCat_flexible 1.0' , '' ).split('vs.')[1].strip()
			if 'MODEL' in line:
				start = True
				i+=1
				coords[i] = ''
			if start == True:
				coords [i] += line

	names = []
	for pdb in pdbs:
		names.append(pdb.split('/')[-1])

	filehandle1 = open( outdir + names[0].replace('pdb', '') + 'vs' + names[1] , 'w')
	filehandle2 = open( outdir + names[1].replace('pdb', '') + 'vs' + names[0] , 'w')
	handles = [filehandle1, filehandle2]
	pdbs = [pdb1, pdb2]

	for i,struct in enumerate(coords):
		handles[i].write(header)
		handles[i].write(coords[struct])
		handles[i].write('ENDMDL')

	for handle in handles:
		handle.close()

	return 0;


def cleanPDBS(pdbs):
	print(pdbs)
	for file in pdbs:
		output = ''
		with open( file , 'r' ) as filein:
			for l in filein:
				if 'SSBOND' in l:
					pass
				else:
					output += l
		with open( file , 'w' ) as fileout:
			fileout.write(output)


def runTMalign_FatcatOut(pargs, overwrite = True):
	TMalignPath , Fatcatstruct, outfile , pdb_pair = pargs
	TMfile = open(outfile, 'w')
	"""
	TITLE  jFatCat_flexible 1.0 domain3_4hj1.pdb vs. domain3_4adi.pdb
	EXPDTA    NMR, 2 STRUCTURES
	MODEL      1
	"""
	coords = {}
	header = ''
	i =0
	start = False

	if os.path.isfile(outfile) or overwrite == True:
		b = os.path.getsize(outfile)
		if b>0:
			return None


	with open( Fatcatstruct , 'r' ) as PDBfile:
		for line in PDBfile:
			if start == False:
				if 'EXPDTA' in line:
					header += line.replace('2','1')
				else:
					header += line
			if 'TITLE' in line:
				pdb2 = line.replace('TITLE  jFatCat_flexible 1.0' , '' ).split('vs.')[0].strip()
				pdb1 = line.replace('TITLE  jFatCat_flexible 1.0' , '' ).split('vs.')[1].strip()
			if 'MODEL' in line:
				start = True
				i+=1
				coords[i] = ''
			if start == True:
				coords [i] += line
	pdbs = [pdb1, pdb2]


	
	try:
		os.mkdir('tempfiles/')
	except:
		pass

	with open('tempfiles/' + pdb1+'vs'+ pdb2 , 'w') as handle:
		handle.write(header)
		handle.write(coords[1])
		handle.write('ENDMDL')
	
	with open( 'tempfiles/' + pdb2+'vs'+ pdb1 , 'w') as handle:
		handle.write(header)
		handle.write(coords[2])
		handle.write('ENDMDL')

	args =  TMalignPath + ' ' + 'tempfiles/' +pdb1+'vs'+pdb2 + ' ' + 'tempfiles/' +pdb2+'vs'+pdb1
	print(args)
	args = shlex.split(args)
	TMfile = open(outfile, 'w')
	p = subprocess.call(args , stdout= TMfile )
	TMfile.close()
	os.remove('tempfiles/' + pdb1+'vs'+pdb2)
	os.remove('tempfiles/'+ pdb2+'vs'+pdb1)
	return outfile , pdb_pair

def mergeAlign(clustalopath , fasta1, fasta2, outfile):
	args =  clustalopath + ' --p1 ' + fasta1 + ' --p2 ' + fasta2 + ' --force -o ' + outfile
	print( args)
	args = shlex.split(args)
	p = subprocess.call(args )

def FatcatToFasta(fatout, filename):

	alns = ['','']
	files = []
	BlockNext = False
	Blocknum = 0
	with open( fatout , 'r' ) as fatfile:
		for line in fatfile:
			try:
				#print( line)
				if 'file from local' in line:
					files.append(line.split('/')[-1])
				#format  = Chain 1:   16 PHC---------------SKTPIVRAQTSQNAMS-------RGMQMQFSIGLHT---------AVC----
				if 'Chain 1:' in line:
					alns[0] += line.split(':')[1].strip().split(' ')[1]
				if 'Chain 2:' in line:
					alns[1] += line.split(':')[1].strip().split(' ')[1]
			except:
				print('line error:')

		try:
			alnstr1 = '>'+files[0][:-1] + 'Fatcat_block' +'\n' +alns[0] + '\n'
			alnstr2 = '>'+files[1][:-1] + 'Fatcat_block' +'\n' +alns[1] + '\n'
			print( alnstr1 + '\n' +alnstr2 )
			handle = open(filename  , 'w')
			handle.write(alnstr1 + '\n')
			handle.write(alnstr2)
			handle.close()
			return filename
		except:
			print([files , alns ])

def FatcatToDistances(fatout):
	alns = ['','']
	files = []
	BlockNext = False
	Blocknum = 0
	pval = 0
	pdbs = []
	cleanlen = 0

	distances = {}
	counts = {}
	numbers = re.compile('[-+]?[0-9]*\.?[0-9]+.')
	integers = re.compile('\s[0-9]+\s')

	with open( fatout , 'r' ) as fatfile:
		for i,line in enumerate(fatfile):
			if 'P-value' in line:
				pval = float(line.split(' ')[1])
				#P-value 1.44e-03 Afp-num 65162 Identity 4.96% Similarity 16.08%
			if 'Align' in line:
				#Align 1SVB.pdb 395 with ememmodelFinalarabidop_chopped.pdb.pdb 467
				files = line.split(' ')
				pdbs = []
				for filename in files:
					if '.pdb' in filename.lower():
						pdbs.append(filename)
				lengths = integers.findall(line)
				cleanlen = []
				for val in lengths:
					cleanlen.append(int(val.strip()))
			if 'Block' in line:
				#Block  0 afp 17 score 201.70 rmsd  3.92 gap 128 (0.48%)
				words = line.split('rmsd')[1]
				rmsd  = float(numbers.findall(words)[0].strip())
				distances[Blocknum] = rmsd
				Blocknum+=1
			if 'Block' not in line and 'Chain' not in line and i > 6:
				for j in range(len(distances)):
					if j not in counts:
						counts[j] = 0
					counts[j]+= line.count(str(j+1))
	return pdbs, cleanlen, distances, counts , pval




def TMaligntoDistance(args , verbose = False):
	#parse TM align file and output parameters
	TMalignfile = args
	"""
	format:
		 **************************************************************************
	 *                        TM-align (Version 20160521)                     *
	 * An algorithm for protein structure alignment and comparison            *
	 * Based on statistics:                                                   *
	 *       0.0 < TM-score < 0.30, random structural similarity              *
	 *       0.5 < TM-score < 1.00, in about the same fold                    *
	 * Reference: Y Zhang and J Skolnick, Nucl Acids Res 33, 2302-9 (2005)    *
	 * Please email your comments and suggestions to: zhng@umich.edu          *
	 **************************************************************************

	Name of Chain_1: ../structural_alignment/effHapViralcontrol/1SVB.pd
	Name of Chain_2: ../structural_alignment/effHapViralcontrol/2ala.pd
	Length of Chain_1:  395 residues
	Length of Chain_2:  384 residues

	Aligned length=  325, RMSD=   5.63, Seq_ID=n_identical/n_aligned= 0.062
	TM-score= 0.56731 (if normalized by length of Chain_1)
	TM-score= 0.57996 (if normalized by length of Chain_2)
	(You should use TM-score normalized by length of the reference protein)

	(":" denotes aligned residue pairs of d < 5.0 A, "." denotes other aligned residues)
	SRCTHLENRDFVTGTQG-TTRVTLVLELGGCVTITAEG-----K-PSMDVW-LDAIYQENPAKTREYCLHAKLSDTKVAARCPTMGPATLAEEHQGGTVCKRDQ-SDRGW-G-N-HC-GLFGKGSIVACVKAAC--EAKKKATGHVYDANKIVYTVKVEPHTGDYVAANET-HSGRKTASFTISSEKTILTMGEYGDVSLLCRVASGVDLAQTVILELDKTVEHLPTAWQVHRDWF-ND----LALPWKHE-GAQNWNNAERLVEFGAPHA--VKMDV-YNLGDQTGVLLKALAGVPV-------AHIEGTKYHLK-SGHVTCEVGLEK---LKMKG--LTYT-MCDKTKFTWKRAPTDSGHDTVVMEVTFS-GT-----KPCRIPVRAVAHGSPDVNVAMLITPNPTIENN-----GGGFIEMQLPPGDNIIY-V-------GELSHQWFQK-
					  .::::. .......::::::     : :::::: .::::.. :::: ::::::::: ::::::..::::.. . .::::.::::.: .:::: : : :: ::...:::::::.:::  ::::::::::...  .::::::::...         .:::::::..  .:::::: :  .:.::.::...:.....::::::   :  :: :::::.... ..    :::::::: .:::::: ::......:::  ::::. ..::.::.::.::.:....       ::::  ::::. .....::::..:   .....  .:.. ..... ...  .:.:.    .::....  .      ........         ...::::..::...      ...::::.:   ::... .       .....::..
	-----------------YEHSTVM-PNVVGFPYKAHIERPGYSPLTLQMQVVETSLEPT-LNLE-YITCEYKTV-VPSPYVKCCGASEC-S-TKEKPDYQCKVYTGVYPFMWGGAYCFCDSENTQLSEAYVDRSDVCRHDHASAYKAHT--ASLKAKVRVMYG--------NVNQTVDVYVN--GDHAVTI-G--GTQFIFGPLSSAWTPFDNKIVVY---K--DE-VFNQDFPPYGSGQPGRFGDIQSRTVESNDLYA-NTALKLARPSPGMVHVPYTQTPSGFKYWLKEKGTALNTKAPFGCQIKTN--PVRAMNCAVGNIPVSMNLPDSAFTRIVEAPTIIDLTCT-VAT--CTHSS----DFGGVLT-LT-YKTNKNGDCSVHS---------HSNVATLQEATAKV-KTAGKVTLHFSTAS---ASPSFVVSLCSARATCSASCEPP-K

	"""
	#maybe incorporate more stuff into this later

	TMscores =[]
	readpdbs =[]
	line1 = False
	line2 = 10**10
	with open( TMalignfile, 'r') as output:
		i = 0
		for i,line in enumerate(output):
			if verbose == True:
				print(line)
			if 'Name of Chain_1:' in line:
				pdb1 = line.split('/')[-1].split('.')[0]
			if 'Name of Chain_2:' in line:
				pdb2 = line.split('/')[-1].split('.')[0]
			if 'Length of Chain_1:' in line:
				chain1len = int(line.split(':')[1].split()[0] )
			if 'Length of Chain_2:' in line:
				chain2len = int(line.split(':')[1].split()[0]	)
			if 'Name of Chain_' in line and '*' not in line :
				readpdbs.append(line.split(':')[1].strip())
			if 'TM-score' in line and '*' not in line and ' normalized by length of the reference structure' not in line:
				score = float(line.split('=')[1].split( )[0])
				TMscores.append(score)


			if line1== True:
				fasta1= line
				line2 = i+2
				line1 = False

			if '(":" denotes residue pairs of d <  5.0 Angstrom, "." denotes other aligned residues)' in line:
				line1 = True

			if i == line2:
				fasta2 = line
		else:
			if i ==0:
				print('err' , TMalignfile )
				return 0 , [None, None] , None 
			else:
				outstr = ''.join(['>' ,pdb1 ,'\n' , fasta2 , '\n' ,  '>' ,pdb2 ,'\n' , fasta1 ,'\n' ])
				return max(TMscores[0], TMscores[1]), [pdb1,pdb2] , outstr
				


def FatcatToDF(fatout):
	print( fatout)
	alns = ['','']
	files = []
	pdbfile =fatout+ 'PDB'
	BlockNext = False
	Blocknum = 0
	pval = 0
	df1 = None
	df2 = None
	chain1_read = False
	chain2_read = False

	distances = {}
	counts = {}
	df1dict = {}
	df2dict = {}
	numbers = re.compile('[-+]?[0-9]*\.?[0-9]+.')
	integers = re.compile('\s[0-9]+\s')
	number = re.compile('[0-9]+')
	sequence = re.compile('[A-Z-]+')

	with open( fatout , 'r' ) as fatfile:
		for i,line in enumerate(fatfile):
			if 'Align' in line:
				#Align 1SVB.pdb 395 with ememmodelFinalarabidop_chopped.pdb.pdb 467
				files = line.split(' ')
				pdbs = []
				for filename in files:
					if '.pdb' in filename.lower():
						pdbs.append(filename)
				lengths = integers.findall(line)
				cleanlen = []
				for val in lengths:
					cleanlen.append(int(val.strip()))
			if 'Block' in line:
				#Block  0 afp 17 score 201.70 rmsd  3.92 gap 128 (0.48%)
				words = line.split('rmsd')[1]
				rmsd  = float(numbers.findall(words)[0].strip())
				distances[Blocknum] = rmsd
				Blocknum+=1

			#Chain 1:   19 TWVDLVLEGDSCVTIMSK----DKPTIDVKMMNMEAANLAEVRSYCYLATVSDLSTKAACPTMGEAHNDK
				#  			   111111111111111111    111111111111111 111111111111111111         11111
			#Chain 2:    3 HVTVIPNTVGVPYKTLVNRPGYSPMVLEMELLSVTLEPTLSLDYITCEYKTVIPSPYVKCCGTAECKDKN

			if chain1_read == True:
				#read block asiignment
				chain1_read = False
				blocks = line[14:]

			if 'Chain 1' in line:
				chain1_read = True
				cleanline = line.split(':')[1].strip()
				start1 = int(number.findall(cleanline)[0])
				alnseq1 = sequence.findall(cleanline)[0]

			if 'Chain 2' in line:
				chain2_read = True
				cleanline = line.split(':')[1].strip()
				start2 = int(number.findall(cleanline)[0])
				alnseq2 = sequence.findall(cleanline)[0]

			if chain2_read == True:
				chain2_read = False
				#parse the two chains and map residues to block, gap or aln positions
				#residue count
				count1 = 0
				count2 = 0
				for i,char1 in enumerate(alnseq1):
					if i < len(alnseq2):
						Blocknum = blocks[i]
						char2 = alnseq2[i]
						skip = False
						if char1 == '-':
							# assign values for the residue that opened the gap
							row2 = [char2,'NONE', 'gap' , 'NONE' , 'NONE']
							df2dict[start2+count2] = row2
							count2 +=1
							skip = True
						if char2 == '-':
							row1 = [char1,'NONE', 'gap' , 'NONE' ,  'NONE']
							df1dict[start1+count1] = row1
							count1 +=1
							skip = True
						if skip == False:
							#spaces will be counted as aligned but with no blocknum
							row1 = [ char1, Blocknum , 'aln' , start2+count2 , char2 ]
							df1dict[start1+count1] = row1
							row2 = [ char2, Blocknum , 'aln' , start1+count1 , char1 ]
							df2dict[start2+count2] = row2
							count2 +=1
							count1 +=1
	df1 = pd.DataFrame.from_dict( df2dict , orient= 'index' )
	df2 = pd.DataFrame.from_dict( df1dict , orient= 'index' )
	df1.columns = ['res','block','status','match','matchres']
	df2.columns = ['res','block','status','match','matchres']
	return df1, df2 , pdbs, pdbfile


def structalign(filedir):
	#filedir = '../Structures/structalign_finalset/'
	outpath = filedir
	if align == True:
		structures = glob.glob(filedir + '*.pdb')
		already = glob.glob(outpath + '*')
		cleanPDBS(structures)
		results = []
		if Fatcatalign == True :
			fatcatfiles = []
			runData = []
			for struct1,struct2 in itertools.combinations(structures, 2):
				outfile = outpath + struct1.split('/')[-1].replace('.pdb','')+ '_'+struct2.split('/')[-1].replace('.pdb','') + '.factcat'
				outPDB = outfile + 'PDB'
				if overwrite == True or outPDB not in already:

					runData.append([FatcatPath ,struct1,struct2, outfile, outPDB])
				
				elif os.path.getsize(outPDB) == 0:

					runData.append([FatcatPath ,struct1,struct2, outfile, outPDB])
				else:
				
					results.append( [ outfile , [ struct1, struct2] ] )
			print(len(runData))
			print('alignments to run')
			print( 'fatcat align allvall')
			pool = MP.Pool()
			results += pool.map_async( runFatcat , runData, MP.cpu_count() ).get()
			pool.close()
			save_obj( results, filedir +'fatcatfiles')
		if TMalign_from_Fatcat == True:
			#run this after fatcat!
			runData = []
			results = []
			already = glob.glob(outpath + '*factcatTMRealign' )

			alignments =  load_obj(filedir +'fatcatfiles' )
						
			for aln,pdb_pair in alignments:
				outfile = aln+'TMRealign'
				if TMoverwrite == True or outfile not in already:
					#TMalignPath , Fatcatstruct, outfile , pdb_pair
					runData.append([TMalignPath , aln+'PDB' ,  outfile , pdb_pair ])
				elif os.path.getsize(outfile) == 0:
					runData.append([TMalignPath , aln+'PDB' ,  outfile , pdb_pair ])
				else:
					results.append( (outfile , pdb_pair) )
			
			print( len(runData), 'alignments to run')
			pool = MP.Pool()
			results += pool.map_async( runTMalign_FatcatOut, runData, MP.cpu_count() ).get()
			save_obj( results, filedir +'TMRealignfiles')
			
			pool.close()
			#for args in runData:
			#	results.append(runTMalign_FatcatOut(args) ) 
			print('done aligning')
	
	if parse == True:
		print( filedir)
		alignments = glob.glob(filedir + '*.factcat')
		print(len(alignments))
		fastas = []
		distances = {}

		print( 'converting fatcat output to fasta')
		for aln in alignments:
			fastas.append(FatcatToFasta(aln,aln+'.fasta'))

		for aln in alignments:
			distances[aln] = FatcatToDistances(aln)

		if maketmfastas == True:
			tmalignments = load_obj(filedir + 'TMRealignfiles')
			print(tmalignments)
			tmfastas = []
			for aln in tmalignments:
				print(aln)
				retval  = TMaligntoDistance(aln)
				if retval:
					distance, pdbs, outstr = retval
					with open(aln[0]+'.fasta' , 'w') as handle:
						handle.write(outstr)
				tmfastas.append(aln[0]+'.fasta')
		print( len(fastas))
		print( 'alignments by fatcat')
		save_obj( tmfastas, filedir+'tmfastas')
		save_obj(  distances, filedir +'distances')
		print( 'done writing fastas')

	if make_distmat == True:
			pdblist = glob.glob(filedir +'*.pdb')
			if SDM == True:
				distances = load_obj(filedir + 'distances')
				maxdist = 0
				pdblist = []
				for out in distances.values():
					pdblist += out[0]
					for dist in out[2]:
						if dist > maxdist:
							maxdist = dist
				pdblist = list(set(pdblist))
				print( 'creating SDM distmat for Fatcat output')
				distmat = np.zeros((len(pdblist),len(pdblist)))
				pvalmat = np.zeros((len(pdblist),len(pdblist)))
				#calculate SDM based distmat
				for out in distances.values():
					lensmallprot = min(out[1])
					SRMSs = []
					count = 0
					pdbs = out[0]
					distances = out[2]
					counts = out[3]
					pval = out[4]
					for block in distances:
						count += counts[block]
						SRMSs.append(counts[block]*(1-distances[block]/maxdist))
					PFTE = count / lensmallprot
					avgSRMS = sum(SRMSs)/count
					w1 = (1 -PFTE + 1 -avgSRMS) /2
					w2 = (PFTE +avgSRMS)/2
					sdm = -100* math.log(w1*PFTE + w2 *avgSRMS)
					distmat[ pdblist.index(pdbs[0]),pdblist.index(pdbs[1])] = sdm
					#low pval -> high score
					pvalmat[ pdblist.index(pdbs[0]),pdblist.index(pdbs[1])] = pval
				distmat += distmat.T
				#distmat /= np.amax(distmat)
				pvalmat = -1 / np.log(pvalmat)
				pvalmat += pvalmat.T
				np.fill_diagonal(pvalmat, 0)


				save_obj(   distmat , filedir + 'SDMmat')
				save_obj(   pvalmat , filedir + 'pvalmat')
				save_obj(   pdblist , filedir + 'structs')

				sdmd=distmat_to_txt( pdblist , distmat, filedir , 'SDMmat')
				pvald=distmat_to_txt( pdblist , pvalmat, filedir , 'pvalmat')

				runFastme(fastmepath , sdmd)
				runFastme(fastmepath , pvald)


	if make_TMdistmat == True:
		print( 'creating TMalign based distmat')

		TMRealignresults = load_obj( filedir +'TMRealignfiles' )
		print(TMRealignresults[0])

		tmlist = []
		namelist = []
		TMdistmat = np.zeros((len(pdblist),len(pdblist)))
		
		print(pdblist)
		for args in TMRealignresults:	
			filename, pdbs2 = args
			try:
				distance,pdbs,outstr =TMaligntoDistance(filename)
			
				TMdistmat[ pdblist.index(pdbs2[0]),pdblist.index(pdbs2[1])] = 1-distance
			except:
				distance,pdbs,outstr =TMaligntoDistance(filename, verbose = True)
				TMdistmat[ pdblist.index(pdbs2[0]),pdblist.index(pdbs2[1])] = 1
				print( pdbs2 )		
		
		TMdistmat += TMdistmat.T
		print( TMdistmat )
		save_obj( TMdistmat , filedir + 'TMdistmat' )

		num2tax = { n:pdblist[n] for n in range(len(pdblist)) }
		save_obj( list(num2tax.keys()) , filedir + 'keys' )
		save_obj( pdblist , filedir + 'structs' )
		tmd  = distmat_to_txt( [ str(i) for i in list(num2tax.keys()) ], TMdistmat, filedir , 'TMRealigndistmat')
		runFastme(fastmepath , tmd)
		print(namelist)



	if merge == True:
	## TODO: clean this up. names not correctly transferring

		#alignments = glob.glob(filedir + '*TMRealign.fasta')
		alignments = load_obj(filedir +'tmfastas')
		print( alignments)
		print( 'merging all alignments')
		nextmerge = []
		i =0
		fastasleft = list(alignments)
		print(len(fastasleft))
		fastasleft = [ check_identical(fasta) for fasta in fastasleft ]
		fastasleft = [ fasta for fasta in fastasleft if fasta ]
		print(fastasleft)
		ids = [ id for fasta in fastasleft for id in  check_ids(fasta) ]
		print(set(ids))
		count = 1
		for i, fasta in enumerate(alignments):
			if i >0:
				if len( set(check_ids(fasta)).difference(set(check_ids(filedir+str(count-1)+'total.fasta'))) ) > 0 :
					mergeAlign(clustaloPath,fasta, filedir+str(count-1)+'total.fasta' , filedir+str(count)+'total.fasta' )
					seen = set([ id for  id in  check_ids(filedir+str(count)+'total.fasta') ])
					count+=1

			else:
				with open( filedir+str(count-1)+'total.fasta' , 'w' ) as out:
					with open(fasta ,'r') as infile:
						out.write(infile.read())

				i+=1
		print('DONE')
	

		tmMat=distmat_to_txt( pdblist , TMdistmat, filedir , 'TMmat')
		runFastme(fastmepath , tmMat)

		save_obj(   TMdistmat , filedir + 'TMRealigndistmat')
		save_obj(   namelist , filedir + 'structs')

if __name__ == '__main__':
	#where a your structures? ... should all be single chain with decent homology.
	#filedir = '../Structures/structalign_finalset/'
	filedir = sys.argv[1]
	structalign(filedir)
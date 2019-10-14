#!/usr/bin/python
import time
import os
import sys
import csv
import numpy
from optparse import OptionParser
import schrodinger.application.canvas.base as canvas
import schrodinger.application.canvas.utils as utils
from Graph import *
from progressbar import *
		
class MatchEngine(object):
	def __init__(self, substructs, types, extended=False, count=False, triplets=False, all=False,linkers=False):
		numpy.set_printoptions(threshold=numpy.nan, linewidth=numpy.nan)
		self.linkers=linkers
		self.types=types
		self.extended=extended
		self.count=count
		self.triplets=triplets
		self.all=all
		self.substructs=substructs
		self.matcher = canvas.ChmQueryMatcher()
		self.unused=[]
		self.onceHit={}
		self.struct=None
		self.subPaths={}
		self.paths=[]
		self.matrix=None
		self.matchTime = 0
		self.graphTime=0
		self.matrixTime = 0

	def readStruct(self, struct):
		self.regular={}
		self.matrix=None
		self.coords=None
		self.line=None
		self.structGraph=None
		self.unused=[]
		self.subPaths={}
		self.paths=[]
		self.onceHit={}
		self.struct=None
		self.struct = struct
		self.structGraph = Graph()
		self.unused=set(self.struct.getAtoms())

		for atom in self.unused:
			self.structGraph.addNode(atom.getLabel(), map(canvas.ChmAtom.getLabel,atom.getNeighbors(False,False)))
		
		self.structGraph.searchRings()

		self.evaluateSmarts()
		for node in self.structGraph.nodes():
			if type(node)==float:
				self.substructDFS(node, path=[node])
				self.subPaths[node]=self.paths
				self.paths=[]
		a = time.clock()
		if self.all==True:
			self.buildExtendedMatrix()
		else:
		    self.buildMatrix()
		b = time.clock()
		self.matrixTime += b-a
		
	def evaluateSmarts(self):
		a = time.clock()
		
		for i in range(len(self.substructs)):
			sub = self.substructs[i]
			appendRings = False
			temp = sub
			name = i
			if '!R0' in temp:
				appendRings=True
			if '!R' in temp:
				temp = ''.join(temp.split('!R'))
			if 'R0' in temp:
				temp = ''.join(temp.split('R0'))
			if 'a' in temp:
				appendRings = True
			if 'R' in temp:
				appendRings = True
			matchedSubs = self.matcher.match(canvas.ChmQuery(sub),self.struct)
			if len(matchedSubs) != 0:
				if self.count == True:
					self.regular[i]=len(matchedSubs)
				else:
				    self.regular[i]=1
				count = 0
				for hit in matchedSubs:
					count +=1
					endname=name+(float(count)/100)
					atoms = set(map(canvas.ChmAtom.getLabel,set(hit.getMatchedAtoms())))
					for n in atoms:
						if 'H' in n:
							atoms = atoms - set([n])
					if appendRings == True:
						for ring in self.structGraph.circles:
							if len(atoms & ring) > 0:
								atoms = atoms | set(ring)
					repeated = False
					for key in self.onceHit.keys():
						if str(key).split('.')[0] == str(name) and self.onceHit[key] == atoms:
							repeated = True
					if repeated == False:
						self.onceHit[endname]=atoms
					else:
						continue
			else:
				self.regular[i]=0
				self.onceHit[name] = set([])
		b = time.clock()
		self.matchTime += b-a
		subs = sorted(self.onceHit.items(), key = lambda x: len(x[1]),reverse=True)
		for node in self.structGraph.nodes():
			for sub in self.onceHit.keys():
				if len(self.structGraph.graph[node] & self.onceHit[sub]) > 0:
					self.structGraph.graph[node] = self.structGraph.graph[node]| set([sub])
		usedAtoms=set([])
		for s in subs:
			if len(s[1]) >0:
				connections = set([])
				for atom in s[1]:
					usedAtoms = usedAtoms | set([atom])
					connections = connections | self.structGraph.graph[atom]
				connections = connections - s[1]
				for conn in connections:
					self.structGraph.addEdge(conn,s[0])
				self.structGraph.addNode(s[0], list(connections))
		for node in usedAtoms:
			self.structGraph.deleteNode(node)
		self.structGraph.removeSelfEdges()
		c = time.clock()
		self.graphTime += c-b
			
	def substructDFS(self, root, path=[], prevnode=None):
		children = self.structGraph.expandNode(root, prevnode)
		if children != []:
			for child in children:
				if child not in path:
					if not type(child)==float:
						newpath = path[:]
						newpath.append(child)
						self.substructDFS(child,newpath,prevnode=root)
					else:
						newpath = path[:]
						newpath.append(child)
						self.paths.append(newpath)
		else:
			pass
			
	
	def buildExtendedMatrix(self):
		subs = set(self.subPaths.keys())		
		indexes = list(sorted(self.onceHit.keys()))
		size= len(self.onceHit)
		matrix = numpy.zeros((size, size), dtype=int)
		coords = {}
		for sub in subs:
			if self.subPaths[sub] is not None:
				for path in self.subPaths[sub]:
					end = path[-1]
					pair = str(sub)+','+str(end)
					if self.extended:
						if len(path) == 2:
							if self.onceHit[sub] & self.onceHit[end] != set([]):
								if self.onceHit[sub] > self.onceHit[end] or self.onceHit[sub] < self.onceHit[end]:
									matrix[indexes.index(sub)][indexes.index(end)]= 1
								else:
									matrix[indexes.index(sub)][indexes.index(end)] = 2
							else:
								matrix[indexes.index(sub)][indexes.index(end)] = 4
						elif len(path) > 2:
							if matrix[indexes.index(sub)][indexes.index(end)] < 3:
								matrix[indexes.index(sub)][indexes.index(end)] = 3
					elif self.linkers:
						matrix[indexes.index(sub)][indexes.index(end)]= len(path)
						if coords.has_key(pair):
							coords[pair].append(len(path)) 
						else:
							coords[pair] = [len(path)]

					else:
						matrix[indexes.index(sub)][indexes.index(end)]= 1
						if coords.has_key(pair):
							if self.count:
								coords[pair] += 1
							else:
								continue
						else:
							coords[pair] = 1
		self.coords=coords
		self.matrix=(matrix,indexes)
	
	def buildMatrix(self):
		subs = list(sorted(self.subPaths.keys()))
		size= len(self.substructs)
		matrix = numpy.zeros((size, size), dtype=int)
		line = ''
		coords={}
		triplets = {}
		for sub in subs:
			if self.subPaths[sub] != None:
				for path in self.subPaths[sub]:
					end = path[-1]
					endindex = int(str(path[-1]).split('.')[0])
					startindex = int(str(path[0]).split('.')[0])
					pair = ','.join(map(str,sorted([startindex+1,endindex+1])))
					if self.triplets:
						for path in self.subPaths[end]:
							end2= path[-1]
							if end2 != sub:
								if int(startindex) < int(end2):
									start = int(startindex+1)
									middle = int(endindex+1)
									end=int(end2+1)
								else:
									start = int(end2+1)
									middle = int(endindex+1)
									end= int(startindex+1)
								if triplets.has_key(start):
									if triplets[start].has_key(middle):
										triplets[start][middle][end]=1
									else:
										triplets[start][middle]={end:1}
								else:
									triplets[start]={middle:{end:1}}
								
					if self.extended == True:
						if len(path) == 2:
							if self.onceHit[sub] & self.onceHit[end] != set([]):
								if self.onceHit[sub] > self.onceHit[end] or self.onceHit[sub] < self.onceHit[end]:
									if matrix[startindex][endindex] < 1:
										matrix[startindex][endindex]= 1
									if coords.has_key(pair):
										continue
									else:
										coords[pair] = 1
								else :
									if matrix[startindex][endindex] < 2:
										matrix[startindex][endindex] = 2
									if coords.has_key(pair):
										if coords[pair] < 2: 
											coords[pair] = 2
									else:
										coords[pair] = 2
							else:
								matrix[startindex][endindex] =4
								coords[pair] = 4

						elif len(path) > 2:
							if matrix[startindex][endindex] < 3:
								matrix[startindex][endindex] = 3
							if coords.has_key(pair):
								if coords[pair] < 3: 
									coords[pair] = 3
							else:
								coords[pair] = 3
					else:
						matrix[startindex][endindex]= 1
						if coords.has_key(pair):
							if self.count == True:
								coords[pair] += 1
							else:
								continue
						else:
							coords[pair] = 1
		for start in sorted(triplets.keys()):
			for middle in sorted(triplets[start]):
				for end in sorted(triplets[start][middle]):
					triplet = ','.join(map(str,[start,middle,end]))
					if coords.has_key(triplet):
						if self.count == True:
							coords[triplet] += 1
						else:
							pass
					else:
						coords[triplet]=1			
		self.matrix=(matrix,subs)
		self.coords=coords
		lineno=0
		for linia in matrix.tolist():
			fpline = linia[lineno:]
			line+=''.join(map(str,fpline))
			lineno+=1
		self.line=line
		
class Writer:
	def __init__(self,fname):
		self.log = open(fname,'w')
		
	def writeIntro(self,smarts,IDs):
		self.log.write('<SMARTS>\n')
		self.log.write(','.join(smarts)+'\n')
		self.log.write('<Substructure IDs>\n')
		self.IDs = IDs
		self.log.write(','.join(IDs)+'\n')
		self.log.write('<Substructure count>\n')
		self.log.write(str(len(smarts))+'\n\n')
	
	def write(self,line):
		self.log.write(line+'\n')
		self.log.flush()
	
	def writeLine(self,name,line):
		self.log.write('<Name>\n')
		self.log.write(name+'\n')
		self.log.write('<start>\n')
		self.log.write(line+'\n')
		self.log.write('<end>\n')
		self.log.write('-'*50+'\n')

	def writeCoord(self,name,coord):
		self.log.write(name+' ')
		line = []
		n = []
		tri=[]
		for a in coord.keys():
			try:
				b = map(int,a.split(','))
			except:
				b = map(float,a.split(','))
			if len(b) == 2:
				b.append(coord[a])
				n.append(b)
			else:
				tri.append(a+':'+str(coord[a]))
		self.log.write(' '.join([','.join(map(str,z)[:-1])+':'+str(z[-1]) for z in sorted(n)])+' ')
		self.log.write(' '.join(tri)+'\n')

	def writeMatrix(self,name,matrix,ext):
		self.log.write('<Name>\n')
		self.log.write(name+'\n')
		self.log.write('<Substructure hits>\n')
		subs = []		
		if ext == True:
			for n in matrix[1]:
				subs.append(self.IDs[int(round(n))])
		else:
			subs = self.IDs
		self.log.write(','.join(subs)+'\n')
		self.log.write('<start>\n')
		for line in matrix[0].tolist():
			z = map(str,line)
			self.log.write(' '.join(z)+'\n')
		self.log.write('<end>\n')
		self.log.write('-'*50+'\n')
			
class Controller:
	def __init__(self, perc, freq, infile, outfile, keys, ext, sil, mat, lin, crd, reg, cnt, tri, all):
		self.sil=sil
		self.towrite=''
		if mat:
			self.towrite+='m'
		if lin:
			self.towrite+='l'
		if crd:
			self.towrite+='c'
		if reg:
			self.towrite+='r'
		if self.towrite == '':
			self.towrite ='m'
		self.a = time.clock()
		self.log=Writer(outfile.split('.')[0].rstrip('_')+'.log')
		self.log.write('2D substructure fingerprint log\nCommand: '+' '.join(sys.argv)+'\nStart time: '+time.ctime())
		self.perc = perc
		self.freq=freq
		self.infile=infile
		self.outfile=outfile
		self.keyfile=keys
		self.structures=None
		self.all = all
		lic = self.checkLicense()
		self.ext=ext
		self.count=cnt
		self.triplet=tri
		if not lic.isValid():
			sys.stderr.write("CANVAS_FULL license is not valid\n")
			self.log.write("CANVAS_FULL license is not valid\n")
			sys.exit(1)
		if self.infile.endswith('.sdf'):
			self.readSDF(self.infile)
		elif self.infile.endswith('.smi') or self.infile.endswith('.smiles'):
			self.readSMI(self.infile)
		else:
			print("Provide proper input file")
			self.log.write('Improper input file\n')
			sys.exit(1)
			
		self.total=len(self.structures)
		self.readCSV(self.keyfile)
		if self.perc != 0 and self.freq != None:
			self.usedKeys=self.parseSubFP(self.freq,self.perc)
		else:
			self.usedKeys='All'
		self.substructs, self.subfields = self.findSubFPInCSV(self.usedKeys, self.keys)
		self.subs = len(self.substructs)
		outs = []
		if 'm' in self.towrite:
			outs.append('Fingerprint Matrix')
		if 'l' in self.towrite:
			outs.append('Linear Fingerprint')
		if 'c' in self.towrite:
			outs.append('Coordinate Fingerprint')
		if 'r' in self.towrite:
			outs.append('Original Fingerprint')	
		self.log.write('Cutoff: '+str(self.perc)+'%\nNumber of substructures: '+str(self.subs)+'\nNumber of unique bits: '+str(((self.subs*(self.subs+1))/2))+'\nNumber of structures: '+str(self.total)+'\n'+'Writing into: '+', '.join(outs)+'\n')
		self.matcher = MatchEngine(self.substructs, self.towrite, self.ext, self.count, self.triplet, self.all)
#############TODO
		if 'm' in self.towrite:
			self.output=Writer(self.outfile)
			self.output.writeIntro(self.substructs,self.subfields)
		if 'l' in self.towrite:
			self.output2=Writer(self.outfile.replace('.dat','_lin.dat'))
			self.output2.writeIntro(self.substructs,self.subfields)
		if 'c' in self.towrite:
			self.output3=Writer(self.outfile.replace('.dat','_coord.dat'))
		if 'r' in self.towrite:
			self.output4=Writer(self.outfile.replace('.dat','_orig.dat'))

		if self.sil == False:
			print('Calculating 2D fingerprint for %s, using %s keys and %s cutoff' %(self.infile.split('/')[-1], self.keyfile.split('/')[-1].split('.')[0],str(self.perc)+'%'))
			self.bar=ProgressBar(widgets=[Percentage(), Bar()], maxval=self.total)
			self.bar.start()
		count = 0
		for struct in self.structures:
			name = struct.getProperty('CMPD_CHEMBLID')
			if not name:
				name = struct.getName()
			self.matcher.readStruct(struct)
			if 'm' in self.towrite:
				self.output.writeMatrix(name,self.matcher.matrix, self.ext)
			if 'l' in self.towrite:
				self.output2.writeLine(name,self.matcher.line)
			if 'c' in self.towrite:
				self.output3.writeCoord(name,self.matcher.coords)
			if 'r' in self.towrite:
				line = ','.join(str(x[1]) for x in sorted(self.matcher.regular.items()))
				self.output4.write(name+':'+line)
			count+=1
			if self.sil==False:
				self.bar.update(count)
		if self.sil==False:
			self.bar.finish()
		self.endTime=time.clock()
		total_time=self.endTime-self.a
		self.log.write('Total time elapsed: '+str(total_time)+' sec.\nSpeed: '+str(round(float(total_time)/self.total,2))+' sec/mol')
		self.log.write('Time elapsed for substructure reading: %i sec. ( %i%% )\nTime elapsed for graph construction: %i sec. ( %i%% )\nTime elapsed for matrix construction: %i sec. ( %i%% )' %\
		(self.matcher.matchTime, int(self.matcher.matchTime*100/float(total_time)),self.matcher.graphTime, int(self.matcher.graphTime*100/float(total_time)),self.matcher.matrixTime, int(self.matcher.matrixTime*100/float(total_time))))
		
			
					
	def readSDF(self,fname):
		self.structures = [a for a in canvas.ChmSDReader(fname)]

	def readSMI(self,fname):
		self.structures = [a for a in canvas.ChmSmilesFileReader(fname)]
	
	def readCSV(self,fname):
		data = []
		with open(fname, 'rb') as finput:
			content = csv.reader(finput, delimiter=',')
			for row in content:
				data.append(row)
		data = data[1:]
		self.keys=data

	def parseSubFP(self, fname, percent):
		subno = []
		for line in open(fname, 'r'):
			line = line.strip().split()
			if float(line[1].strip('%')) >= float(percent):
				subno.append(line[0])
		return subno

	def findSubFPInCSV(self,subfp, csv):
		found = []
		substr = []
		for j in csv:
			if subfp == 'All':
				found.append(j[-1])
				substr.append(j[0])				
			else:
				for i in subfp:
					if j[0] == i:
						found.append(j[-1])
						substr.append(j[0])
		return found, substr
		
	def checkLicense(self):
		try:
			lic = utils.get_license()
		except:
			sys.stderr.write("Unable to checkout CANVAS_FULL license\n")
			sys.exit(1)
		return lic
		

def main():
	parser = OptionParser(usage="usage: %prog -p [PERCENT] -f [FREQ] -k [KEYS] (-o [OUTPUT] -e [EXT]) [args]")
	parser.add_option("-p", "--perc", dest="perc", default=0,help="INT or FLOAT substructure occurence cutoff")
	parser.add_option("-f","--freq", dest="freq", default=None,help=".freq file containing substructure occurence")
	parser.add_option("-o", "--output", dest="outfile", help="path to the output file")
	parser.add_option("-k", "--keys", dest="keys", help="CSV file with substructure SMARTS keys")
	parser.add_option("-e", "--extended", dest="ext", action="store_true", default=False, help="Construct extended substructure matrix")
	parser.add_option("-s", "--silent", dest="sil", action="store_true", default=False, help="Force a silent run")
	parser.add_option("-m", "--matrix", dest="mat", action="store_true", default=False, help="Write full connection matrix")	
	parser.add_option("-l", "--linear", dest="lin", action="store_true", default=False, help="Write linearized fingerprit")
	parser.add_option("-c", "--coord", dest="crd", action="store_true", default=False, help="Write coordinates fp")
	parser.add_option("-r", "--regular", dest="reg", action="store_true", default=False, help="Write regular fingerprints")
	parser.add_option("--count", dest="cnt", action="store_true", default=False, help="Write counted substructures")
	parser.add_option("-t","--triplets", dest="tri", action="store_true", default=False, help="Write triplet connections")
	parser.add_option("-a","--all", dest="all", action="store_true", default=False, help="Write all substructures")

	(options, args) = parser.parse_args()
	if not args:
		print("Provide proper input file")
		sys.exit(1)
	elif len(args)>1:
		print("Too many arguments")
		sys.exit(1)
	if not options.perc:
		if options.sil == False:
			print("No cutoff defined, calculating unbiased fingerprint")
	if not options.keys:
		print("No CSV with structure keys provided")
		sys.exit(1)
	if not options.freq:
		if options.sil == False:
			print("No substructure occurence file defined, calculating unbiased fingerprint")
	if not options.outfile:
		name = ''
		name += args[0].split('/')[-1].split('.')[0]
		#perc = str(options.perc).replace('.','')
		name += '_'+options.keys.split('/')[-1].split('.')[0]
		if options.tri:
			name += '_tri'
		if options.cnt:
			name += '_count'
		options.outfile = name+'.dat'
	controller = Controller(options.perc, options.freq, args[0], options.outfile, options.keys, options.ext,options.sil,options.mat,options.lin,options.crd,options.reg,options.cnt,options.tri, options.all)
	c = time.clock()
		
	

if __name__ == "__main__":
	main()


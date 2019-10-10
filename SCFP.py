#!/usr/bin/python
import time
import os
import sys
import csv
import argparse
import pandas as pd
import numpy
from rdkit import Chem
from rdkit.Chem import PandasTools
from tqdm import tqdm
from functools import partial
from Graph import *

def GetProperLabel(struct,index):
	return Chem.Atom.GetSymbol(Chem.Mol.GetAtomWithIdx(struct,index))+str(index)


class MatchEngine(object):
    def __init__(self, substructs, types, extended=False, count=False, triplets=False, all=False,linkers=False):
        self.linkers=linkers
        self.types=types
        self.extended=extended
        self.count=count
        self.triplets=triplets
        self.all=all
        self.substructs=substructs
        self.unused=[]
        self.onceHit={}
        self.struct=None
        self.subPaths={}
        self.paths=[]
        self.matrix=None
        self.matchTime = 0
        self.graphTime=0
        self.matrixTime = 0

    def readStruct(self, struct, calc_matrix):
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
        self.unused=set(self.struct.GetAtoms())
        self.calc_matrix = calc_matrix

        for atom in self.unused:
            self.structGraph.addNode(GetProperLabel(self.struct,Chem.Atom.GetIdx(atom)), [GetProperLabel(self.struct,Chem.Atom.GetIdx(x)) for x in Chem.Atom.GetNeighbors(atom)])

        self.structGraph.searchRings()

        self.evaluateSmarts()
        if calc_matrix:
            for node in self.structGraph.nodes():
                if type(node)==float:
                    self.substructDFS(node, path=[node])
                    self.subPaths[node]=self.paths
                    self.paths=[]
        a = time.process_time()
        if self.all==True:
            self.buildExtendedMatrix()
        else:
            self.buildMatrix(calc_matrix)
        b = time.process_time()
        self.matrixTime += b-a
		
    def evaluateSmarts(self):
        a = time.process_time()

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
            patt = Chem.MolFromSmarts(temp)
            matchedSubs = self.struct.GetSubstructMatches(patt)
            if len(matchedSubs) != 0:
                if self.count == True:
                    self.regular[i]=len(matchedSubs)
                else:
                    self.regular[i]=1
                if self.calc_matrix:
                    count = 0
                    for hit in matchedSubs:
                        count +=1
                        endname=name+(float(count)/100)
                        atoms = set([Chem.Atom.GetSymbol(Chem.Mol.GetAtomWithIdx(self.struct,a))+str(a) for a in hit])
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
        b = time.process_time()
        self.matchTime += b-a
        if self.calc_matrix:
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
        c = time.process_time()
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
            if self.subPaths[sub] != None:
                for path in self.subPaths[sub]:
                    end = path[-1]
                    pair = str(sub)+','+str(end)
                    if self.extended==True:
                        if len(path) == 2:
                            if self.onceHit[sub] & self.onceHit[end] != set([]):
                                if self.onceHit[sub] > self.onceHit[end] or self.onceHit[sub] < self.onceHit[end]:
                                    matrix[indexes.index(sub)][indexes.index(end)]= 1
                                else :
                                    matrix[indexes.index(sub)][indexes.index(end)] = 2
                            else:
                                matrix[indexes.index(sub)][indexes.index(end)] = 4
                        elif len(path) > 2:
                            if matrix[indexes.index(sub)][indexes.index(end)] < 3:
                                matrix[indexes.index(sub)][indexes.index(end)] = 3
                    elif self.linkers == True:
                        matrix[indexes.index(sub)][indexes.index(end)]= len(path)
                        if pair in coords.keys():
                            coords[pair].append(len(path)) 
                        else:
                            coords[pair] = [len(path)]

                    else:
                        matrix[indexes.index(sub)][indexes.index(end)]= 1
                        if pair in coords.keys():
                            if self.count == True:
                                coords[pair] += 1
                            else:
                                continue
                        else:
                            coords[pair] = 1
        self.coords=coords
        self.matrix=(matrix,indexes)
	
    def buildMatrix(self,calc_matrix):
        subs = list(sorted(self.subPaths.keys()))
        size= len(self.substructs)
        matrix = numpy.zeros((size, size), dtype=int)
        
        coords={}
        triplets = {}
        if calc_matrix:
            for sub in subs:
                if self.subPaths[sub] is None:
                    continue
                for path in self.subPaths[sub]:
                    end = path[-1]
                    endindex = int(str(path[-1]).split('.')[0])
                    startindex = int(str(path[0]).split('.')[0])
                    pair = ','.join(list(map(str,sorted([startindex+1,endindex+1]))))
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
                                if start in triplets.keys():
                                    if middle in triplets[start].keys():
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
                                    if pair in coords.keys():
                                        continue
                                    else:
                                        coords[pair] = 1
                                else :
                                    if matrix[startindex][endindex] < 2:
                                        matrix[startindex][endindex] = 2
                                    if pair in coords.keys():
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
                            if pair in coords.keys():
                                if coords[pair] < 3: 
                                    coords[pair] = 3
                            else:
                                coords[pair] = 3
                    else:
                        matrix[startindex][endindex]= 1
                        if pair in coords.keys():
                            if self.count == True:
                                coords[pair] += 1
                            else:
                                continue
                        else:
                            coords[pair] = 1
            for start in sorted(triplets.keys()):
                for middle in sorted(triplets[start]):
                    for end in sorted(triplets[start][middle]):
                        triplet = ','.join(list(map(str,[start,middle,end])))
                        if triplet in coords.keys():
                            if self.count == True:
                                coords[triplet] += 1
                            else:
                                pass
                        else:
                            coords[triplet]=1			
            self.matrix=(matrix,subs)
            self.coords=coords
        line = numpy.array([])
        for lineno in range(matrix.shape[0]):
            line = numpy.append(line,matrix[lineno][lineno:])
        self.line=','.join(map(str,line.tolist()))
		
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
		self.log.write(line)
	
	def writeTestLine(self,name,line):
		self.log.write('<Name>\n')
		self.log.write(name+'\n')
		self.log.write('<start>\n')
		self.log.write(line+'\n')
		self.log.write('<end>\n')
		self.log.write('-'*50+'\n')

	def writeLine(self,name,line,csv):
		if csv:
			line = ','.join([a for a in line])
			self.log.write(name.strip()+','+line+'\n')
		else:
			self.log.write(name.strip()+' '+line+'\n')

	def writeCoord(self,name,coord):
		self.log.write(name+' ')
		line = []
		n = []
		tri=[]
		for a in coord.keys():
			try:
				b = list(map(int,a.split(',')))
			except:
				b = list(map(float,a.split(',')))
			if len(b) == 2:
				b.append(coord[a])
				n.append(b)
			else:
				tri.append(a+':'+str(coord[a]))
		self.log.write(' '.join([','.join(list(map(str,z))[:-1])+':'+str(z[-1]) for z in sorted(n)])+' ')
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
			z = list(map(str,line))
			self.log.write(' '.join(z)+'\n')
		self.log.write('<end>\n')
		self.log.write('-'*50+'\n')
			
class Controller:
    def __init__(self, perc, freq, infile, outfile, keys, ext, sil, mat, lin, crd, reg, cnt, tri, all, csv):
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
        self.a = time.process_time()
        self.log=Writer(outfile.split('.')[0].rstrip('_')+'.log')
        self.log.write('2D substructure fingerprint log\nCommand: '+' '.join(sys.argv)+'\nStart time: '+time.ctime())
        self.perc = perc
        self.freq=freq
        self.infile=infile
        self.outfile=outfile
        self.keyfile=keys
        self.structures=None
        self.all = all
        self.ext=ext
        self.count=cnt
        self.triplet=tri
        self.data_frame = None

        if self.infile.endswith('.sdf'):
            self.readSDF(self.infile)
        elif self.infile.endswith('.smi') or self.infile.endswith('.smiles'):
            self.readStructsFromCSV(self.infile)
            
        elif self.infile.endswith('.csv'):
            self.readStructsFromCSV(self.infile)
        else:
            print("Provide proper input file")
            self.log.write('Improper input file\n')
            sys.exit(1)
        self.total=self.data_frame.shape[0]
        self.readCSV(self.keyfile)
        if self.perc != 0 and self.freq != None:
            self.usedKeys=self.parseSubFP(self.freq,self.perc)
        else:
            self.usedKeys='All'
        self.substructs, self.subfields = self.findSubFPInCSV(self.usedKeys, self.keys)
        self.subs = len(self.substructs)
        outs = []
        self.calc_matrix = False
        if 'm' in self.towrite:
            outs.append('Fingerprint Matrix')
            self.calc_matrix = True
        if 'l' in self.towrite:
            outs.append('Linear Fingerprint')
            self.calc_matrix = True
        if 'c' in self.towrite:
            outs.append('Coordinate Fingerprint')
            self.calc_matrix = True
        if 'r' in self.towrite:
            outs.append('Original Fingerprint')	
            
        self.log.write('Cutoff: '+str(self.perc)+'%\nNumber of substructures: '+str(self.subs)+'\nNumber of unique bits: '+str(((self.subs*(self.subs+1))/2))+'\nNumber of structures: '+str(self.total)+'\n'+'Writing into: '+', '.join(outs)+'\n')

        if 'm' in self.towrite:
            self.output=Writer(self.outfile)
            self.output.writeIntro(self.substructs,self.subfields)
        if 'l' in self.towrite:
            if csv:
                self.output2=Writer(self.outfile.replace('.dat','_lin.csv'))
                
            else:
                self.output2=Writer(self.outfile.replace('.dat','_lin.dat'))
            self.output2.write("SMILES"+","+",".join(["2D"+keys.split("/")[-1].split(".")[0]+"_{number}".format(number = num+1) for num in range(int((self.subs*self.subs-self.subs)/2+self.subs))])+"\n")
        if 'c' in self.towrite:
            self.output3=Writer(self.outfile.replace('.dat','_coord.dat'))
        if 'r' in self.towrite:
            self.output4=Writer(self.outfile.replace('.dat','_orig.dat').replace('2D',''))

        if self.sil == False:
            print('Calculating 2D fingerprint for %s, using %s keys and %s cutoff' %(self.infile.split('/')[-1], self.keyfile.split('/')[-1].split('.')[0],str(self.perc)+'%'))
        count = 0
        
        self.matcher = MatchEngine(self.substructs, self.towrite, self.ext, self.count, self.triplet, self.all)        
        for struct in tqdm(self.data_frame["Structure"]):
            if struct.HasProp("Name"):
                name = struct.GetProp('Name')
            else:
                name = Chem.MolToSmiles(struct)
            self.matcher.readStruct(struct, self.calc_matrix)
            if 'm' in self.towrite:
                self.output.writeMatrix(name,self.matcher.matrix, self.ext)
            if 'l' in self.towrite:
                self.output2.writeLine(name,self.matcher.line,csv)
            if 'c' in self.towrite:
                self.output3.writeCoord(name,self.matcher.coords)
            if 'r' in self.towrite:
                line = ','.join(str(x[1]) for x in sorted(self.matcher.regular.items()))
                self.output4.write(name+':'+line)
            count+=1
        self.endTime=time.process_time()
        total_time=self.endTime-self.a
        self.log.write('Total time elapsed: '+str(total_time)+' sec.\nSpeed: '+str(round(float(total_time)/self.total,2))+' sec/mol')
	
    def _worker_func(self,df):
        func = partial(self._calculate_fp, calc_matrix=self.calc_matrix)
        df["Structure"].apply(func)

    def _calculate_fp(self, struct, calc_matrix):
        matcher = MatchEngine(self.substructs, self.towrite, self.ext, self.count, self.triplet, self.all)
        matcher.readStruct(struct, calc_matrix)
        name = struct.GetProp('Name')
        if 'm' in self.towrite:
            self.output.writeMatrix(name,matcher.matrix, self.ext)
        if 'l' in self.towrite:
            self.output2.writeLine(name,matcher.line,csv)
        if 'c' in self.towrite:
            self.output3.writeCoord(name,matcher.coords)
        if 'r' in self.towrite:
            line = ','.join(str(x[1]) for x in sorted(matcher.regular.items()))
            self.output4.write(name+':'+line)
        self.t.update(n=1)
        del matcher
					
    def readSDF(self,fname):
        self.data_frame = PandasTools.LoadSDF(fname)
        self.data_frame.rename(columns = {"ROMol":"Structure"}, inplace = True)

    def readStructsFromCSV(self,fname):
        column = None
        self.data_frame = pd.read_csv(fname)
        for col in self.data_frame.columns:
            if col in ["SMILES","smiles","Smiles"]:
                column = col
        if column is None:
            raise Exception("No SMILES column found!")
        self.data_frame["Structure"] = self.data_frame[column].apply(self._get_structure)

    def _get_structure(self,smiles):
        struct = Chem.MolFromSmiles(smiles)
        struct.SetProp('Name',smiles)
        return struct


    def readCSV(self,fname):
        data = []
        with open(fname, 'r') as finput:
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

def main():
	parser = argparse.ArgumentParser(usage="usage: -p [PERCENT] -f [FREQ] -k [KEYS] (-o [OUTPUT] -e [EXT]) [args]")
	parser.add_argument("-p", "--perc", dest="perc", default=0,help="INT or FLOAT substructure occurence cutoff")
	parser.add_argument("-f","--freq", dest="freq", default=None,help=".freq file containing substructure occurence")
	parser.add_argument("-o", "--output", dest="outfile", default = None, help="path to the output file")
	parser.add_argument("-i", "--input", dest="infile", help="path to the input file")
	parser.add_argument("--csv", dest="csv", action="store_true", help="Write linear output as CSV")
	parser.add_argument("-k", "--keys", dest="keys", help="CSV file with substructure SMARTS keys")
	parser.add_argument("-e", "--extended", dest="ext", action="store_true", default=False, help="Construct extended substructure matrix")
	parser.add_argument("-s", "--silent", dest="sil", action="store_true", default=False, help="Force a silent run")
	parser.add_argument("-m", "--matrix", dest="mat", action="store_true", default=False, help="Write full connection matrix")	
	parser.add_argument("-l", "--linear", dest="lin", action="store_true", default=False, help="Write linearized fingerprit")
	parser.add_argument("-c", "--coord", dest="crd", action="store_true", default=False, help="Write coordinates fp")
	parser.add_argument("-r", "--regular", dest="reg", action="store_true", default=False, help="Write regular fingerprints")
	parser.add_argument("--count", dest="cnt", action="store_true", default=False, help="Write counted substructures")
	parser.add_argument("-t","--triplets", dest="tri", action="store_true", default=False, help="Write triplet connections")
	parser.add_argument("-a","--all", dest="all", action="store_true", default=False, help="Write all substructures")

	options = parser.parse_args()
	if not options.infile:
		print("Provide proper input file")
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
	name = ''
	if not options.outfile:
		name = ''
		name += options.infile.split('/')[-1].split('.')[0]
		name += '_2D'+options.keys.split('/')[-1].split('.')[0]
		if options.tri:
			name += '_tri'
		if options.cnt:
			name += '_count'
		name+='.dat'
	else:
		name = options.outfile
	controller = Controller(options.perc, options.freq, options.infile, name, options.keys, options.ext,options.sil,options.mat,options.lin,options.crd,options.reg,options.cnt,options.tri, options.all, options.csv)
	c = time.process_time()
		
	

if __name__ == "__main__":
	main()


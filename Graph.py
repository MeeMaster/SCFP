def prepareRing(atoms):
	pass

class Graph(object):
	def __init__(self):
		self.graph={}
		self.circles = []
		self.gotPaths=False
		self.paths={'End':[],'Ring':[],'Complete':[]}
		self.truncated = False
		self.DLSFoundStatus=False
		self.DLSresult = None
	def addNode(self, node, edges=None):
		if not node in self.graph.keys():
			self.graph[node]=set([])
			if edges != None:
				#eligible = set(edges) & set(self.nodes())
				self.graph[node]=set(edges)#eligible
			else: 
				self.graph[node]=set([])
		else:
			pass
			#print 'Node already exists!'
				
	def addEdge(self, node1, node2, oneWay=False):
		if node1 in self.graph.keys() and node2 in self.graph.keys():
				if oneWay==False:
					self.graph[node1] | set(node2)
					self.graph[node2] | set(node1)
				else:
					self.graph[node1] | set(node2)
	
	def removeSelfEdges(self):
		for node in self.nodes():
			if node in self.graph[node]:
				self.graph[node]=self.graph[node] - set([node])
	
	def truncate(self):
		for node in self.nodes():
			self.graph[node] = (set(self.nodes()) & self.graph[node])
		self.truncated = True
		
	def deleteNode(self,node):
		if node in self.nodes():
			del self.graph[node]
	def size(self):
		return len(self.nodes())
	def nodes(self):
		return self.graph.keys()
	
	def neighbors(self, node):
		if node in self.nodes():
			return list(self.graph[node])
		else: return []
	
	def edges(self):
		inboundEdges=[]
		outboundEdges=[]
		verts = self.nodes()
		for vert in verts:
			for node in self.graph[vert]:
				if node in verts:
					if not set([vert,node]) in inboundEdges:
						inboundEdges.append(set([vert,node]))
				else:
					if not set([vert,node]) in outboundEdges:
						outboundEdges.append(set([vert,node]))
		return map(tuple,inboundEdges),map(tuple,outboundEdges)
		
	def intersect(self, queryGraph):
		if len(set(self.nodes()) & set(queryGraph.nodes())) > 0:
			return set(self.nodes()) & set(queryGraph.nodes())
		else: return False
	
	#def merge(self, queryGraph):
	#	temp = self.graph.nodes()
	#	for node in queryGraph.nodes():
	#		if node not in temp:
	#			temp.append(node)
	#	return molGraph(temp)
		
	def isWithin(self, queryGraph):
		if set(self.nodes()) < set(queryGraph.nodes()):
			return True
		else: return False
	
	def contains(self, queryGraph):
		if set(self.nodes()) > set(queryGraph.nodes()):
			return True
		else: return False
	
	def equals(self, queryGraph):
		if set(self.nodes()) == set(queryGraph.nodes()):
			return True
		else: return False
	
	def expandNode(self, node, prevnode=None):
		if prevnode !=None:
			children = list(set(self.neighbors(node)) - set([prevnode]))
			#children.remove(prevnode)
		else:
			children = self.neighbors(node)
		return children
		
	def DLS(self, node, maxdepth, goal, path=[],prevnode=None):
		newpath=path[:]
		newpath.append(node)
		if node == goal:
			self.DLSFoundStatus=True
			self.DLSresult=newpath
		if self.DLSFoundStatus==True:
			pass
		elif maxdepth > 0:
			children = self.expandNode(node)
			for child in children:
				self.DLS(child,maxdepth-1,goal,newpath)
				
	def iteratedDLS(self, root, goal):
		self.DLSFoundStatus=False
		self.DLSresult=None
		depth = 0
		while True:
			self.DLS(root,depth,goal)
			if self.DLSFoundStatus==True:
				return self.DLSresult
			depth +=1

	def DFS(self, root, path=[], prevnode=None):
		children = self.expandNode(root, prevnode)
		if children != []:
			for child in children:
				if child not in path:
					newpath = path[:]
					newpath.append(child)
					self.DFS(child,newpath,prevnode=root)
				else:
					ringpath = path[:]
					ringpath.append(child)
					if not ringpath in self.paths['Ring']:
						self.paths['Ring'].append(ringpath)
		else:
			if not path in self.paths['End']:
				self.paths['End'].append(path)
			
	
	def getPaths(self):
		startNode = min(self.graph.items(), key = lambda x: x[1])[0]
		self.DFS(startNode,path=[startNode])
		self.gotPaths=True
	
	def checkIntegrity(self):
		temp = set(self.nodes())
		if self.gotPaths==False:
			self.getPaths()
		#temp = self.graph
		for a in self.paths['End']:
			temp = temp-set(a)
		if len(temp)>0:
			return False,temp
		else:
			return True, None
	
	def searchRings(self):
		if self.gotPaths==False:
			self.getPaths()
		for path in self.paths['Ring']:
			ring = set(path[path.index(path[-1]):])
			if not ring in self.circles:
				self.circles.append(ring)


	
	#	if type(input)==canvas.ChmMol:
	#		atoms = input.getAtoms()
			#for bond in bonds:
			#	atom1 = bond.atom1()
			#	atom2 = bond.atom2()
			#	if graph.has_key(atom1):
			#		graph[atom1].append(atom2)
			#	else:
			#		graph[atom1]=[atom2]
			#	if graph.has_key(atom2):
			#		graph[atom2].append(atom1)
			#	else:
			#		graph[atom2]=[atom1]
	#	elif type(input)==list:
	#		atoms = input
	#	for atom in atoms:
	#		if type(atom)!=canvas.ChmAtom:
	#			raise TypeError
	#		graph[atom]=set(atom.getNeighbors())

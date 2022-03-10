#########################################################################
# File Name: graph.py
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Tue 22 Sep 2020 04:57:12 PM AEST
#########################################################################
#!/bin/bash
import os
import sys
#import tools
import networkx as nx

from libprism.local import tools
#type(nx.Graph()) <class 'networkx.classes.graph.Graph'>
# add node or edge
# 'add_cycle', 'add_edge', 'add_edges_from', 'add_node', 'add_nodes_from', 'add_path', 'add_star', 'add_weighted_edges_from'

# adj and factory 
# 'adj', 'adjacency', 'adjlist_inner_dict_factory', 'adjlist_outer_dict_factory', 'edge_attr_dict_factory', 'graph_attr_dict_factory'
#'node_attr_dict_factory', 'node_dict_factory'

# basic 
# 'clear', 'copy', 'degree', 'edges',  'get_edge_data', 'graph',  'neighbors', 'node', 'nodes',  'order', 'remove_edge', 'remove_edges_from', 'remove_node', 'remove_nodes_from',  'size', 'subgraph', 

# stat infor
# 'number_of_edges', 'number_of_nodes',

# selfloop
# 'nodes_with_selfloops',  'number_of_selfloops', 'selfloop_edges',

# check edge node
# 'has_edge', 'has_node', 'is_directed', 'is_multigraph',

#unknown or unuse
# 'name', 'nbunch_iter', 'edge_subgraph', 'fresh_copy', 'to_directed', 'to_directed_class', 'to_undirected', 'to_undirected_class', 'update'

class myGraph(nx.classes.graph.Graph):
    #def __init__(self, G):
    #    self.graph = G
     
    def displayGraph(self):
        print ("graph info:*******************")
        print ("number of nodes: ", format(self.number_of_nodes(), ',') )
        print ("number of edges: ", format(self.number_of_edges(), ',') )
        num = nx.number_connected_components(self)
        print ("number of components: ", num)       
        #self.show_subGraph() 
        objScore = self.get_total_weight_in_self()
        print ("total weight score is", format(objScore, ',') )
        print ("******************************")
        return objScore, num
   
    def simple_show(self):
        print ("graph info:*******************")
        print ("number of nodes: ", format(self.number_of_nodes(), ',') )
        print ("number of edges: ", format(self.number_of_edges(), ',') )
        return

    def get_total_weight_in_self(self):
        objScore = self.size(weight='weight') # weight sum in g
        return objScore


    def get_subgraphs(self):

        #graphs = list(nx.connected_component_subgraphs(self)) # <=2.3 networkx

        graphs = list(self.subgraph(c) for c in nx.connected_components(self))
        #graphs= list( (self.subgraph(c)).copy() for c in nx.connected_components(self))
        return graphs

    def show_subGraph(self):
        
        graphs = self.get_subgraphs()
        for g in graphs:
            gNodesNum = g.number_of_nodes()
            print (gNodesNum, end =" ")
        print ("")   
        return graphs


    def get_weight(self, n1, n2): 

        if self.has_edge(n1, n2):
            w = self.get_edge_data(n1, n2)['weight']
        else:
            w = 0
        return w    


    def update_weight_for_conflict_edges(self):
        removeEdges = set()
        for n in list(self.nodes):
            neig = set() 
            symN = tools.get_symmetrical_node(n)
            for ele in list(self.neighbors(n) ):
                symEle = tools.get_symmetrical_node(ele)
                if symEle in neig: # symEle and ele both connect to n
                    w1, w2 = self.get_weight(n, ele), self.get_weight(n, symEle)
                    w3, w4 = self.get_weight(symN, ele), self.get_weight(symN, symEle)
                    if w1 == 'A' or w2 == 'A' or w3 == 'A' or w4 == 'A':
                        continue
                    if (w1+w4)==0 or (w2+w3) == 0:
                        continue
                    if w1+w4 < 0.5*(w2+w3):  # read error is 10%,   10+10<90+90(r); 10+10<11+11(w)
                        if (ele, n) not in removeEdges:
                            removeEdges.add( (n, ele) )    
                            removeEdges.add( (symN, symEle ) )
                            self.add_edge(n, ele, weight=0)
                            self.add_edge(symN, symEle, weight=0)
                            #  add one edge can cover origal edges

                            wtemp = (w1+w4)/2.0
                            self.add_edge(n, symEle, weight= max(w2-wtemp,0 ) )
                            self.add_edge(symN, ele, weight= max(w3-wtemp,0 ) )
                            
                            # for high coverage and high hete, can use following code  
                            #self.add_edge(n, symEle, weight=0)
                            #self.add_edge(symN, ele, weight=0)
                            
                    elif 0.5*(w1+w4) > w2+w3:
                        if (symEle, n) not in removeEdges:
                            removeEdges.add( (n, symEle) )
                            removeEdges.add( (symN, ele) )
                            self.add_edge(n, symEle, weight=0)
                            self.add_edge(symN, ele, weight=0)

                            wtemp = (w2+w3)/2.0
                            self.add_edge(n, ele, weight= max(w1-wtemp,0 ) )
                            self.add_edge(symN, symEle, weight= max(w4-wtemp,0 ) )

                            #way 2
                            #self.add_edge(n, ele, weight=0)
                            #self.add_edge(symN, symEle, weight=0)
                    else:
                        self.add_edge(n, ele, weight=0)
                        self.add_edge(symN, symEle, weight=0)
                        self.add_edge(n, symEle, weight=0)
                        self.add_edge(symN, ele, weight=0)
                    
                        if (ele, n) not in removeEdges:
                            removeEdges.add( (n, ele) )    
                            removeEdges.add( (symN, symEle ) )
                        if (symEle, n) not in removeEdges:
                            removeEdges.add( (n, symEle) )
                            removeEdges.add( (symN, ele) )
                if ele not in neig:
                    neig.add(ele)

        self.remove_edges(removeEdges)
        return


    def check_symmetrical_node(self): # using the neg of symmetical node

        removeNode = set()
        for n in list(self.nodes):
            degreeNode = nx.degree(self, n)
            if n[-2:]=='_0' and degreeNode >= 2:
                symN = tools.get_symmetrical_node(n)
                if self.has_node(symN):
                    symDegreeNode = nx.degree(self, symN)          
                    if symDegreeNode >=2 and symDegreeNode + degreeNode >= 18:
                        neig1 = set()
                        for ele in list(self.neighbors(n) ):
                            neig1.add(ele[:-2])
                        if n[:-2] in neig1:
                            neig1.remove(n[:-2])
                        neig2 = set()
                        for ele in list(self.neighbors(symN) ):
                            neig2.add(ele[:-2])
                        inter = neig1.intersection(neig2)
                        # y<x<3y n and symN degree can be quite different
                        if len(inter) == 0: 
                            # or abs(degreeNode-symDegreeNode)>max(20, 2*min(degreeNode, symDegreeNode)):
                            #print ("remove node because intersection ==0", n, symN, degreeNode, symDegreeNode, len(inter) )
                            removeNode.add(n)
                            removeNode.add(symN)
        print ("remove node number: ", len(removeNode) )
        self.remove_nodes(removeNode)
        return

    def get_node_degree_distribution(self):

        degreeCnt = {}
        for n in list(self.nodes):
            degree = nx.degree(self, n)
            if degree not in degreeCnt:
                degreeCnt[degree] = 0
            degreeCnt[degree] += 1
        sortedCnt = sorted(degreeCnt.items(), key=lambda d: d[1], reverse=True)

        print ("Debug: node degree distribution")
        print ("Debug:", sortedCnt[0:20])
        print ("Debug: highest node degree:", sorted(degreeCnt, reverse=True)[0:20])
        
        if sortedCnt[0][0] != 0:
            threshold = 3*sortedCnt[0][0] # degree may 0, 1
        else:
            threshold = 3*sortedCnt[2][0]
        return threshold, threshold/3 # high degree, normal
    
    def remove_high_degree_node(self, threshold):

        removeNode = set()
        for n in list(self.nodes):
            degreeNode = nx.degree(self, n)
            if degreeNode >threshold:
                removeNode.add(n)

        print ("remove degree more than %s,  node number: %s " % (threshold ,len(removeNode)))
        self.remove_nodes(removeNode)
        return

    def remove_nodes(self, nodeSet):

        for node in nodeSet:
            if self.has_node(node):
                self.remove_node(node)
        print ("")        
        return

    def remove_type_A_edges(self): 

        for (n1, n2) in list(self.edges):
            w = self.get_weight(n1, n2)
            if w == 'A': # remove edges between symmetrical nodes
                self.remove_edge(n1, n2)
        return        

    def add_type_A_edges(self): # add edge bewtween 1_0 and 1_1

    
        for n in list(self.nodes):
            symN = tools.get_symmetrical_node(n)
            if self.has_node(symN) == False:
                self.add_node(symN)
            
        edges = []
        #selfWeightList = []
        for node in list(self.nodes):
            if node[-2:] == '_0':
                symN = tools.get_symmetrical_node(node)
                
                #w =  self.get_weight(node, symN) 
                #if w > 0:
                    #selfWeightList.append( (w, node) )
                
                #assert (symN in self.nodes)
                edges.append((symN, node, 'A'))
        self.add_weighted_edges_from(edges)

        #print ( sorted(selfWeightList, reverse = True) )
        #sys.exit()
        print ("type A edges (between symmetrical node)  number:", len(edges))        
        print ("add weight 'A' edges (between symmetrical node) number of nodes: ", self.number_of_nodes())
        return #G
    def calculate_edge_weight_distribution(self):

        weightCnt = {}
        for (n1, n2) in list(self.edges):
            w = self.get_edge_data(n1, n2)['weight'] 
            if w == 'A':
                continue
            if w not in weightCnt:
                weightCnt[w] = 0
            weightCnt[w] += 1
        
        weightList = sorted(weightCnt.items(), key=lambda item:item[1], reverse=True)
        print ("weight distribution" , weightList[1:10] )
        #return sorted(weightCnt.items())
        return weightList

    def remove_edges_according_weight(self, w): # filter edges without symmetrical edge and weigt <= wTemp

        removeEdges = set()
        for (n1,n2) in list(self.edges):
            w1 = self.get_weight(n1, n2)
            if w1 == 'A':
                continue
            symN1 = tools.get_symmetrical_node(n1)
            symN2 = tools.get_symmetrical_node(n2)
            w2 = self.get_weight(symN1, symN2) 
            if w1 <= w and w2 == 0:
                removeEdges.add( (n1, n2)  )
            #if( (w1>w2 and w2>0 and w1/w2>20)
            #or (w2>w1 and w1>0 and w2/w1>20) ) :
            #    removeEdges.add( (n1, n2)  )
            #    removeEdges.add( (symN1, symN2)  )
        print ("remove edge number", len(removeEdges))
        self.remove_edges(removeEdges)
        return 
    
    def remove_edges(self, removeEdges):

        #print("debug remove edge number:", len(removeEdges) )
        rD = {} # rDegree
        for (ele1, ele2) in removeEdges:
            #print (ele1, ele2, end =' ')
            if self.has_edge(ele1, ele2):
                self.remove_edge(ele1, ele2)
                '''
                if ele1 in rD:
                    rD[ele1] =rD[ele1] + 1
                else:
                    rD[ele1] = 1

                if ele2 in rD:
                    rD[ele2]= rD[ele2] + 1
                else:
                    rD[ele2] = 1            
                '''    
        #print ("")
        #for ele in rD:
        #    if rD[ele]>10:
        #        print (ele, rD[ele])        
                
        return
    
    def calculate_objScore_in_groupList(self, group):

        scoreSum, cnt = 0, 0
        for ele in group:
            #cnt+=1
            #print ("group %s generate subgraph" % (cnt))
            temp = self.subgraph(ele) 
            score = temp.get_total_weight_in_self() 
            scoreSum += score
        return scoreSum   


    def calculate_objScore(self, g1, g2):

        objScore, confilictScore = 0, 0
        for (n1, n2) in list(self.edges):
            if (n1 in g1) and (n2 in g1):
                objScore += self.get_edge_data(n1, n2)['weight']
            elif (n1 in g2) and (n2 in g2):
                objScore += self.get_edge_data(n1, n2)['weight']
            else:
                confilictScore += self.get_edge_data(n1, n2)['weight']

        ##assert graph.size(weight="weight") == objScore + confilictScore 
        return objScore, confilictScore

    def calculate_objScore2(self, g1, g2):

        scoreSum = 0
        
        temp = self.subgraph(g1) 
        score = temp.get_total_weight_in_self() 
        scoreSum += score
        
        temp = self.subgraph(g2) 
        score = temp.get_total_weight_in_self() 
        scoreSum += score

        return scoreSum
    def calculate_objScore_in_groupMap(self, group):

        scoreSum, cnt = 0, 0
        for ele in group:
            #cnt+=1
            #print ("group %s generate subgraph" % (cnt))
            temp = self.subgraph(group[ele]) 
            score = temp.get_total_weight_in_self() 
            scoreSum += score
        return scoreSum   


def calculate_weight_sum_between_node_and_group(node, group, g):

    wSum = 0
    edges = []
    if node not in g.nodes(): 
        return wSum , edges
    if len(g.adj[node]) < len(group):
        for adj, w in g.adj[node].items():
            if adj in group:
                if node == adj:
                    continue
                wtemp = int (w['weight'])
                if wtemp > 0:
                    edges.append( (node, adj) )
                wSum += wtemp        
    else:
        for n in group:
            if n == node:
                continue
            wtemp = g.get_weight(node, n)
            if wtemp > 0:
                edges.append( (node, n) )
            wSum += wtemp        

    return wSum, edges


def calculate_weight_sum_between_group_and_group(group1, group2, g):
    wSum = 0
    edges = []
    ##assert isinstance(group1, set) and  isinstance(group2, set)
    if len(group1) < len(group2):
        for node in group1:
            wtemp, etemp = calculate_weight_sum_between_node_and_group(node, group2, g)
            wSum += wtemp  
            edges.extend( etemp )
    else:
        for node in group2:
            wtemp, etemp = calculate_weight_sum_between_node_and_group(node, group1, g)
            wSum += wtemp  
            edges.extend( etemp )
    return wSum,  edges

def get_one_group_nei(group, phasedNodes, graph, n):

    neigh = list()
    for node in group[n]:
        if node not in graph.nodes:
            continue
        d = nx.degree(graph, node)
        if d == 0:
            continue
        temp = list( graph.neighbors(node) )
        for l in temp:
            if l in group:
                if isinstance(group[l], set):
                    neigh.append( l )
                else: 
                    temp = group[l]
                    neigh.append( temp )
            elif l in phasedNodes:
                temp = phasedNodes[l]
                while temp not in group:
                    temp = phasedNodes[temp]
                while isinstance(group[temp], set) ==  False:
                    temp = group[temp]
                neigh.append( temp )
            else:
                #print ("in myGraph.get_group_neighbors unknow reason")
                #print (l)
                group[l]=set()
                group[l].add(l)
                neigh.append( l )
                symL = tools.get_symmetrical_node(l)
                group[ symL ] =set()
                group[symL].add(symLi)
    return neigh         
                #sys.exit()
def get_group_neighbors(group, phasedNodes, graph, n):
   
    neigh =[]
    if isinstance(group[n], set) == False:
        return neigh

    neigh = get_one_group_nei(group, phasedNodes, graph, n)
    
    symN1 = tools.get_symmetrical_node(n)
    neigh2 =get_one_group_nei(group, phasedNodes, graph, symN1)
    
    neigh.extend(neigh2)    
    res = set( neigh )
    if n in res:
        res.remove(n)

    if symN1 in res:
        res.remove(symN1)

    res_0 = set()
    for ele in res:
        if ele.endswith("_0"):
            res_0.add(ele)
        elif ele.endswith("_1"):
            sym = tools.get_symmetrical_node(ele)
            if sym not in res:
                res_0.add(sym)
    res_0 = sorted(res_0) 
    return sorted(res_0, key=lambda d: len(group[d]) ) 

#!/usr/bin/env python3.4

import os
import sys
import random
import copy
import pdb

from libprism.local import myGraph
from libprism.local import tools

def inital_group(graph):

    group = {}
    for node in list(graph.nodes):
        if node in group:
            continue
        symN = tools.get_symmetrical_node( node )
        group[node] = set()
        group[node].add(node)
        group[symN] = set()
        group[symN].add(symN)
    return group

def keep_group_set(group): 
    # remove merged node in group set

    #unphasedNode = []
    removeNode = []
    for node in group:
        if isinstance(group[node], set):
            temp = group[node]
             # 18 May
            if len(temp) == 1:
                #unphasedNode.append(node)
                removeNode.append(node)
                continue
        else:
            removeNode.append(node)  
    print ("[Debug] inital group size:", len(group))        
    for n in removeNode:
        group.pop(n)
    print ("[Debug] final group size:", len(group))        
    return #unphasedNode


def keep_group_set2(group, phasedNode={}): 
    # remove merged node in group set

    removeNode = []
    for node in group:
        if isinstance(group[node], set) == False:
            removeNode.append(node) 
            phasedNode[node] = group[node]
    print ("[Debug] inital group size, phased group size:", len(group), len(removeNode))        
    for n in removeNode:
        group.pop(n)
    print ("[Debug] final group size:", len(group))        
    return #unphasedNode


def remove_small_group(group): 
    # remove merged node in group set
    removeNode = []
    newNodes = []
    for node in group:
        if isinstance(group[node], set):
            temp = group[node]
            if len(temp) <= 10 and len(temp)>1:
                removeNode.append(node)
                newNodes.extend(list(temp) )
                continue
    print ("[Debug] inital group size:", len(group))        
    for n in removeNode:
        group.pop(n)
    print ("[Debug] final group size (remove small group) :", len(group))       

    for node in newNodes:
        ''' 
        # debug can be remove
        if node in group:
            #print (node, group[node])
            assert len(group[node]) == 1
        '''    
        if node not in group:
            group[node] = set()
            group[node].add(node)
    return #unphasedNode


def print_group_set(group): 
    # remove merged node in group set
 
    print ("group detail")  
    for node in group:
        if isinstance(group[node], set) == True:
            if len(group[node]) >1:
                print (node, len(group[node]))        
    
    return #unphasedNode


def get_phased_nodes(group): 

    phasedNode = set()
    for node in group:
        if isinstance(group[node], set):
            temp = group[node]
            if len(temp) > 1: # 18 May
                phasedNode = phasedNode.union( group[node] )
    return phasedNode

def get_group_set(group, n1, n2):
    
    r1, r2 = n1, n2 
    if isinstance(group[n1], set):
        g1 = group[n1] 
    else:
        r1 = group[n1] 
        g1 = group[ r1 ] 
    if isinstance(group[n2], set):
        g2 = group[n2] 
    else:
        r2 = group[n2]
        g2 = group[ r2 ] 
    assert isinstance(g1, set)
    assert isinstance(g2, set)
    return g1, g2, r1, r2

def transfer_group_map_into_list(groupMap):
    

    finalGroups = []
    for node in groupMap:
        #assert isinstance(groupMap[node], set)
        temp = groupMap[node]
        #assert len(temp) > 1
        finalGroups.append( sorted(temp) )

    groupList = sorted(finalGroups, key=lambda ele: (len(ele), ele) )

    return groupList


def update_two_groups(r1, r2,  symR1, symR2, group, conflictSet):

     
    if (r1, r2) in conflictSet or (r2, r1) in conflictSet:
        #print ("(1) already in conflictSet")
        return False
  
    if check_groups_conflict(group[r1], group[r2]):
        #print ("add in conflictSet")
        conflictSet.add( (r1, r2) )
        conflictSet.add( (symR1, symR2) )
        return False
        
    lLabel = update_group(group, r1, r2, conflictSet)
    symLabel = update_group(group, symR1, symR2, conflictSet)
    return True

def update_group(group, n1, n2, conflictSet):

    #print (n1,n2)
    #print (group[n1], group[n2])
    #assert isinstance(group[n1], set) and isinstance(group[n2], set)

    #assert n2 != n1        
    if len(group[n1]) < len(group[n2]): # reduce running time
        temp = n1
        n1 = n2
        n2 = temp
    group[n1] = group[n1].union(group[n2]) # union may change id

    if len(group) < len(group[n2]):
        for ele in group:
            if ele != n1 and ele != n2 and ele in group[n2]:
                group[ele] = n1
    else:            
        for ele in group[n2]:
            if ele != n1 and ele != n2 and ele in group:
                group[ele] = n1
    group[n2] = n1
    
    return True


def check_groups_conflict(g1, g2):
    
    #print ("check g1, g2", g1, g2)
    l = []
    if len(g1) <= len(g2): #for reduce time
        for node in g1:
            l.append( tools.get_symmetrical_node(node) )
        lset = set(l)

        if len(lset.intersection(g2)) > 0: # conflict
            #print ("groups conflict", g1, g2) 
            return True
    else:
        for node in g2:
            l.append( tools.get_symmetrical_node(node) )
        lset = set(l)

        if len(lset.intersection(g1)) > 0: # conflict
            #print ("groups conflict", g1, g2) 
            return True
    return False


def check_group(group):
    for node in group:
        assert node.count('_') == 1
        if isinstance(group[node], set):
            for n in group[node]:
                assert n.count('_') == 1
    return True          


def write_groupList(groups, filename):
    
    fout = open(filename, "w")
    ID = 0
    groupsLen = len(groups)
    for j in range(0, groupsLen-1, 2):
        if len(groups[j]) != len(groups[j+1]):
            print ("two group lengths not equal")
            print (groups[j])
            print (groups[j+1])
            sys.exit()

        fout.write('group_%s_A_groupsize_%s\n' % (ID, len(groups[j])) )
        for node in groups[j]:
            fout.write(node+" ")
        fout.write("\n")    

        fout.write('group_%s_B_groupsize_%s\n' % (ID, len(groups[j+1] ) ) )
        for node in groups[j+1]:
            fout.write(node+" ")
        fout.write("\n")
        ID += 1
        
    fout.close()
    return 

def get_groups_sum(n1, n2, symN1, symN2, group, graph, groupSum):
    
    #g1, g2 = group[n1], group[n2]
    g1,g2, r1, r2 = get_group_set(group, n1, n2)
    l1, l2 = len(g1), len(g2) 
    if (n1, n2, l1, l2) in groupSum:
        w1, w2, w3, w4, e1, e2, e3, e4 = groupSum[ (n1, n2, l1, l2)  ]
    else:
        symG1, symG2 = group[ symN1 ], group[ symN2 ]
        w1, e1 = myGraph.calculate_weight_sum_between_group_and_group(g1, g2, graph)
        w2, e2 = myGraph.calculate_weight_sum_between_group_and_group(g1, symG2, graph)
        w3, e3 = myGraph.calculate_weight_sum_between_group_and_group(symG1, g2, graph)
        w4, e4 = myGraph.calculate_weight_sum_between_group_and_group(symG1, symG2, graph)
        groupSum[ (n1, n2, l1, l2)  ] = w1, w2, w3, w4, e1, e2, e3, e4

    return w1, w2, w3, w4, l1, l2, e1, e2, e3, e4 


def is_two_nodes_merged(n1, n2, symN1, symN2, group):
    
    if isinstance(group[n1], set) and (n2 in group[n1]) :
        return True
    if isinstance(group[n2], set) and (n1 in group[n2]) :
        return True
    if group[n1] == group[n2] :
        return True
    return False

# merge group

def merge_two_groups(n1, n2, symN1, symN2, group, conflictSet):
    
    if is_two_nodes_merged(n1, n2, symN1, symN2, group):
        return False  
    g1, g2, r1, r2 = get_group_set(group, n1, n2)
    symG1, symG2, symR1, symR2 = get_group_set(group, symN1, symN2)
    label = update_two_groups(r1, r2, symR1, symR2, group, conflictSet)
    
    return label


def merge_group_according_threshold_strict(graph, group, phasedNodes, groupSum, conflictSet, w, label):
    
    temp = sorted( sorted(group.items()), key=lambda item:len(item[1]))
    #index = tools.binarySearch_ListValue(temp, 1) # first group size >1
    groupKeys = [ele[0] for ele in temp]
    groupLen = len(group)
    #print ("for order", groupKeys) 
    for i in range(groupLen-1):
    #for i in range(index, groupLen-1):
        if groupKeys[i].endswith('_0'):
            continue 
        n1 = groupKeys[i]
        if n1 not in group: # or n1 not in graph.nodes:
            continue    
        #assert n1 in group
        l = myGraph.get_group_neighbors(group, phasedNodes, graph, n1)
        for n2 in l:
            symN1 = tools.get_symmetrical_node( n1 )
            symN2 = tools.get_symmetrical_node( n2 ) 
            w1, w2, w3, w4, l1, l2, e1, e2, e3, e4 = get_groups_sum(n1, n2, symN1, symN2, group, graph, groupSum)
            if l1 <=1 and l2 <= 1:
                continue
            if label==True: # label == True for realDara
                if l1 <=1 or l2 <= 1:
                    continue
            s1 = w1 + w4 - w2 - w3
            s2 = -1*s1 #w2 + w3 - w1 - w4
            rate1 = float(s1) / (l1 + l2)
            rate2 = float(s2) / (l1 + l2)
            if ( s1> max(0.2*min(l1, l2), 5) and w1>2 and w4>2 and rate1>w 
                and (w2+w3)<=min(20, (w1+w4)*0.1)  #min( 2, (w1+w4)*0.2 ) 
                and max(w1,w4)/min(w1,w4)<10 ):
                r1 = (w1+w4)/(len(e1)+len(e4))
                if r1>1.0:
                    merge_two_groups(n1, n2, symN1, symN2, group, conflictSet)
                    groupSum.pop( (n1, n2, l1, l2) )
                    #print ("[Debug] merge strict", n1, n2,  w1, w2, w3, w4, l1, l2)
            elif (s2> max(0.2*min(l1, l2), 5) and w2>2 and w3>2 and rate2>w
                   and (w1+w4)<=min(20, (w2+w3)*0.1) #min( 2, (w2+w3)*0.2 ) 
                   and max(w2,w3)/min(w2,w3)<10 ):
                r2 = (w2+w3)/(len(e2)+len(e3))
                if r2>1.0:
                    merge_two_groups(n1, symN2, symN1, n2, group, conflictSet)
                    groupSum.pop( (n1, n2, l1, l2) )
                    #print ("[Debug] merge strict", n1, n2,  w1, w2, w3, w4, l1, l2)
            else:
                continue
            #break
    keep_group_set2(group, phasedNodes)
    label = True
    if len(group) == groupLen:
        label = False
    return label #, set(newKeys) 

def merge_node_group_strict(graph, group, phasedNodes, groupSum, conflictSet, w): #, groupKeys):
   
    temp = sorted( sorted(group.items()), key=lambda item:len(item[1]), reverse=True)
    #index = tools.binarySearch_ListValue(temp, 1) # first group size >1
    groupKeys = [ele[0] for ele in temp]
    groupLen = len(groupKeys)
    
    for i in range(groupLen-1):
    #for i in range(index, groupLen-1):
        if groupKeys[i].endswith('_0'):
            continue 
        n1 = groupKeys[i]
        if n1 not in group: #or n1 not in graph.nodes:
            continue
        l = myGraph.get_group_neighbors(group, phasedNodes, graph, n1)    
        for n2 in l:
            symN1 = tools.get_symmetrical_node( n1 )
            symN2 = tools.get_symmetrical_node( n2 ) 
            w1, w2, w3, w4, l1, l2, e1, e2, e3, e4 = get_groups_sum(n1, n2, symN1, symN2, group, graph, groupSum)
            
            if (l1 <=1 and l2 <= 1) or (l1>2 and l2>2):
                continue
            s1 = w1 + w4 - w2 - w3
            s2 = -1*s1
            rate1 = float(s1) / (l1 + l2)
            rate2 = float(s2) / (l1 + l2)
            if (rate1>w and (w2+w3)<=min(20, (w1+w4)*0.1)
                    and s1 >5 and w1>2 and w4>2): #min( 2, (w1+w4)*0.2 ) ):
                r1 = (w1+w4)/(len(e1)+len(e4))
                if (r1>1.0):
                    merge_two_groups(n1, n2, symN1, symN2, group, conflictSet)
                    groupSum.pop( (n1, n2, l1, l2) )
                    #print ("[Debug] merge node group", n1, n2,  w1, w2, w3, w4, l1, l2, r1)
                    #if w2+w3 > 0 :
                    #    graph.remove_edges(e2+e3)
            elif (rate2>w and (w1+w4)<=min(20, (w2+w3)*0.1)
                      and s2 >5 and w2>2 and w3>2):  #min( 2, (w2+w3)*0.2 ) ):
                r2 = (w2+w3)/(len(e2)+len(e3))
                if( r2>1.0):
                    merge_two_groups(n1, symN2, symN1, n2, group, conflictSet)
                    groupSum.pop( (n1, n2, l1, l2) )
                    #print ("[Debug] merge node group", n1, n2,  w1, w2, w3, w4, l1, l2,  r2)
                    #if w1+w4 > 0 :
                    #    graph.remove_edges(e1+e4)
            else:
                continue

    #print ("3", group)
    keep_group_set2(group, phasedNodes)
    #print ("4", group)
    label = True
    if len(group) >= groupLen:
        label = False
    return label #, set(newKeys) 





def merge_group_according_threshold(graph, group, phasedNodes, threshold, groupSum, conflictSet, label):
     
    temp = sorted(sorted(group.items()), key=lambda item:len(item[1]))
    index = tools.binarySearch_ListValue(temp, 1) # first group size >1
    groupKeys = [ele[0] for ele in temp]
    groupLen = len(group)

    #for i in range(groupLen-1):
    for i in range(index, groupLen-1):
        if groupKeys[i].endswith('_1') == False:
            continue 
        n1 = groupKeys[i]
        if n1 not in group:# or n1 not in graph.nodes:
            continue   
        l = myGraph.get_group_neighbors(group, phasedNodes, graph, n1)
        for n2 in l:
            symN1 = tools.get_symmetrical_node( n1 )
            symN2 = tools.get_symmetrical_node( n2 ) 
            w1, w2, w3, w4, l1, l2, e1, e2, e3, e4 = get_groups_sum(n1, n2, symN1, symN2, group, graph, groupSum)
            rate = 0.02
            if label == True:   
                if l1 <=1 or l2 <= 1:
                    continue
                rate = 0.085
            s1 = w1 + w4 - w2 - w3
            s2 = -1*s1
            if (s1 > max(threshold*min(l1, l2), max(5, rate*(l1+l2)) ) and w1>2 and w4>2   
                 and (w2+w3)<min(100, 0.2*(w1+w4)) and min(w1,w4) >max(w2,w3)
                 and max(w1,w4)/min(w1,w4)<5 ):
                r1 = (w1+w4)/(len(e1)+len(e4))
                if r1>1.0:
                    merge_two_groups(n1, n2, symN1, symN2, group, conflictSet)
                    print ("[Debug] merge group", n1, n2,  w1, w2, w3, w4, l1, l2, r1)
                    groupSum.pop( (n1, n2, l1, l2) )
                    #if w2+w3 > 0 :
                        #graph.remove_edges(e2+e3)
            if (s2 > max(threshold*min(l1, l2), max(5, rate*(l1+l2)) ) and w2>2 and w3>2 
                  and (w1+w4)<min(0.2*(w2+w3), 100) and min(w2, w3) >max(w1, w4) 
                  and max(w2,w3)/min(w2,w3)<5 ): 
                r2 = (w2+w3)/(len(e2)+len(e3))
                if r2>1.0:
                    merge_two_groups(n1, symN2, symN1, n2, group, conflictSet)
                    print ("[Debug] merge group", n1, n2,  w1, w2, w3, w4, l1, l2, r2 )
                    groupSum.pop( (n1, n2, l1, l2) )
                    #if w1+w4 > 0 :
                        #graph.remove_edges(e1+e4)
    keep_group_set2(group, phasedNodes)
    label = True
    if len(group) == groupLen:
        label = False
    return label 


def merge_groups_only_choice(graph, group, phasedNodes, groupSum, conflictSet):
 
    temp = sorted( sorted(group.items()), key=lambda item:len(item[1]))
    groupKeys = [ele[0] for ele in temp]
    groupLen = len(group)
    print ("before merge only choice group: %s pairs of group" %  (groupLen/2) )

    onlyChoice = {}
    for i in range(groupLen-1):
        if groupKeys[i].endswith('_0'):
            continue
        n1 = groupKeys[i]
        if n1 not in group:# or n1 not in graph.nodes:
            continue    
        l = myGraph.get_group_neighbors(group, phasedNodes, graph, n1)
        for n2 in l:
            symN1 = tools.get_symmetrical_node( n1 )
            symN2 = tools.get_symmetrical_node( n2 ) 
            w1, w2, w3, w4, l1, l2, e1, e2, e3, e4 = get_groups_sum(n1, n2, symN1, symN2, group, graph, groupSum)
            if w1+w2+w3+w4 == 0:
                continue
            if n1 not in onlyChoice:
                onlyChoice[n1] = []
            if n2 not in onlyChoice:
                onlyChoice[n2] = []
            
            onlyChoice[n1].append( (n1, n2, symN1, symN2, w1, w2, w3, w4, l1, l2, e1, e2, e3, e4) )
            onlyChoice[n2].append( (n1, n2, symN1, symN2, w1, w2, w3, w4, l1, l2, e1, e2, e3, e4) )
            
    for n in onlyChoice:
        if len(onlyChoice[n]) == 1:
            n1, n2, symN1, symN2, w1, w2, w3, w4, l1, l2, e1,e2,e3,e4 = onlyChoice[n][0]
            if isinstance(group[n1], set) == False or isinstance(group[n2], set) == False:
                continue 
            s1 = w1 + w4 - w2 - w3
            s2 = (-1)*s1
            if (s1 > max(2, 0.01*(l1+l2)) and (w1>0 and w4>0)  #or (w2+w3==0 )) 
                    and w2+w3<min(100,0.2*(w1+w4) ) ):
                if len(group[n1]) >10 and len(group[n2])>10:
                    if w1==0 or w4==0 or max(w1,w4)/min(w1,w4)>5:
                        continue 
                r1 = (w1+w4)/(len(e1)+len(e4))
                if r1>1.0:
                    merge_two_groups(n1, n2, symN1, symN2, group, conflictSet)
                    #if w2+w3 > 0 :
                    #    graph.remove_edges(e2+e3)
                    print ("[Debug] merge only choice", n1, n2,  w1, w2, w3, w4, l1, l2)
                    groupSum.pop( (n1, n2, l1, l2) )
            elif (s2 > max(2, 0.01*(l1+l2)) and (w2>0 and w3>0) # or (w1+w4==0)) 
                      and w1+w4<min(100,0.2*(w2+w3) ) ):

                if len(group[n1]) >10 and len(group[n2])>10: #give small group more chance
                    if w2==0 or w3 ==0 or max(w2,w3)/min(w2,w3)>5:
                        continue 

                r2 = (w2+w3)/(len(e2)+len(e3))
                if r2>1.0:
                    merge_two_groups(n1, symN2, symN1, n2, group, conflictSet)
                    #if w1+w4 > 0 :
                    #    graph.remove_edges(e1+e4)
                    print ("[Debug] merge only choice", n1, n2,  w1, w2, w3, w4, l1, l2)
                    groupSum.pop( (n1, n2, l1, l2) )
    
    keep_group_set2(group, phasedNodes)
    label = True
    if len(group) == groupLen:
        label = False
    return label


def show_group_connection(group, graph, groupSum):
   
    print ("show group connection") 
    temp = sorted( sorted(group.items()), key=lambda item:len(item[1]))
    groupKeys = [ele[0] for ele in temp]
    groupLen = len(group)
    groupDegree = {}
    #onlyChoice = {}
    vipPair= [] 
    for i in range(groupLen-1):
        for j in range(i+1, groupLen):
            if groupKeys[i].endswith('_1') and groupKeys[j].endswith('_1'):

                n1, n2 = groupKeys[i], groupKeys[j] 
                symN1 = tools.get_symmetrical_node( n1 )
                symN2 = tools.get_symmetrical_node( n2 ) 
                w1, w2, w3, w4, l1, l2, e1, e2, e3, e4 = get_groups_sum(n1, n2, symN1, symN2, group, graph, groupSum)
                if l1 <=1 or l2 <=1:
                    continue
                if w1+w2+w3+w4 == 0:
                    continue
                
                if n1 not in groupDegree:
                    groupDegree[n1] = 0
                groupDegree[n1] += 1
                if n2 not in groupDegree:
                    groupDegree[n2] = 0
                groupDegree[n2] += 1
               
                if (w2+w3==0) or (w1+w4==0):
                    vipPair.append((n1, n2))
        
                if w1+w2+w3+w4>1:
                    if (w1+w4)>(w2+w3):
                        print ("[Debug] connection", n1, n2,  w1, w2, w3, w4, l1, l2, (w1+w4)/(len(e1)+len(e4) ) )
                    if (w2+w3)>(w1+w4):
                        print ("[Debug] connection", n1, n2,  w1, w2, w3, w4, l1, l2, (w2+w3)/(len(e2)+len(e3) ) )

    #temp = sorted(groupDegree.items(), key=lambda item:item[1], reverse=True)
    #print ("groupDegree", temp)
    return vipPair, groupDegree            

def merge_group_according_threshold_loose_only(graph, group, phasedNodes, groupSum, conflictSet):
    
    temp = sorted(group.items(), key=lambda item:len(item[1]))
    groupKeys = [ele[0] for ele in temp]
    groupLen = len(group)

    #print ("merge group loose")
    onlyChoice = {}
    for i in range(groupLen-1):
        if groupKeys[i].endswith('_0'):
            continue 
        n1 = groupKeys[i]
        if n1 not in group:
            continue    
        l = myGraph.get_group_neighbors(group, phasedNodes, graph, n1)
        for n2 in l:
            symN1 = tools.get_symmetrical_node( n1 )
            symN2 = tools.get_symmetrical_node( n2 ) 
            w1, w2, w3, w4, l1, l2, e1, e2, e3, e4 = get_groups_sum(n1, n2, symN1, symN2, group, graph, groupSum)
            s1 = w1 + w4 - w2 - w3
            s2 = w2 + w3 - w1 - w4
            if ( s1 >abs(w1-w4)+abs(w2-w3) 
                    and ( (w1>1 and w4>1) and (w2+w3==0 or (w1+w4)/(w2+w3)>3) )  ): 
                rate2 = float(s1) / (l1 + l2)
                if rate2 < 0.01:
                    continue
                if n1 not in onlyChoice:
                    onlyChoice[n1] = []
                if n2 not in onlyChoice:
                    onlyChoice[n2] = []
                
                onlyChoice[n1].append( (n1, n2, symN1, symN2, w1, w2, w3, w4, l1, l2) )
                onlyChoice[n2].append( (n1, n2, symN1, symN2, w1, w2, w3, w4, l1, l2) )
            elif  ( s2> abs(w1-w4)+abs(w2-w3) 
                    and ( (w2>1 and w3>1) and (w1+w4==0 or (w2+w3)/(w1+w4)>3 ) ) ):
                rate2 = float(s2) / (l1 + l2)
                if rate2 < 0.01:
                    continue

                if n1 not in onlyChoice:
                    onlyChoice[n1] = []
                if n2 not in onlyChoice:
                    onlyChoice[n2] = []
                
                onlyChoice[n1].append( (n1, n2, symN1, symN2, w1, w2, w3, w4, l1, l2) )
                onlyChoice[n2].append( (n1, n2, symN1, symN2, w1, w2, w3, w4, l1, l2) )

    for n in onlyChoice:
        if len(onlyChoice[n]) == 1:
            n1, n2, symN1, symN2, w1, w2, w3, w4, l1, l2 = onlyChoice[n][0]
            if len(onlyChoice[n1]) > 1 or len(onlyChoice[n2]) >1:
                continue
            if isinstance(group[n1], set) == False or isinstance(group[n2], set) == False:
                continue   
            if (w1 + w4 - w2 - w3 > 0):                
                merge_two_groups(n1, n2, symN1, symN2, group, conflictSet)
                print ("[Debug] merge loose only", n1, n2,  w1, w2, w3, w4, l1, l2)
                groupSum.pop( (n1, n2, l1, l2) )
            elif (w2 + w3 - w1 - w4 > 0):
                merge_two_groups(n1, symN2, symN1, n2, group, conflictSet)
                print ("[Debug] merge loose only", n1, n2,  w1, w2, w3, w4, l1, l2)
                groupSum.pop( (n1, n2, l1, l2) )
    
    keep_group_set2(group, phasedNodes)
    label = True
    if len(group) == groupLen:
        label = False
    return label 


def merge_group_final_loose(graph, group, groupSum, conflictSet, threshold):
    
    temp = sorted(group.items(), key=lambda item:len(item[1]))
    index = tools.binarySearch_ListValue(temp, 1) # first group size >1
    #print (index, temp[index-1:index+2])
    groupKeys = [ele[0] for ele in temp]
    groupLen = len(group)

    print ("merge group final loose")
    for i in range(groupLen-1):
        #for j in range(max(i+1, index), groupLen):
        for j in range(i+1, groupLen):
            if groupKeys[i].endswith('_1') and groupKeys[j].endswith('_1'):

                n1, n2 = groupKeys[i], groupKeys[j] 
                symN1 = tools.get_symmetrical_node( n1 )
                symN2 = tools.get_symmetrical_node( n2 ) 
                w1, w2, w3, w4, l1, l2, e1, e2, e3, e4 = get_groups_sum(n1, n2, symN1, symN2, group, graph, groupSum)
                if l1>10000 or l2>10000:
                    continue
                s1 = w1 + w4 - w2 - w3
                s2 = w2 + w3 - w1 - w4
                rate1 = float(s1) / (l1 + l2)
                rate2 = float(s2) / (l1 + l2)
                #k=31 ok
                #k=21 a lot of err
                if (s1 > 2 and w1>1 and w4>1
                       and  min(w1,w4) > max(w2,w3) and rate1> threshold and max(w1,w4)/min(w1,w4)<4):
                    if ( ( (w2+w3)<0.1*(w1+w4) ) #and s1>(l1+l2)*0.02 )
                    #or   ( (w2+w3)<0.3*(w1+w4) and (rate1>2*threshold or s1>=max(0.5*min(l1,l2), 3)) ) 
                    or ( w1>2 and w4>2 and w2+w3==0 ) 
                    or (min(l1,l2)==1 and l1+l2>2) ):
                        
                        r1 = (w1+w4)/(len(e1)+len(e4))
                        if r1>1.2: # or (min(l1,l2)<=10 and w2+w3==0 and s1>=max(3, 0.5*min(l1,l2))):
                            merge_two_groups(n1, n2, symN1, symN2, group, conflictSet)
                            print ("[Debug] merge final group", n1, n2,  w1, w2, w3, w4, l1, l2, r1)
                            groupSum.pop( (n1, n2, l1, l2) )
                elif (s2 > 2 and w2>1 and w3>1
                        and min(w2,w3) > max(w1,w4) and rate2>threshold and max(w2,w3)/min(w2,w3)<4 ):
                    if ( ( (w1+w4)<0.1*(w2+w3) ) #and s2>(l1+l2)*0.02)
                    #or   ( (w1+w4)<0.3*(w2+w3) and (rate2>2*threshold or s1>=max(0.5*min(l1,l2), 3)) ) 
                    or ( w2>2 and w3>2 and w1+w4==0) 
                    or (min(l1, l2)==1 and l1+l2>2) ):

                        r2 = (w2+w3)/(len(e2)+len(e3))
                        if r2>1.2: #or (min(l1,l2)<=10 and w1+w4==0 and s2>=max(3, 0.5*min(l1,l2)) ) :
                            merge_two_groups(n1, symN2, symN1, n2, group, conflictSet)
                            print ("[Debug] merge final group", n1, n2,  w1, w2, w3, w4, l1, l2,r2)
                            groupSum.pop( (n1, n2, l1, l2) )
    keep_group_set2(group)
    label = True
    if len(group) == groupLen:
        label = False
    return label 


def merge_group_on_reliable_graph(graph, group, phasedNodes, groupSum, conflictSet, threshold, l):
    
    temp = sorted( sorted(group.items()), key=lambda item:len(item[1]))
    #index = tools.binarySearch_ListValue(temp, 1) # first group size >1
    groupKeys = [ele[0] for ele in temp]
    groupLen = len(group)

    print ("merge group on reliable graph")
    for i in range(groupLen-1):
        #for j in range(max(i+1, index), groupLen):
        for j in range(i+1, groupLen):
            if groupKeys[i].endswith('_1') and groupKeys[j].endswith('_1'):

                n1, n2 = groupKeys[i], groupKeys[j] 
                symN1 = tools.get_symmetrical_node( n1 )
                symN2 = tools.get_symmetrical_node( n2 ) 
                w1, w2, w3, w4, l1, l2, e1, e2, e3, e4 = get_groups_sum(n1, n2, symN1, symN2, group, graph, groupSum)
                if l1>10000 or l2>10000: 
                    continue
                if l == True:
                    if l1<=1 or l2<=1:
                        continue
                s1 = w1 + w4 - w2 - w3
                s2 = w2 + w3 - w1 - w4
                rate1 = float(s1) / (l1 + l2)
                rate2 = float(s2) / (l1 + l2)
                if (s1 >= 3 and w1>0 and w4>0 and 
                        (w2+w3==0 or (w1+w4)/(w2+w3)>5)
                        and rate1>threshold):
                    if min(l1, l2)>10 and min(w1,w4)<=1:
                        continue

                    r1 = (w1+w4)/(len(e1)+len(e4))
                    if r1>1.0:
                        merge_two_groups(n1, n2, symN1, symN2, group, conflictSet)
                        print ("[Debug] merge reliable graph", n1, n2,  w1, w2, w3, w4, l1, l2, r1)
                        #print (len(e1), len(e2), len(e3), len(e4))
                        groupSum.pop( (n1, n2, l1, l2) )
                elif (s2 >= 3 and w2>0 and w3>0 and 
                        (w1+w4==0 or (w2+w3)/(w1+w4)>5)
                            and rate2>threshold):
                    #(s1>(l1+l2)*0.01 or s1>=max(0.5*min(l1,l2), 3))
                    if min(l1, l2)>10 and min(w2,w3)<=1:
                        continue
                    r2 = (w2+w3)/(len(e2)+len(e3))
                    if r2>1.0:
                        merge_two_groups(n1, symN2, symN1, n2, group, conflictSet)
                        print ("[Debug] merge reliable graph", n1, n2,  w1, w2, w3, w4, l1, l2,r2)
                        #print (len(e1), len(e2), len(e3), len(e4))
                        groupSum.pop( (n1, n2, l1, l2) )
    keep_group_set2(group, phasedNodes)
    label = True
    if len(group) == groupLen:
        label = False
    return label 


def merge_group_on_original_graph(graph, group, phasedNodes, groupSum, conflictSet, threshold): # groupDegree):
    
    temp = sorted( sorted(group.items()), key=lambda item:len(item[1]))
    groupKeys = [ele[0] for ele in temp]
    groupLen = len(group)

    print ("merge group on original graph")
    for i in range(groupLen-1):
        conn={}
        for j in range(i+1, groupLen):
            if groupKeys[i].endswith('_1') and groupKeys[j].endswith('_1'):
                n1, n2 = groupKeys[i], groupKeys[j]
                symN1 = tools.get_symmetrical_node( n1 )
                symN2 = tools.get_symmetrical_node( n2 )
                w1, w2, w3, w4, l1, l2, e1, e2, e3, e4 = get_groups_sum(n1, n2, symN1, symN2, group, graph, groupSum) 
                s1 = w1 + w4 - w2 - w3
                if abs(s1) >0:
                    conn[(n1, n2, symN1, symN2)] = abs(s1)
        
        temp = sorted( sorted(conn.items()), key=lambda item:item[1], reverse=True)
        print ("debug", temp)
        pairKey= [ele[0] for ele in temp]
        for (n1, n2, symN1, symN2) in pairKey:       
                w1, w2, w3, w4, l1, l2, e1, e2, e3, e4 = get_groups_sum(n1, n2, symN1, symN2, group, graph, groupSum) 
                if l1>10000 or l2>10000 or l1<=1 or l2<=1:
                    continue
                s1 = w1 + w4 - w2 - w3
                s2 = w2 + w3 - w1 - w4
                rate1 = float(s1) / (l1 + l2)
                rate2 = float(s2) / (l1 + l2)
                if (s1 >= 3 and w1>2 and w4>2 
                        and min(w1,w4)>max(w2,w3)
                        and (w2+w3==0 or (w1+w4)/(w2+w3)>5)
                        and  ( (max(w1,w4)/min(w1,w4) <4 and rate1>threshold)
                            or s1>0.5*min(l1,l2) ) ):
                    r1 = (w1+w4)/(len(e1)+len(e4))
                    if r1>1.0:
                        merge_two_groups(n1, n2, symN1, symN2, group, conflictSet)
                        print ("[Debug] merge original graph", n1, n2,  w1, w2, w3, w4, l1, l2, r1)
                        #print (len(e1), len(e2), len(e3), len(e4))
                        groupSum.pop( (n1, n2, l1, l2) )
                elif (s2 >= 3 and w2>2 and w3>2
                        and min(w2,w3)>max(w1,w4)
                        and (w1+w4==0 or (w2+w3)/(w1+w4)>5)
                        and ( (max(w2,w3)/min(w2,w3) <4 and rate2>threshold) 
                        or s2>0.5*min(l1,l2) ) ):
                    #(s1>(l1+l2)*0.01 or s1>=max(0.5*min(l1,l2), 3))
                    r2 = (w2+w3)/(len(e2)+len(e3))
                    if r2>1.0:
                        merge_two_groups(n1, symN2, symN1, n2, group, conflictSet)
                        print ("[Debug] merge original graph", n1, n2,  w1, w2, w3, w4, l1, l2,r2)
                        #print (len(e1), len(e2), len(e3), len(e4))
                        groupSum.pop( (n1, n2, l1, l2) )
    keep_group_set2(group, phasedNodes)
    label = True
    if len(group) == groupLen:
        label = False
    return label 


def merge_vip_pair(graph, group, phasedNodes, groupSum, conflictSet, vipPair):
    
    temp = sorted( sorted(group.items()), key=lambda item:len(item[1]))
    #index = tools.binarySearch_ListValue(temp, 1) # first group size >1
    groupKeys = [ele[0] for ele in temp]
    groupLen = len(group)

    print ("merge group on original graph")
    #for i in range(groupLen-1):
        #for j in range(i+1, groupLen):
        #    if groupKeys[i].endswith('_1') and groupKeys[j].endswith('_1'):
    for (n1, n2) in vipPair:
        #n1, n2 = groupKeys[i], groupKeys[j]
        symN1 = tools.get_symmetrical_node( n1 )
        symN2 = tools.get_symmetrical_node( n2 )

        w1, w2, w3, w4, l1, l2, e1, e2, e3, e4 = get_groups_sum(n1, n2, symN1, symN2, group, graph, groupSum)
        if l1>10000 or l2>10000 or l1<=1 or l2<=1:
            continue
        #if l1 > 100 and l2 >100:
        #    continue
        s1 = w1 + w4 - w2 - w3
        s2 = w2 + w3 - w1 - w4
        rate1 = float(s1) / (l1 + l2)
        rate2 = float(s2) / (l1 + l2)
        if (s1 >= 3 and w1>0 and w4>0 
                and (w2+w3==0 or (w1+w4)/(w2+w3)>=5) ):
            r1 = (w1+w4)/(len(e1)+len(e4))
            if r1>1.0:
                merge_two_groups(n1, n2, symN1, symN2, group, conflictSet)
                print ("[Debug] merge vip pair", n1, n2,  w1, w2, w3, w4, l1, l2, r1)
            #print (len(e1), len(e2), len(e3), len(e4))
            groupSum.pop( (n1, n2, l1, l2) )
        elif (s2 >= 3 and w2>0 and w3>0
                and (w1+w4==0 or (w2+w3)/(w1+w4)>=5) ):
            r2 = (w2+w3)/(len(e2)+len(e3))
            if r2>1.0:
                merge_two_groups(n1, symN2, symN1, n2, group, conflictSet)
                print ("[Debug] merge vip pair", n1, n2,  w1, w2, w3, w4, l1, l2,r2)
            groupSum.pop( (n1, n2, l1, l2) )
        #else:  
        #    print ("[Debug] not merge vip pair", n1, n2,  w1, w2, w3, w4, l1, l2)
    keep_group_set2(group, phasedNodes)
    label = True
    if len(group) == groupLen:
        label = False
    return label 



#useless
def merge_group_according_decision_tree(graph, group, groupSum, conflictSet, cfg):
   

    temp = sorted(group.items(), key=lambda item:len(item[1]))
    groupKeys = [ele[0] for ele in temp]
    index = tools.binarySearch_ListValue(temp, 1) # first group size >1
    #print (index, temp[index-1:index+2])
    
    groupLen = len(group)
    print ("merge group according decision tree")
    print ("[Debug] group size", len(group))

    #print (group.keys())
    #print_group_set(group)
    for i in range(groupLen-1):
        for j in range(max(i+1, index), groupLen):
            if groupKeys[i].endswith('_1') and groupKeys[j].endswith('_1'):
                n1, n2 = groupKeys[i], groupKeys[j] 
                symN1 = tools.get_symmetrical_node( n1 )
                symN2 = tools.get_symmetrical_node( n2 ) 
                w1, w2, w3, w4, l1, l2, e1, e2, e3, e4 = get_groups_sum(n1, n2, symN1, symN2, group, graph, groupSum)
                temp = [w1, w2, w3, w4, l1, l2]
                #pdb.set_trace()
                if w1+w2+w3+w4 ==0:
                    continue
                res = cfg.predict( [temp] )
                
                if (res[0] == 1):
                    print (temp)
                
                if (res[0] == 1 and  w1 + w4 - w2 - w3 > abs(w1-w4)+abs(w2-w3) ):
                    merge_two_groups(n1, n2, symN1, symN2, group, conflictSet)
                    print ("[Debug] merge group", n1, n2,  w1, w2, w3, w4, l1, l2)
                    check_group(group)
                    groupSum.pop( (n1, n2, l1, l2) )
                elif (res[0] == 1 and w2 + w3 - w1 - w4 >abs(w1-w4)+abs(w2-w3) ):
                    merge_two_groups(n1, symN2, symN1, n2, group, conflictSet)
                    print ("[Debug] merge group", n1, n2,  w1, w2, w3, w4, l1, l2)
                    check_group(group)
                    groupSum.pop( (n1, n2, l1, l2) )
    keep_group_set2(group)
    print ("[Debug] group size", len(group))
    
    #print (group.keys())
    #print_group_set(group)
    check_group(group)
    label = True
    if len(group) == groupLen:
        label = False
    return label 

def merge_three_groups(group, phasedNodes, graph, groupSum, conflictSet):
    
    temp = sorted(group.items(), key=lambda item:len(item[1]))
    groupKeys = [ele[0] for ele in temp]
    groupKeysSet = set(groupKeys) 
    groupLen = len(group)

    #print ("merge group loose")
    onlyChoice = {}
    for i in range(groupLen-1):
        if groupKeys[i].endswith('_0'):
            continue 
        n1 = groupKeys[i]
        if n1 not in group:
            continue    
        l = myGraph.get_group_neighbors(group, phasedNodes, graph, n1)
        for n2 in l:
            symN1 = tools.get_symmetrical_node( n1 )
            symN2 = tools.get_symmetrical_node( n2 ) 
            w1, w2, w3, w4, l1, l2, e1, e2, e3, e4 = get_groups_sum(n1, n2, symN1, symN2, group, graph, groupSum)
            if w1+w2+w3+w4 == 0:
                continue
            if w1 + w4 - w2 - w3 < 0 and w2 > 0 and w3 >0 and w1+w4==0:
                w1, w2, w3, w4 = w2, w1, w4, w3
                n2, symN2 = symN2, n2

            if w1 + w4 - w2 - w3 > 0 and w1 > 0 and w4 >0 and w2+w3==0:
                neig1 = set( myGraph.get_group_neighbors(group, phasedNodes, graph, n1) )
                neig2 = set( myGraph.get_group_neighbors(group, phasedNodes, graph, n2) )
                if len(neig1) == 0 or len(neig2) == 0:
                    continue
                L1 = neig1.intersection(neig2)

                for n3 in L1:
                    symN3 = tools.get_symmetrical_node( n3 )
                    #print ("[Debug] three groups", n1, n2,  w1, w2, w3, w4, l1, l2)
                    
                    N13_w1, N13_w2, N13_w3, N13_w4, N13_l1, N13_l2 = get_groups_sum(n1, n3, symN1, symN3, group, graph, groupSum)
                    #print ("[Debug] three groups", n1, n3, N13_w1, N13_w2, N13_w3, N13_w4, N13_l1, N13_l2)
                    
                    N23_w1, N23_w2, N23_w3, N23_w4, N23_l1, N23_l2 = get_groups_sum(n2, n3, symN2, symN3, group, graph, groupSum)
                    #print ("[Debug] three groups", n2, n3, N23_w1, N23_w2, N23_w3, N23_w4, N23_l1, N23_l2)

                    if (N13_w2>0 and N13_w3>0) and (N23_w2 >0 and N23_w3>0) and N13_w1 + N13_w4==0 and N23_w1 + N23_w4 == 0:
                        N13_w1, N13_w2, N13_w3, N13_w4 = N13_w2, N13_w1, N13_w4, N13_w3
                        N23_w1, N23_w2, N23_w3, N23_w4 = N23_w2, N23_w1, N23_w4, N23_w3
                        n3, symN3 = symN3, n3

                    if (N13_w1>0 and N13_w4>0) and (N23_w1>0 and N23_w4>0) and N13_w2 + N13_w3==0 and N23_w2 + N23_w3 == 0:
                        merge_two_groups(n1, n2, symN1, symN2, group, conflictSet)
                        print ("[Debug] merge three groups", n1, n2,  w1, w2, w3, w4, l1, l2)
                        #groupSum.pop( (n1, n2, l1, l2) ) #may have swap n2
                        if l2 > l1:
                            merge_two_groups(n2, n3, symN2, symN3, group, conflictSet)
                            print ("[Debug] merge three groups", n2, n3,  N23_w1, N23_w2, N23_w3, N23_w4, N23_l1, N23_l2)
                            #groupSum.pop( (n2, n3, N23_l1, N23_l2) )
                        else:
                            merge_two_groups(n1, n3, symN1, symN3, group, conflictSet)
                            print ("[Debug] merge three groups", n1, n3, N13_w1, N13_w2, N13_w3, N13_w4, N13_l1, N13_l2)
                            #groupSum.pop( (n1, n3, N13_l1, N23_l2) )

    keep_group_set2(group, phasedNodes)
    label = True
    if len(group) == groupLen:
        label = False
    return label 


def merge_two_groups_once(graph, group, phasedNodes, threshold, groupSum, conflictSet):
    
    temp = sorted(group.items(), key=lambda item:len(item[1])) #reverse=True
    groupKeys = [ele[0] for ele in temp]
    groupLen = len(group)
    #index = tools.binarySearch_ListValue(temp, 1) # first group size >1
    #print ("before merge group: %s pairs of group" %  (groupLen/2) )
 
    for i in range(groupLen-1):   
        if groupKeys[i].endswith('_0'):
            continue
        bestRate = threshold
        bestN1, bestN2 = "", ""
        bestSymN1, bestSymN2 = "", ""
        n1 = groupKeys[i]
        if n1 not in group:
            continue    
        l = myGraph.get_group_neighbors(group, phasedNodes, graph, n1)
        for n2 in l:
            symN1 = tools.get_symmetrical_node( n1 )
            symN2 = tools.get_symmetrical_node( n2 ) 
            w1, w2, w3, w4, l1, l2, e1, e2, e3, e4= get_groups_sum(n1, n2, symN1, symN2, group, graph, groupSum)
            if l1 >1000 and l2 >1000:
                continue
            s1 = (w1 + w4 - w2 - w3)
            s2 = (w2 + w3 - w1 - w4)
            if (w1+w4)>0:
                r1 = (w1+w4)/(len(e1)+len(e4))
            if (w2+w3)>0:    
                r2 = (w2+w3)/(len(e2)+len(e3))
            if (s1  > 2 and r1>1.0 and w1>0 and w4>0 
                 and (w2+w3)<min(100, 0.3*(w1+w4)) and min(w1,w4) >max(w2,w3) ):
                
                if l1>1 and l2>1 and (w1<=2 or w4<=2):
                    continue
                rate = float(s1) / (l1+l2)
                if rate > bestRate:
                    bestRate = rate
                    bestN1, bestN2 = n1, n2
                    bestSymN1, bestSymN2 = symN1, symN2
            elif s2  > 2 and r2>1.0 and min(w2, w3)>max(w1, w4):
                if w2+w3 <= 3*(w1+w4) or w2==0 or w3==0 or (w1+w4)>100:
                    continue
                if l1>1 and l2>1 and (w1<=2 or w2<=2):
                    continue
                rate = float(s2) / (l1 + l2)
                if rate > bestRate:
                    bestRate = rate
                    bestN1, bestN2 = n1, symN2
                    bestSymN1, bestSymN2 = symN1, n2

        if bestN1 != "":               
            l1, l2= len(group[bestN1]), len(group[bestN2])
            print ("[Debug] merge group once", bestRate)
            if (bestN1, bestN2, l1, l2) in groupSum:
                print (bestN1, bestN2, l1, l2, groupSum[(bestN1, bestN2, l1, l2) ][:4])
                groupSum.pop( (bestN1, bestN2, l1, l2) )
            elif (bestN1, bestSymN2, l1, l2) in groupSum:
                print (bestN1, bestSymN2, l1, l2, groupSum[(bestN1, bestSymN2, l1, l2) ][:4])
                groupSum.pop( (bestN1, bestSymN2, l1, l2) )
            else: 
                print ("debug of groupSum", bestN1, bestN2, l1, l2)
                print (groupSum)
              
            merge_two_groups(bestN1, bestN2, bestSymN1, bestSymN2, group, conflictSet)
            keep_group_set2(group, phasedNodes)
    label = True
    if len(group) == groupLen:
        label = False
    return label 

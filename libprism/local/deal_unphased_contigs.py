#########################################################################
# File Name: deal_unphased_contigs.py
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Wed 22 Sep 2021 04:26:47 PM AEST
#########################################################################
#!/bin/bash
import os
import sys
import tools

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def read_blasr(filename):
    ######
    # if contigs align two sym refs, keep both alignment
    # otherwise, only keep best alignment
    #####
    contigMap = {}
    refMap = {}
    dealedC =set()
    partAlign = set() # contigID set
    #removeContig =set()
    with open(filename, "r") as f:
        for line in f:
            words = line.strip().split(' ')
            contigID=words[0]
            if contigID in dealedC:
                continue
            #print words[9], words[10], words[11]
            cLen = int(words[1])
            if cLen < 5000:
                #print (contigID, cLen, "too short, donot merge")
                continue
            cStart = int (words[2])
            cEnd = int (words[3])
            refID  = words[6]
            rLen = int(words[7])
            rStart = int (words[8])
            rEnd = int (words[9])
            numMatch = int(words[12]) 
            
            if numMatch < 2000 or numMatch <0.1*cLen:
                continue

            if rStart> 0.01*rLen and rEnd< 0.99*rLen and cStart> 0.01*cLen and cEnd< 0.99*cLen:
                partAlign.add(contigID)
                continue

            score = int(words[11]) 
            if contigID not in contigMap:
                contigMap[ contigID ] = []
                contigMap[ contigID ].append(  words )
                if refID not in refMap:
                    refMap[ refID ] = []
                refMap[ refID ].append(words)    
            elif len(contigMap[ contigID] ) == 1:
                rName0=contigMap[ contigID ][0][6]
                if rName0 == refID:
                    print ("best two alignment align to same refID")
                    score0 = int (contigMap[ contigID ][0][11])
                    print (contigID, "two scores:", score, score0)
                    assert score >= score0 
                    continue
                dealedC.add(contigID)
                rName1 = refID
                if rName0.split('_')[1] != rName1.split('_')[1]:
                    # best two not symm, only keep one 
                    print ("unphased contig align to two different phased contig")
                    print (rName0, rName1)
                    #print (contigMap[ contigID ][0][0:12] )
                    #print (words[0:12] )
                    #removeContig.add(contigID)
                else:
                    contigMap[ contigID ].append(  words )
                    if refID not in refMap:
                        refMap[ refID ] = []
                    refMap[ refID ].append(words)    
            else:    
                if score <  int(contigMap[ contigID][1][11]):
                    print (contigID, "error: the first two not with highest score")
                    #sys.exit()
                assert score >= int(contigMap[ contigID][1][11]) 
    '''  
    for c in contigMap:
        print c, len(contigMap[c])
        for words in contigMap[c]:
            print words[:9]
    sys.exit()
    '''
    return contigMap, refMap, partAlign

def read_contigs(filename):
    contigs = {} 
    length = 0
    shortC = set()
    for seq_re in SeqIO.parse( filename, "fasta"):
        contigs[ seq_re.id ] = seq_re
        l = len(seq_re.seq)
        if l <5000:
            shortC.add(seq_re.id)
        length += l
    return contigs, shortC, length


def get_extend_contig_info2(words):
    cName = words[0]
    cLen = int(words[1])
    cStart = int (words[2])
    cEnd = int (words[3])
    cD = words[4] # direction

    rName = words[6]
    rLen = int (words[7])
    rStart = int (words[8])
    rEnd = int (words[9])
    rD = words[10] # direction
    numMatch, numMismatch, numIns, numDel = (
            int(words[12]), int(words[13]), int(words[14]), int(words[15]) )
    score = int(words[11])  
    print (cName, rName, cD, rD, cLen, cStart, cEnd, rLen, rStart, rEnd, numMatch, numMismatch, numIns, numDel)
    
    mergeList = []
    if cD == '+' and rD=='+':
        if cLen - cEnd < max(100, 0.1*cLen) and rStart<100 and numMatch>5000:
            mergeList.append((score, (cName, 0, cStart, rName, rStart, rLen) ) )
        elif rLen -rEnd < 100 and cStart< max(100, 0.1*cLen) and numMatch>5000:
            mergeList.append((score,  (rName, 0, rEnd, cName, cEnd, cLen )) )    
        else:
            #unknowContig.add(cName)
            # wait to do; deal with overhang
            flag = False
            if numMatch > 10000:
                if cStart > rStart:
                    overhang = rStart
                    #assert cLen - cEnd < rLen - rEnd
                    if cLen - cEnd < rLen - rEnd:
                        overhang += cLen - cEnd
                        print (cLen - cEnd)
                        print ("overhang", overhang)
                        if numMatch > overhang:
                            mergeList.append((score,  (cName, 0, cStart, rName, rStart, rLen)) )   
                            flag = True
                else:
                    overhang = cStart
                    #assert cLen - cEnd > rLen - rEnd 
                    if cLen - cEnd > rLen - rEnd:
                        overhang += rLen - rEnd
                        print ("overhang", overhang)
                        if numMatch > overhang:
                            mergeList.append((score,  (rName, 0, rEnd, cName, cEnd, cLen )) )     
                            flag = True
            if flag != True:                
                print ("++ unknow")
    elif cD == '+' and rD=='-':
        if cStart < max(100, 0.1*cLen) and rStart<100 and numMatch>5000:
            mergeList.append( (score, ("reverse", cName, 0, cLen-cEnd, rName, rStart, rLen) )) 
        elif cLen-cEnd < max(100, 0.1*cLen) and rLen - rEnd < 100 and numMatch>5000:
            mergeList.append( (score, (rName, 0, rEnd, cName, cLen-cStart, cLen ,"reverse")) )
        else:
            #unknowContig.add(cName)
            flag = False
            if numMatch > 10000:
                if cStart > rLen-rEnd:
                    overhang = rLen-rEnd
                    #assert rStart > cLen - cEnd
                    if rStart > cLen - cEnd:
                        overhang += cLen - cEnd
                        print ("overhang", overhang)
                        if numMatch > overhang:
                            mergeList.append( (score, (rName, 0, rEnd, cName, cLen-cStart, cLen ,"reverse")) )
                            flag = True
                else:
                    overhang = cStart
                    #assert rStart < cLen - cEnd
                    if rStart < cLen - cEnd:
                        overhang += rStart
                        print ("overhang", overhang)
                        if numMatch > overhang:
                            mergeList.append( (score, ("reverse", cName, 0, cLen-cEnd, rName, rStart, rLen) )) 
                            flag = True
            if flag != True:                
                print ("+- unknow")
            '''    
            print ("+- unknow")
            if numMatch > 2000 and ((cLen-cEnd<0.1*cLen and rStart<0.1*rLen) 
                    or (rLen-rEnd<0.1*rLen and cStart<0.1*cLen)):
                print ("double check")
            '''     
    else:
        #unknowContig.add(cName)
        print ("-- unknow")
    #print "\n"
    return mergeList



def get_extend_contig_info(words):
    cName = words[0]
    cLen = int(words[1])
    cStart = int (words[2])
    cEnd = int (words[3])
    cD = words[4] # direction

    rName = words[6]
    rLen = int (words[7])
    rStart = int (words[8])
    rEnd = int (words[9])
    rD = words[10] # direction
    numMatch, numMismatch, numIns, numDel = (
            int(words[12]), int(words[13]), int(words[14]), int(words[15]) )
    score = int(words[11])
    
    print (cName, rName, cD, rD, cLen, cStart, cEnd, rLen, rStart, rEnd, numMatch, numMismatch, numIns, numDel)
    
    mergeList = []
    if cD == '+' and rD=='+':
        if cLen - cEnd < max(100, 0.1*cLen) and rStart<100 and numMatch>5000:
            mergeList.append((score, (cName, 0, cStart, rName, rStart, rLen) ) )
        elif rLen -rEnd < 100 and cStart< max(100, 0.1*cLen) and numMatch>5000:
            mergeList.append((score,  (rName, 0, rEnd, cName, cEnd, cLen )) )    
        else:
            #unknowContig.add(cName)
            print ("++ unknow")
            if numMatch > 2000 and ((cLen-cEnd<0.1*cLen and rStart<0.1*rLen) 
                    or (rLen-rEnd<0.1*rLen and cStart<0.1*cLen)):
                print ("double check")

    elif cD == '+' and rD=='-':
        if cStart < max(100, 0.1*cLen) and rStart<100 and numMatch>5000:
            #print ("two equal")
            #print ("reverse", rName, 0, rLen-rStart, cName, cEnd, cLen)
            #mergeList.append( ("reverse", rName, 0, rLen-rStart, cName, cEnd, cLen ))
            #print ("reverse", cName, 0, cLen-cEnd, rName, 0, rLen)
            mergeList.append( (score, ("reverse", cName, 0, cLen-cEnd, rName, 0, rLen) )) 
        elif cLen-cEnd < max(100, 0.1*cLen) and rLen - rEnd < 100 and numMatch>5000:
            #print ("two equal")
            #print (cName, 0, cStart, rName, rLen-rEnd, rLen, "reverse")
            #mergeList.append( (cName, 0, cStart, rName, rLen-rEnd, rLen, "reverse"))
            #print (rName, 0, rLen, cName, cLen-cStart, cLen ,"reverse")
            mergeList.append( (score, (rName, 0, rLen, cName, cLen-cStart, cLen ,"reverse")) )
        else:
            #unknowContig.add(cName)
            print ("+- unknow")
            if numMatch > 2000 and ((cLen-cEnd<0.1*cLen and rStart<0.1*rLen) 
                    or (rLen-rEnd<0.1*rLen and cStart<0.1*cLen)):
                print ("double check")
    else:
        #unknowContig.add(cName)
        print ("-- unknow")
    #print "\n"
    return mergeList

def write_final_contigs(finalC, usedC, contigs, unMapContigs, unKnowContigs):
    Acontigs=[]
    Bcontigs=[]
    for re in finalC:
        if re.id.count('groupA_')>=1:
            Acontigs.append(re)
        elif re.id.count('groupB_')>=1:
            Bcontigs.append(re)
        else:
            print ("error")
            sys.exit()

    for re in contigs:
        if re in usedC:
            continue
        if contigs[re].id.count('groupA_')>=1:
            Acontigs.append( contigs[re] )
        elif contigs[re].id.count('groupB_')>=1:
            Bcontigs.append( contigs[re] )

    print ("final group size", len(Acontigs), len(Bcontigs))
    SeqIO.write(Acontigs, 'final_groupA.fasta', 'fasta')
    SeqIO.write(Bcontigs, 'final_groupB.fasta', 'fasta')

    Ccontigs=[]

    for re in unMapContigs:
        Ccontigs.append( contigs[re] )
        Acontigs.append( contigs[re] )
        Bcontigs.append( contigs[re] )
    #print

    print ("final group (u)  size", len(Acontigs), len(Bcontigs))
    SeqIO.write(Acontigs, 'final_groupA_u.fasta', 'fasta')
    SeqIO.write(Bcontigs, 'final_groupB_u.fasta', 'fasta')

    for re in unKnowContigs:
        Ccontigs.append( contigs[re] )
        Acontigs.append( contigs[re] )
        Bcontigs.append( contigs[re] )
    #print

    print ("final group (u_u)  size", len(Acontigs), len(Bcontigs))
    SeqIO.write(Acontigs, 'final_groupA_u_u.fasta', 'fasta')
    SeqIO.write(Bcontigs, 'final_groupB_u_u.fasta', 'fasta')

    SeqIO.write(Ccontigs, 'unknow.fasta', 'fasta')
    return

def check_direction(r, eles):
    #inter =
    print ("debug", r, eles)
    newEles = []
    e1= eles[0]
    e2= eles[1]
    if e2[0] != r:
        e1, e2 =e2, e1
    if e2[0] == r:
        #print (e1[-3:])
        #print (e2[0:3])
        if e1[-3] == r:
            start = e1[:-3]
            end = e2[3:]
            newEles.append(start)
            mid = (r, e1[-2], e2[2])
            newEles.append(mid)
            newEles.append(end)
            print (newEles)
            return newEles
    print ("direction problem")
    newEles = eles[0]
    # can be improved
    #sys.exit()
    return newEles




def merge_contig(mergeList, un_contigs, contigs):

    print (len(mergeList), " pairs need to be merged")
    contigs.update(un_contigs)
    finalC = []
    AusedC, BusedC = set(), set()
    for (r, eles) in mergeList:
        #print (eles)
        l = len(eles) 
        if l == 1:
            merge_two_contigs(r, eles[0], AusedC, BusedC, finalC)
        elif l == 2:
            eles = check_direction(r, eles)
            if len(eles) == 3:
                merge_three_contigs(r, eles, AusedC, BusedC, finalC)
                #print ("3 used size", len(AusedC), len(BusedC))
            else:
                print ("error 2")
                sys.exit()
        else:
            print (r, l)
            print ("error 3")
            sys.exit()
     
    usedC = AusedC.union(BusedC)
    return finalC, usedC


def merge_three_contigs(r, eles, AusedC, BusedC, finalC):
       
    assert len(eles) == 3
    seq=""
    names=[]

    #print ("used size", len(AusedC), len(BusedC))
    for i in range(3):
        if len(eles[i]) == 3:
            cName, cStart, cEnd = eles[i]
            seq = seq  + contigs[cName].seq[cStart:cEnd]
        elif len(eles[i]) == 4:
            if eles[i][0] == "reverse":
                cName, cStart, cEnd = eles[i][1:]
            elif eles[i][-1] == "reverse":
                cName, cStart, cEnd = eles[i][:-1]
            temp = tools.reverse3(contigs[cName].seq)
            seq = seq + temp[cStart:cEnd]
        names.append(cName) 
        print (cName, len( contigs[cName].seq ) )

    name = '_'.join(names)
    print (names)
    if r.count("groupA_")>=1: 
        if (names[0] in AusedC) or (names[1] in AusedC) or (names[2] in AusedC):
            print (names, "already used")
            #sys.exit()
            return
        #AuesdC = AusedC.union(set(names))
        for n in names:
            AusedC.add(n)
    else:
        if (names[0] in BusedC) or (names[1] in BusedC) or (names[2] in BusedC):
            print (names, "already used")
            #sys.exit() 
            return
        #BusedC = BusedC.union(set(names))
        for n in names:
            BusedC.add(n)
        
    print ("after merge")
    #print ("AusedC", AusedC)
    #print ("used size", len(AusedC), len(BusedC))
    print (name, len(seq), "\n" )
    n = SeqRecord(seq, id=name)
    finalC.append(n)
    return 


def merge_two_contigs(r, eles, AusedC, BusedC, finalC):

    if len(eles) == 6:
        cName, cStart, cEnd, rName, rStart, rEnd = eles
        seq = contigs[cName].seq[cStart:cEnd] + contigs[rName].seq[rStart:rEnd]
    elif len(eles) == 7:
        if eles[0] == "reverse":
            cName, cStart, cEnd, rName, rStart, rEnd = eles[1:]
            temp = tools.reverse3(contigs[cName].seq)
            seq = temp[cStart:cEnd] + contigs[rName].seq[rStart:rEnd]
        elif eles[-1] == "reverse":
            cName, cStart, cEnd, rName, rStart, rEnd = eles[:-1]
            temp = tools.reverse3(contigs[rName].seq)
            seq = contigs[cName].seq[cStart:cEnd] + temp[rStart:rEnd]
    a = cName.count('groupA_')       
    b = cName.count('groupB_')
    if a >=1 or b>=1:
        newId = cName + '_' + rName
    else:
        newId = rName + '_' + cName

    if r.count('groupA_') >= 1:
        if (rName in AusedC) or (cName in AusedC) :
            print (rName, "or", cName, "already used")
            #sys.exit()
            return
        AusedC.add(rName)
        AusedC.add(cName)
    else:  
        if (rName in BusedC) or (cName in BusedC):
            print (rName, "or", cName, "already used")
            #sys.exit()
            return
        BusedC.add(rName)
        BusedC.add(cName)



    print ( cName, cStart, cEnd, len(contigs[cName].seq), 
            rName, rStart, rEnd, len( contigs[rName].seq ) )
    print ( newId, len(seq), "\n" )
    n = SeqRecord(seq, id=newId)
    finalC.append(n)
    return 


def deal_long_unphased_contigs(refMap):
    ######### 
    # remove contigs covered by phased refs
    ##########
    inSize, inNum = 0, 0
    refSet = set() # covered unphased contigID
    for r in refMap:
        l = len( refMap[r] )
        #print (r, l)
        remove=[] # only remove phased contig be covered situation
        # can be improve by extend covered phased contig
        flag = False 
        for i in range(l):
            words = refMap[r][i]
            cID = words[0] 
            cLen = int(words[1])
            cStart = int (words[2])
            cEnd = int (words[3])
            cD = words[4] # direction

            rName = words[6]
            rLen = int (words[7])
            rStart = int (words[8])
            rEnd = int (words[9])
            rD = words[10] # direction
            
            if (rStart <= 0.1*rLen and rEnd >= rLen*0.9):
                #print ("ref be covered") print (words[:11])
                flag = True
                remove.append(  words )
            elif (cD==rD and cStart>=rStart-100 and cLen-cEnd>=rLen-rEnd):
                flag = True
                remove.append(  words )
                #if (cEnd-cStart<cLen*0.5):
                #    print ("++ ref be covered")
                #    print (words[:11])
            elif (cD!=rD and cStart>=rLen-rEnd and rStart<=cLen-cEnd): 
                flag = True
                remove.append(  words )
                #if (cEnd-cStart<cLen*0.5):
                #    print ("+- ref be covered")
                #    print (words[:11])
        if flag == True:
            inSize += rLen
            inNum += 1
            refSet.add(r)
            
        for words in remove:
            refMap[r].remove(words)
    print (inNum, "phased contigs are covered by unphased contigs", len(refSet))
    print (inSize, "length phased contigs are covered by unphased contigs")
    return refSet



def deal_covered_unphased_contigs(refMap):
    ######### 
    # remove contigs covered by phased refs
    ##########
    inSize, inNum = 0, 0
    dealed = set() # dealed covered unphased contigID
    for r in refMap:
        l = len( refMap[r] )
        #print (r, l)
        remove=[]
        for i in range(l):
            words = refMap[r][i]
            cID = words[0] 
            cLen = int(words[1])
            cStart = int (words[2])
            cEnd = int (words[3])
            cD = words[4] # direction

            rName = words[6]
            rLen = int (words[7])
            rStart = int (words[8])
            rEnd = int (words[9])
            rD = words[10] # direction
            numMatch = int(words[12])
            #if cStart <= 100 and cEnd >= cLen-100: #ignore extend length smaller than 200
            if (cStart <= 0.1*cLen and cEnd >= cLen*0.9):
               #or (cEnd-cStart>cLen*0.85 and (rStart>0.1*cLen or rLen-rEnd>0.1*cLen))): #ignore extend length smaller than 200
                #print words[:15]
                #print ("covered")
                #print (words[:11])
                if cID not in dealed:
                    inSize += cLen
                    inNum += 1
                    #print (cID, "alreay covered by other phased contigs")
                    dealed.add(cID)
                remove.append(  words )
            elif (cD==rD and rStart>=cStart and rLen-rEnd>=cLen-cEnd): 
                #if (cEnd-cStart<cLen*0.5):
                #    print ("++ covered")
                #    print (words[:11])

                overhang = cStart
                overhang += cLen - cEnd
                if overhang > numMatch:
                    continue
                if cID not in dealed:
                    inSize += cLen
                    inNum += 1
                    #print (cID, "alreay covered by other phased contigs")
                    dealed.add(cID)
                remove.append(  words )
            elif (cD!=rD and rStart>=cLen-cEnd and cStart<=rLen-rEnd): 
                #if (cEnd-cStart<cLen*0.5):
                #    print ("+- covered")
                #    print (words[:11])
                overhang = cStart
                overhang += cLen - cEnd
                if overhang > numMatch:
                    continue
                if cID not in dealed:
                    inSize += cLen
                    inNum += 1
                    #print (cID, "alreay covered by other phased contigs")
                    dealed.add(cID)
                remove.append(  words )
        for words in remove:
            refMap[r].remove(words)

    print (inNum, "contigs are covered by phased contigs")
    print (inSize, "length unphased contigs are covered by phased contigs")
    return dealed

def choose_contig(r, refMap):
    mt,res = [],[]
    for words in refMap:            
        #mergeTemp = get_extend_contig_info(words[:16])
        mergeTemp = get_extend_contig_info2(words[:16])
        mt.extend(mergeTemp)
        #print ("1")
        #print (words[:16])

    #print ("2")
    #print (mt)

    if len(mt) <=1:
        for (a, b) in mt:
            res.append(b)
    else:
        mt_pre, s_pre="", 0
        mt_sub, s_sub="", 0
        for (a, b) in mt:
            assert b[0]==r or b[-3]==r
            if b[0] == r and a<s_pre:
                mt_pre = b
            elif b[-3] == r and a<s_sub:
                mt_sub = b
        if mt_sub != "":
            res.append(mt_sub)
        if mt_pre != "":    
            res.append(mt_pre)

    #print ("3")
    #print (res)
    return res

if __name__ == "__main__":

    #input: blasr, unphased_contig.fasta, phased_contigs.fasta 
    filename = sys.argv[1]
    contigMap, refMap, partAlign = read_blasr(filename)
    un_contigs, un_s, l1 = read_contigs(sys.argv[2]) # un_s: unphased short contigs
    print ('unphased contig size:', len(un_contigs) )
    print ('unphased contig length:', l1)

    contigs, s, l2 = read_contigs(sys.argv[3]) 
    print ('phased contig size:', len(contigs) )
    print ('phased contig length:', l2)

    print ("part align contigs number", len(partAlign)) 

    unMapContig = set(un_contigs.keys()) - un_s - set(contigMap.keys()) - partAlign 
    print ("unMap contigs number", len(unMapContig) )
    #print (sorted(unMapContig))
    unMapSize = 0
    for c in unMapContig:
        unMapSize += len(un_contigs[c].seq)

    print (unMapSize, "length contigs still unmaped")

    #unknowContig.update(unMapContig)
    coveredSet = deal_covered_unphased_contigs(refMap) #unphased contigs are covered by phased contig
    longSet = deal_long_unphased_contigs(refMap) #phased contigs are covered by unphased contig
    #unknowContig = set()
    mergeList = []
    for r in refMap:
        l = len( refMap[r] )
        #print (r, l)
        mt = choose_contig(r, refMap[r])
        l2 = len(mt) 
        assert l2 <= 2
        if l2 >0:
            mergeList.append( (r, mt))
    finalC, usedSet = merge_contig(mergeList, un_contigs, contigs) # during merge, all used contigs

    unknowContig = set(un_contigs.keys()) - un_s - coveredSet - usedSet -unMapContig 
    print ("unKnow contigs number", len(unknowContig) )
    print (sorted(unknowContig))

    unknowSize = 0
    for c in unknowContig:
        unknowSize += len(un_contigs[c].seq)

    print (unknowSize, "length contigs still unknow (unphased)")

    write_final_contigs(finalC, usedSet, contigs, unMapContig, unknowContig)    

'''    
if __name__ == "__main__":
    filename = sys.argv[1]
    contigMap, refMap = read_blasr(filename)
    un_contigs = read_contigs(sys.argv[2])
    contigs = read_contigs(sys.argv[3])
    
    unMapContig = set(un_contigs.keys()) - set(contigMap.keys())  
    print "unMap contigs number", len(unMapContig)
    print unMapContig
    unknowContig = set()
    unknowContig.update(unMapContig)
    inSize, inNum = 0, 0
    unknowSize = 0
    mergeList = []
    for c in contigMap:
        words = contigMap[c][0]
        cLen = int(words[1])
        cStart = int (words[2])
        cEnd = int (words[3])
        if cStart <= 100 and cEnd >= cLen-100: #ignore extend length smaller than 200
            #print words[:15]
            inSize += cLen
            inNum += 1
            print c, "alreay in phased contigs"
        else:
            print words[:15]
            unknowSize, mergeTemp = get_extend_contig_info(words[:16], unknowSize, unknowContig)
            mergeList.extend(mergeTemp)
            if len(contigMap[c]) >1:
                words = contigMap[c][1]
                print words[:15]
                unknowSize, mergeTemp = get_extend_contig_info(words[:16], unknowSize, unknowContig)
                mergeList.extend(mergeTemp)
    print inNum, "contigs are covered by phased contigs"        
    print inSize, "length unphased contigs are covered by phased contigs"      
    print "unkonw contigs number", len(unknowContig)
    print unknowContig
    print unknowSize, "length contigs still unknow (unphased)"       
  
    merge_contig(mergeList, un_contigs, contigs)
'''

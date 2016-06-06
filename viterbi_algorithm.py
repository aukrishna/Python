#Assignment 4
#Krishna Vinnakota
#UW NetID: VINNAK

##To Print:
##(a) the HMM emission/transition parameters,
##(b) the log probability (natural log, base-e) of the overall Viterbi path,
##(c) the total number of "hits," found, and
##(d) the lengths and locations (starting and ending positions) of the first k
##"hits," i.e., the first k (contiguous) subsequences that are assigned to state 2
##in the Viterbi path. (Print all hits, if there are fewer than k of them. Recall
##that by convention genomic positions are 1-based, not 0-based, indices.)

##Let's also have a language bake-off: also print
##(a) the language in which you coded
##(b) the total CPU time taken by your algorithm for the first 9 Viterbi iterations.
##--Do NOT include the time taken to read the genome sequence, nor the time taken by the 10th iteration, which may be dominated by printing hits
##Include basic info about the computer on which you ran the timing test ("6 bazigahertz Intelerola winter olympathalon with 33 bytes of RAM").
##(c) please feel free to post this information, and the first 10 hits of length >= 50 from your final iteration, to the catalyst GoPost site, so you can (partially) check your work against each other
##--and see whether language matters.
    
import sys
import datetime
import math


#HMM Parameters
intialTransitionProbabilities=[0.9999,0.0001]
transitionProbabilities=[[0.9999,0.0001],[0.01,0.99]]
EmissionProbabilities= [{'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}, {'A': 0.20, 'C': 0.30, 'G': 0.30, 'T': 0.20}]
numberOfHiddenStates = 2

filePath='NC_000909.fna.txt'
logProbability=0.00
emittedSequence = ""

iterations=10


##Functions

#This function takes the filepath as parameter and
#returns the genome sequence in the file as a string and replaces not ACGT charecters with T
def getSequence(filePath):
    gSequence=""
    gList=[]
    #tlist =[]    
    with open(filePath) as f:
        content = f.readlines()
    for i in range(1, len(content)):
        gSequence=gSequence+content[i].rstrip()
    #gList=[]*len(gSequence)
    #print(len(gSequence))
    gSequence=gSequence.strip()
    gSequence=gSequence.rstrip()
    gList = list(gSequence)
    #print(len(gSequence),len(gList))
    for i in range(0,len(gList)):
        if gList[i] not in ['A','C','G','T']:
            #tlist.append(gList[i])
            gList[i]='T'                    
    #print(len(tlist))
    #print(tlist)       
    gSequence = ''.join(gList)        
    return gSequence

#This function prints the language information, total CPU time taken in seconds, basic computer information.
def languageBakeOff():
    print("Language Bake Off Information:")
    print("------------------------------------------------------------------------------")
    print("Language in which I coded: Python: Python 3.5.1")
    #print("CPUTimeTaken:",CPUTimeTaken.total_seconds(),"Seconds")
    print("Computer Configuration: Intel Core i7-4600U CPU @ 2.10 GHz with 8 GB of RAM")
    print("------------------------------------------------------------------------------")

def CPUTimeTaken(CPUTimeTaken):
    print("CPUTimeTaken:",CPUTimeTaken.total_seconds(),"Seconds")

#function to print HMM Parameters
def printHMMParameters():
    print("TransistionProbabilities:")
    for i in range(0,len(transitionProbabilities)):
        for j in range(0,len(transitionProbabilities[i])):
            print(transitionProbabilities[i][j])
    print("EmissionProbabilities:")
    for i in range(0,len(EmissionProbabilities)):
        #print('\n')
        print("State:",i)
        for key,value in EmissionProbabilities[i].items():
            print(key,value)

#function to update HMM Parameters
def UpdateHMMParameters(hiddenStates, seq):
    transitions=[[0.00,0.00],[0.00,0.00]]
    totalTransitionCount=[0,0]
    totalEmissionCount=[0,0]
    emissions=[{'A':0.00,'C':0.00,'G':0.00,'T':0.00},{'A':0.00,'C':0.00,'G':0.00,'T':0.00}]
    emissions[hiddenStates[0]][seq[0]]=1
    for i in range(1,len(seq)):
        transitions[hiddenStates[i]][hiddenStates[i-1]] = transitions[hiddenStates[i]][hiddenStates[i-1]]+1
        totalTransitionCount[hiddenStates[i-1]]=totalTransitionCount[hiddenStates[i-1]]+1
        if seq[i] not in emissions[hiddenStates[i]]:
            emissions[hiddenStates[i]][seq[i]]=0
        emissions[hiddenStates[i]][seq[i]]=emissions[hiddenStates[i]][seq[i]]+1
        totalEmissionCount[hiddenStates[i]]=totalEmissionCount[hiddenStates[i]]+1
    for i in range(0,2):
        for key, value in emissions[i].items():
            #EmissionProbabilities[value][i]=value/totalEmissionCount[i]
            EmissionProbabilities[i][key]=value/totalEmissionCount[i]

    for i in range(0,2):
        for j in range(0,2):
            transitionProbabilities[j][i]=transitions[j][i]/totalTransitionCount[i]          

  
#Implementation of viterbi 
def viterbiAlgorithm(sequence,intialTransitionProbabilities,transitionProbabilities,EmissionProbabilities):
    #print("in printHMMParams")
    seqLength = len(sequence)
    finalState = -1
    #create two dimensional array with [seqLength,numberOfHiddenStates] size
    maxStates = [[0 for x in range(seqLength)] for y in range(numberOfHiddenStates)]    
    prevSeqProbabilities = [0.0,0.0]
    currSeqProbabilities = [0.0,0.0]
    for i in range(0,numberOfHiddenStates):
        maxStates[i][0]=i
        prevSeqProbabilities[i]=math.log(intialTransitionProbabilities[i])+math.log(EmissionProbabilities[i][sequence[0]])
    for i in range(1,seqLength):
        currEmission = sequence[i]
        for k in range(0,numberOfHiddenStates):
            maxSeqProbability=-sys.float_info.max 
            maxState = -1
            for j in range(0,numberOfHiddenStates):
                seqProbability = prevSeqProbabilities[j]+math.log(transitionProbabilities[k][j])+math.log(EmissionProbabilities[k][currEmission])
                if seqProbability > maxSeqProbability:
                    maxSeqProbability = seqProbability
                    maxState = j
            currSeqProbabilities[k]=maxSeqProbability
            maxStates[k][i]=maxState
        prevSeqProbabilities=currSeqProbabilities
        currSeqProbabilities=[0.00,0.00]
    maxTotalProbability = -sys.float_info.max
    for k in range(0,numberOfHiddenStates):
        if prevSeqProbabilities[k]>maxTotalProbability:
            maxTotalProbability=prevSeqProbabilities[k]
            finalState = k
            
    #logProbability = maxTotalProbability
    print("Log Probability:",maxTotalProbability)
    return traceBack(finalState, maxTotalProbability, maxStates, sequence)
            
#Implementation of traceback
def traceBack(finalState, maxTotalProbability, maxStates, sequence):
    #print("in Traceback")
    emittedSequence = sequence
    #logProbability = maxTotalProbability
    seqLength = len(emittedSequence)
    currState = finalState
    while (seqLength>0):
        hiddenStates[seqLength-1]=currState
        seqLength=seqLength-1
        #currState = maxStates[seqLength][currState]
        currState = maxStates[currState][seqLength]
    #return maxTotalProbability



#Function to print CPG Islands in the sequence - Updated function
def printCPGIslandInfo1(hiddenStates,hitsToPrint):
    #print("printing")
    curPos=1
    cpgIslandBegin = -1
    cpgIslandEnd = -1
    totalHits = 0;
    for i in range(0, len(hiddenStates)):
        if hiddenStates[i] == 1:
            if cpgIslandBegin > 0:
                cpgIslandEnd = curPos
            else:
                cpgIslandBegin = curPos
                cpgIslandEnd = curPos
        else:
            if cpgIslandBegin > 0:
                totalHits = totalHits+1
                if totalHits <= hitsToPrint:
                    print("Start:",cpgIslandBegin,"End:",cpgIslandEnd,"Length:",cpgIslandEnd-cpgIslandBegin+1)
                cpgIslandBegin = -1
                cpgIslandEnd = -1
        curPos = curPos+1
    if cpgIslandBegin > 0:
        totalHits = totalHits+1
        if totalHits <= hitsToPrint:
            print("Start:",cpgIslandBegin,"End:",cpgIslandEnd,"Length:",cpgIslandEnd-cpgIslandBegin+1)
        cpgIslandBegin = -1
        cpgIslandEnd = -1
    print("Total Hits:",totalHits)

                
                


    
##Main Program    
languageBakeOff()
seq = getSequence(filePath)
#print(seq[0],seq[1],seq[2],seq[3])
hiddenStates=[0]*len(seq)
totalTimeTaken = 0
#startTime = datetime.datetime.now()
for i in range(1,iterations+1):
    print("------------------------------------------------------")
    print("iteration:",i)
    print("------------------------------------------------------")    
    printHMMParameters()
    if i < iterations:
        startTime = datetime.datetime.now()
        viterbiAlgorithm(seq,intialTransitionProbabilities,transitionProbabilities,EmissionProbabilities)
        endTime = datetime.datetime.now()
        totalTimeTaken = totalTimeTaken + (endTime-startTime).total_seconds()
    #print("unique hiddenStates",set(hiddenStates))
        print("First Ten Hits:")
        printCPGIslandInfo1(hiddenStates,10)
        #print("LogProbability:",logProbability)
    else:
        viterbiAlgorithm(seq,intialTransitionProbabilities,transitionProbabilities,EmissionProbabilities)
        print("All Hits:")
        printCPGIslandInfo1(hiddenStates,sys.maxsize)
        #print("LogProbability:",logProbability)
    #printHMMParameters()
    #print("LogProbability:",logProbability)
    UpdateHMMParameters(hiddenStates,seq)
    print("\n")
print("\n")
print("Execution Time for first 9 iterations:",totalTimeTaken,"Seconds")






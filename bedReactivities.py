#!/usr/bin/env python

'''
Script to find the 5'coverage and overall coverage of given bed regions using assigned .coverage files.
Bed regions assumed to describe one transcript or regions. Bed regions should be sorted by position and non-overlapping.
Coverage is returned as unbroken columns describing each nucleotide in the bed regions in the 5'-3' direction.
Overall coverage and 5'coverage are written to separate files.
(- strand bed regions get reversed)
USAGE: python bedCoverage file.bed [SLBPuv | SLBPform | SLBPnai | FOSIR]
'''    
import sys
import math
import numpy as np
import getopt
import os
from normalizeSHAPE import *
import pybedtools as bedtools
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

tnxDict = {}
SCRIPTPATH = __file__.rstrip("bedReactivities.py")
genefile = open(SCRIPTPATH+"gene2name.txt",'r')
genes = genefile.read().splitlines()
genefile.close()
for line in genes:
    thisline = line.split()
    gene = thisline[1]  # name
    transcript = thisline[0]  # geneid
    tnxDict[transcript] = gene

global GENOME, DATAPATH, BASENAME
BASENAME=""

def usage():
    print('''python bedCoverage.py -i input_regions.bed [...]
        Options:
        -h, --help

        REQUIRED
        -i, --input     The input text file of regions in BED format. Requires strand information.
        
        -g, --genome    The genome file (.fa) for the given data. Should match the genome build reads were aligned to.

        -a, --treated   Comma-separated list of paths to directories contaiting coverage files generated for treated (or in vitro) samples.

        -b, --untreated Comma-separated list of paths to directories contaiting coverage files generated for untreated (or in vivo) samples.

        OPTIONAL
        -n, --name    Specify that regions should be output by BED name instead of BED region. Useful for outputting a reactivity profile across a spliced transcript.

        -o, --output   The output files' basename. Outputs include .map files for each region with reactivity per base, and a .bigwig file describing reactivties across all regions.
        ''')

def getAverage(reactivities): #return the average SHAPE reactivity of a region
    #print("getAverage")
    total = 0.
    counter = 0.
    for i in range(0,len(reactivities[0])):
        values = reactivities[reactivities[:,i]!=-999]
        if len(values)!=0: #all the values were -999
            total+=np.mean(values[:,i])
            counter+=1
    if counter==0:
        total = -999
        counter = 1
    average = total/counter
    return average

def recSearch(LHS,RHS,lastLHS,SORTED,TARGET): #SORTED is the sorted array to search for TARGET in
    mid = (RHS+LHS)//2
    if TARGET == SORTED[mid]:
        return mid
    elif LHS >= RHS:
        return lastLHS
    elif TARGET < SORTED[mid]:
        return recSearch(LHS,mid-1,LHS,SORTED,TARGET)
    else:
        return recSearch(mid+1,RHS,mid,SORTED,TARGET)

def binarySearchClosest(SORTED,TARGET):
    if len(SORTED)==0:
        return 0
    else:
        return recSearch(0,len(SORTED)-1,0,SORTED,TARGET)

#the cov files only list coverge for genomic pos with > 0 coverage

def cov_from_files(covFiles,bedPosIn,bedSum,chrom,add1): #returns np matrix with bed region's coverage from each of files in covfiles
    cov = np.zeros((len(covFiles),bedSum), dtype='int')
    for i in range(0,len(covFiles)):
        bedPos = [] #make a copy of the bed regions
        for element in bedPosIn:
            bedPos.append(list(element))
        thiscov = np.zeros(bedSum,dtype='int')+add1 #will contain coverages for each nucleotide position
        counter = 0 #keeps track of how many bp of the combined bed regions we've gone through
        fh = open(covFiles[i],'r')
        fhlist = fh.read().splitlines() #assuming coverage files small enough to read into a list
        fh.close()
        positions = np.zeros(len(fhlist), dtype="int")
        coverage = np.zeros(len(fhlist),dtype='int')
        for line in range(0,len(fhlist)): #cov file format: chrom position coverageVal
            thisline = fhlist[line].split()
            positions[line] = int(thisline[1])
            coverage[line] = int(thisline[2])
          
        #Going through the coverage listed in the file. The coverage files don't list pos with 0 coverage
        #So I have to deal with the case where bed regions (or parts thereof) are not in the cov file.
        #First use binary search to find a starting point in the file:
        startPos = binarySearchClosest(positions,int(bedPos[-1][0])) #assuming bedPos isn't empty now
        #print("START POS in COV FILE:")
        #print(startPos, positions[startPos])
        for line in range(startPos,len(positions)):
            thispos = positions[line]
            if len(bedPos)==0:
                break
            elif thispos>=int(bedPos[-1][0]):
                while len(bedPos)>0 and thispos>=int(bedPos[-1][0]):
                    if thispos<int(bedPos[-1][1]):
                        counter += thispos-int(bedPos[-1][0])
                        thiscov[counter] = coverage[line]
                        bedPos[-1][0] = thispos+1
                        counter+=1
                    else:
                        counter+=int(bedPos[-1][1])-int(bedPos[-1][0])
                        bedPos.pop()
        cov[i] = thiscov
    return np.transpose(cov)

def getCoverage(bedlines,prefixes): #given sorted matrix of bed regions
    chrom = bedlines[0][0] #assuming all bed regions on same chromosome
    if len(bedlines[0])<6:
        print("Bed regions need strand information.")
        sys.exit()
    strand = str(bedlines[0][5]) #assuming all same strand too
    #put bed regions into a stack (list)
    bedPos = []
    bedSum = 0
    for i in range(len(bedlines)-1,-1,-1): #listing bed regions in reverse order, like a stack
        bedPos.append(bedlines[i][1:3])
        bedSum += int(bedlines[i][2]) - int(bedlines[i][1]) #the total number of bp in the bed regions
    #print(prefixes)
    covFiles = [] #coverage files are split up by chromosome, only need to parse the one corresponding to this bed file
    fileDNE = False
    for i in prefixes:
        strand_str = ".pos" if strand == "+" else ".neg"
        covFiles.append(i+chrom+strand_str+".cov")
        if not os.path.isfile(covFiles[-1]):
            fileDNE = True
    if fileDNE: #for the case of a chromosome that has no coverage file. e.g. weird chromosomes
        return np.zeros((bedSum,len(prefixes)),dtype='int')+1,np.zeros((bedSum,len(prefixes)),dtype='int')
    
    cov = cov_from_files(covFiles,bedPos,bedSum,chrom,1) #starting at a base coverage of 1 because later have to divide by overall coverage

    covFiles = [] #now get the 5' coverage
    for i in prefixes:
        strand_str = ".pos" if strand == "+" else ".neg"
        covFiles.append(i+chrom+strand_str+".mut")
        if not os.path.isfile(covFiles[-1]):
            fileDNE = True
    if fileDNE: #for the case of a chromosome that has no coverage file. e.g. weird chromosomes
        return np.zeros((bedSum,4),dtype='int')+1,np.zeros((bedSum,4),dtype='int')
    cov5 = cov_from_files(covFiles,bedPos,bedSum,chrom,0)
    if strand=="-":
        return cov[::-1], cov5[::-1]
    else:
        return cov,cov5
    
#mergeSort of bed files (sort -k 1,1 -k 2,2n)
def mergeSort(array1, array2): #merge 2d arrays sorted based on chromosome then strand then start position
    merged = [[0]]*(len(array1)+len(array2))
    counter1 = 0
    counter2 = 0
    while counter1 < len(array1):
        if counter2 >= len(array2):
            break
        else:
            sameChrom = (array1[counter1][0]==array2[counter2][0])
            sameStrand = (array1[counter1][5]==array2[counter2][5])
            if(array1[counter1][0] < array2[counter2][0]):
                merged[counter1+counter2] = array1[counter1]
                counter1+=1
            elif(sameChrom and array1[counter1][5]<array2[counter2][5]):
                merged[counter1+counter2] = array1[counter1]
                counter1+=1
            elif(sameChrom and sameStrand and int(array1[counter1][1])<int(array2[counter2][1])):
                merged[counter1+counter2] = array1[counter1]
                counter1+=1
            else:
                merged[counter1+counter2] = array2[counter2]
                counter2+=1
    #add remaining array1 or array2 values to merged (one of them will have leftover values)
    while counter1 < len(array1):
        merged[counter1+counter2] = array1[counter1]
        counter1+=1
    while counter2 < len(array2):
        merged[counter1+counter2] = array2[counter2]
        counter2+=1
    return merged

def binarySortBed(array2d): #2d bed array. Each bed interval per row. Chrom-col1, start-col2, strand-col6
    if len(array2d)<=1:
        return array2d
    else:
        sortedLHS = binarySortBed(array2d[0:len(array2d)//2]) #// = floor division in python3
        sortedRHS = binarySortBed(array2d[len(array2d)//2:])
        return mergeSort(sortedLHS,sortedRHS)
    
def mergeBedSameChrom(beds): #beds is a 2d array of sorted bed regions from the same chrom (and strand)
    merged = [[-1,-1,-1,-1,-1,-1]]
    for i in beds:
        if int(i[1]) < int(merged[-1][2]):
            if int(i[2]) > int(merged[-1][2]):
                merged[-1][2] = i[2]
        else:
            icopy = [0]*len(i) #to avoid modifying original array, make deep copy
            for j in range(0,len(i)):
                icopy[j] = i[j]
            merged.append(icopy)
    merged.pop(0)
    return merged

def mergeRelPos(relPos):
    merged = [[-1,-1]]
    for i in relPos:
        if int(i[0]) < int(merged[-1][1]):
            if int(i[1]) > int(merged[-1][1]):
                merged[-1][1] = i[1]
        else:
            icopy = [0]*len(i) #to avoid modifying original array, make deep copy
            for j in range(0,len(i)):
                icopy[j] = i[j]
            merged.append(icopy)
    merged.pop(0)
    return merged

def relPosBedSameChrom(beds):
    counter= 0
    lastStop = [-1]
    i=0
    relPos = np.zeros((len(beds),2), dtype="int")
    while i < len(beds):
        maxLast = max(lastStop)
        if maxLast!=lastStop[-1]: #the current interval starts inside a large interval that encompasses yet another smaller intervals
            if int(beds[i][1])>maxLast:
                counter+=maxLast-lastStop[-1]
            else:
                counter+=int(beds[i][1])-lastStop[-1]
        else:
            lastStop = lastStop[len(lastStop)-1:] #the maxLastStop is the latest interval, just keep that one
            counter+=min(0,int(beds[i][1])-lastStop[-1])
        relPos[i] = [counter,counter+(int(beds[i][2])-int(beds[i][1]))]
        counter+=(int(beds[i][2])-int(beds[i][1]))
        lastStop.append(int(beds[i][2]))
        #print(lastStop)
        i+=1
    return relPos

#******Modify so that a region containing all -999 values isn't written to file (waste of space)
def write_to_file(cov,name,extension,ignore=-999):
    #the "ignore" argument tells write_to_file to not write output if all values in cov==ignore
    for i in range(0,len(cov[0])):
        values = cov[cov[:,i]!=ignore]
        if len(values)==0: #all the values were -999
            return
    print(name+extension)
    outfile = open(name+extension, 'w') #writing given matrices to output files
    for line in cov:
        outline = ""
        for i in line:
            outline += str(i)+"\t"
        outline = outline.rstrip('\t')
        outfile.write(outline+"\n")
    outfile.close()

def average_replicates(rx): #rx = [rep1val,rep2val,etc]
    #ignoring -999 values
    no999 = rx[rx!=-999]
    if len(no999)==0:
        average = -999.
        var = 0.
    else:
        average = np.mean(rx[rx!=-999])
        var = np.var(rx[rx!=-999])
    return average,var

def write_map_file(rx,name,sequence,ignore=-999,basename=BASENAME):
    #the "ignore" argument tells write_to_file to not write output if all values in cov==ignore
    extension = ".map"
    for i in range(0,len(rx[0])):
        values = rx[rx[:,i]!=ignore]
        if len(values)==0: #all the values were -999
            return False
    
    strand = name[-1]
    if strand == '+':
        name = name[:-1] + "-pos"
    elif strand == '-':
        name = name[:-1] + "-neg"
    outfile = open("{}{}{}".format(basename, name.replace(':','-'), extension), 'w') #writing given matrices to output files
    counter = 1
    for line in rx:
        average,var = average_replicates(line)
        outfile.write(str(counter)+"\t"+str(average)+"\t"+str(var)+"\t"+sequence[counter-1]+"\n")
        counter+=1
    outfile.close()
    return True

def write_to_bedgraph(rx_in,chrom,start,strand,ignore=-999,basename=BASENAME):
    extension = ".bedgraph"
    rx = rx_in
    if strand=="-":
        rx = rx_in[::-1]
    for i in range(0,len(rx[0])):
        values = rx[rx[:,i]!=ignore]
        if len(values)==0: #all the values were -999, dont write anything to file
            return
    #since bedgraph is strandless format, will output to separate + and - strand bedgrph files
    formatted_strand = "pos" if strand == '+' else "neg"
    outfile = open(basename+"{}.strand".format(formatted_strand)+extension, 'a') #writing given matrices to output files
    
    counter = int(start)-1
    for line in rx:
        counter+=1
        average,var = average_replicates(line)
        if average==ignore:
            continue
        outfile.write(chrom+"\t"+str(counter)+"\t"+str(counter+1)+"\t"+str(average)+"\n")
    outfile.close()    


def combineCoverage(cov,cov5,relPos_in,beds_in,strand,sequence,USE_BED_NAME = False): #input: coverage 5'coverage arrays and a 2darray of bed regions to extract. relPos format: [[relStart,relStop],...]
        #cov and cov5 profiles are returned reversed by getCoverage() if strand is negative
    if strand=='-': #reversing ordering/ranges of relPos (and beds) if strand is negative
        numRegions = len(relPos_in)
        relPos = np.zeros((numRegions,2), dtype="int")
        length = len(cov) #the number of nucleotides that the regions cover
        for i in range(numRegions-1,-1,-1):
            relPos[numRegions-i-1] = [length-relPos_in[i][1],length-relPos_in[i][0]]
        beds = beds_in[::-1]
    else:
        beds = beds_in
        relPos = relPos_in
    reactivities = []
    if USE_BED_NAME: 
        nameDict = {}
        for i in range(0,len(beds)):
            name = beds[i][3]
            if name in nameDict:
                nameDict[name].append(relPos[i])
            else:
                nameDict[name] = [relPos[i]]
        
        for name in nameDict: #piece together coverages for each seg of regions*
            unmerged_regions = nameDict[name] #Regions that have the same name could overlap, giving redundant cov values
            regions = mergeRelPos(unmerged_regions) #solve this by first merging relPos regions, then combining cov profiles
            regionCov = np.zeros((1,len(cov[0])), dtype='int') #len cov[0] is the number of samples
            regionCov5 = np.zeros((1,len(cov[0])), dtype='int')
            #print(regionCov)
            regionSeq = ""
            for region in regions:
                regionSeq+=sequence[region[0]:region[1]]
                regionCov = np.concatenate((regionCov,cov[region[0]:region[1]]))
                regionCov5 = np.concatenate((regionCov5,cov5[region[0]:region[1]]))
                #print(region[0], region[1])
            #print(regionSeq)
            if name in tnxDict: #to get gene name in output file name
                name += "."+tnxDict[name]
            reactivities = normalizeSHAPEmain(regionCov[1:],regionCov5[1:]+0.,regionSeq,trimEnds=True)
            #write_to_file(reactivities,name,".rx",-999,BASENAME)
            enoughdata = write_map_file(reactivities,name,regionSeq,-999,BASENAME)
            write_to_bedgraph(rx_in=reactivities,chrom='chr1',start=1,strand='+',ignore=-999,basename=BASENAME)
            ## uncomment to output coverages of	individual transcripts output in addition to map and rx	files
            #if enoughdata:
            #    write_to_file(regionCov[1:],name,".cov",1) #last argurment tells write_to_file to not write if all values==1, for ex
            #    write_to_file(regionCov5[1:],name,".5cov",0) 
    else: #simply write the cov of each bed region to individual file
        #print(beds)
        c = 0
        for region in relPos: #relPos and beds are in the same order right?
            #regionSeq = bedSeqDict[beds[c][0]+":"+beds[c][1]+"-"+beds[c][2]+"("+strand+")"]
            regionSeq = sequence[region[0]:region[1]]
            #print(regionSeq)
            reactivities = normalizeSHAPEmain(cov[region[0]:region[1]],cov5[region[0]:region[1]]+0.,regionSeq,trimEnds=True)
            name = beds[c][0]+":"+beds[c][1]+"-"+beds[c][2]+beds[c][5] #chr:start-stop strand
            name = beds[c][3]
            if name in tnxDict:
                name += "."+tnxDict[name]
            #write_to_file(reactivities,name,".rx",-999,BASENAME)
            write_to_bedgraph(reactivities,beds[c][0],beds[c][1],beds[c][5],-999,BASENAME)
            enoughdata = write_map_file(reactivities,name,regionSeq,-999,BASENAME)
            ## uncomment to output coverages of individual transcripts output in addition to map and rx files 
            #if enoughdata:
            #    write_to_file(cov[region[0]:region[1]],name,".cov",1)
            #    write_to_file(cov5[region[0]:region[1]],name,".5cov",0)
            c+=1
     
    
def bedCoverageMain(bedfile,USE_BED_NAME=False):
    #pre = ["untreat1","untreat2","treat1","treat2"]
    #pre = ["vivo1","vivo2","vitro1","vitro2"]
    #filePrefixes = [DATAPATH+"/"+p+"/coverage/" for p in pre]    
    filePrefixes  = [p+"/" for p in DATAPATH.split(',')]
    #The input needs to be sorted. The best way to ensure that is to do it here.   
    sortedBed = binarySortBed(bedfile) #sorts by chromosome then start position in ascending order
    #The input can cover multiple chromosomes. Need to split up by chromosome AND strand for cov calculations. Binary sort bed does this. Just need to iterate through and define where chrom+strand starts and stops 
    lastKey = sortedBed[0][0]+sortedBed[0][5]
    start = 0
    stop = 1
    i = 1
    if len(sortedBed)==1:
        i=0
    while i < len(sortedBed):
        thisKey = sortedBed[i][0]+sortedBed[i][5]
        if thisKey!=lastKey or i==(len(sortedBed)-1):
            if thisKey!=lastKey and i==(len(sortedBed)-1): #bedfile has reached the end AND new chrom
                i-=1 #need to analyze the current chrom/strand but then the last region as well
            elif i==(len(sortedBed)-1):
                stop+=1
            lastKeyBed = sortedBed[start:stop]
            #Some of the regions could be overlapping. Need to merge them to get correct coverages,
            #but keep track of regions' relative positions in the original input with relPos.
            mergedBed = mergeBedSameChrom(lastKeyBed)
            relPos = relPosBedSameChrom(lastKeyBed)
            #print("MERGED BED:")
            #print(mergedBed)
            #print("REL POS:")
            #print(relPos)
            sequence = getSequence(mergedBed,lastKey[-1])
            cov, cov5 = getCoverage(mergedBed,filePrefixes) #cov and cov5 are 2darrays
            combineCoverage(cov,cov5,relPos,lastKeyBed,lastKey[-1],sequence,USE_BED_NAME) #lastKey[-1] gives strand info
            
            lastKey = thisKey
            start=stop
        i+=1
        stop+=1

def getSequence(mergedBed,strand):
    #write mergedBed regions to a single string. Give to BedTool
    #combine all sequences together to be the same length as cov
    seqstring = ""
    bedString = ""
    if strand=="-":
        for i in range(len(mergedBed)-1,-1,-1):
            bed = mergedBed[i]
            bedString += bed[0]+" "+bed[1]+" "+bed[2]+" "+bed[3]+" "+bed[4]+" "+bed[5]+"\n"
    else:
        for i in range(0,len(mergedBed)):
            bed = mergedBed[i]
            bedString += bed[0]+" "+bed[1]+" "+bed[2]+" "+bed[3]+" "+bed[4]+" "+bed[5]+"\n"        
    bedString.rstrip('\n') #does bedtools always return sequences in the order the bed regions are provided? I think so...
    bedObject = bedtools.BedTool(bedString, from_string=True)
    seqs = bedObject.sequence(fi=GENOME, s=True)
    fastas = open(bedObject.seqfn).read()
    fastalines = fastas.split('\n')
    for i in range(1,len(fastalines),2):
        seqstring+=fastalines[i].upper()
     
    #print(seqstring)
    return seqstring
    #Need ref sequence info because DMS only reacts with A and Cs
    #Putting the sequence of each bed region in INPUT into seqdict
    #-strand sequences get reverse complemented by pybedtools
        
if __name__ == "__main__":
    #######Command line options
    INPUT = ""
    USE_BED_NAME = False
    GENOMEFILE = 'hg38.fa'
    DATAPATHa = ""
    DATAPATHb = ""
    BASENAME = ""
    parser = ArgumentParser(
        description='''bedReactivities.py'''
    )
    parser.add_argument(
        "--input",
        "-i",
        help="histfile",
        required=True
    )
    parser.add_argument(
        "--genome",
        "-g",
        help="genome FASTA file",
        required=True
    )
    parser.add_argument(
        "--treated",
        "-a", 
        dest='treated',
        action='append', 
        default=[], 
        help="treated"
    )
    parser.add_argument(
        "--untreated",
        "-b", 
        dest='untreated',
        action='append', 
        default=[], 
        help="untreated"
    )
    parser.add_argument(
        "--output",
        "-o",
        help="output prefix (eg. SLBP)",
        required=True
    )
    parser.add_argument(
        "--name",
        "-n",
        help="USE_BED_NAME flag",
        required=False,
        default=False,
        action='store_true'
    )
    args = parser.parse_args()
    INPUT = args.input
    GENOMEFILE = args.genome
    DATAPATHa = ','.join(args.treated)
    DATAPATHb = ','.join(args.untreated)
    BASENAME = args.output + '.'
    USE_BED_NAME = args.name
    print("ARGUMENTS: ")
    print("INPUT: {}".format(INPUT))
    print("GENOMEFILE: {}".format(GENOMEFILE))
    print("DATAPATHa: {}".format(DATAPATHa))
    print("DATAPATHb: {}".format(DATAPATHb))
    print("BASENAME: {}".format(BASENAME))
    print("USE_BED_NAME: {}".format(USE_BED_NAME))
    
    DATAPATH = DATAPATHb+","+DATAPATHa
    if len(DATAPATHa.split(','))<1 or len(DATAPATHb.split(','))<1:
        print("Provided paths to treated/untreated coverage data are either not comma-separated or too few samples.")
        sys.exit()

    # if len(args)>0 and len(args)<initialArgLen:
    #     print("WARNING: Unused options", args)
    # elif len(args)>0: #options were supplied without switches - old version: bedCoverge.py file.bed sample_type
    #     INPUT = args[0]
    #     if len(args)>1:
    #         dataType = args[1].lower()
            
            
            
    ################# END command line options code
    GENOME = bedtools.example_filename(GENOMEFILE)
    bedObject = bedtools.BedTool(INPUT)
    seqs = bedObject.sequence(fi=GENOME, s=True)
    fastas = open(bedObject.seqfn).read()
    fastalines = fastas.split('\n')
    bedSeqDict = {} #building a dictionary with every input bed sequence
    #Need ref sequence info because DMS only reacts with A and Cs
    #Putting the sequence of each bed region in INPUT into seqdict
    #-strand sequences get reverse complemented by pybedtools
    for i in range(1,len(fastalines),2):
        name = fastalines[i-1].lstrip('>') #looks like chr6:26027123-26027480(-)
        bedSeqDict[name] = fastalines[i].upper()

    bedfile = open(INPUT, 'r')
    bedlines = bedfile.read().splitlines()
    bedfile.close()
    bed = np.zeros((len(bedlines),6), dtype='object')
    for i in range(0, len(bedlines)):
        if (len(bedlines[i].split()))<6:
            print("Need strand information for each bed region",bed[-1])
            sys.exit()
        #maybe should also check that start < stop
        bed[i] = bedlines[i].split()[0:6]
    if len(bedlines)<1: #input was empty
        print("Input bed file is empty.")
        sys.exit()
    bedCoverageMain(bed,USE_BED_NAME)

    

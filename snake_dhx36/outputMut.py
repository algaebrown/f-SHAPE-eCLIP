import sys
import re
def countMutationsRead(startAt,qualstring,MDtag,mutArray,covArray,strand):
    #read: [chrom,strand,pos,quals,MDtag]
    #mutArray will hold mutation counts at each pos in given segment of genome
    #posInArray tells me where in array I am with given read
    #but requires another variable to say where in genome mutArray starts
    #Based on how frequently read are mutated at the 5'end (and sometimes 3')
    #I should clip the 5'ends of + read, -3 ends of - reads################
    if strand=="+":
        trim5 = 5
        trim3 = len(qualstring)-1
    elif strand=="-":
        trim5 = 1
        trim3 = len(qualstring)-5
    #trim5 and trim3 tells me how many bases to ignore counting at each end of the given read
    item = re.compile(r"[0-9]+|\^[A-Z]+|[A-Z]") #split MD tag into sections
    delItem = re.compile(r"\^[A-Z]+")
    #originalStarAt = startAt
    #originalMutLen = len(mutArray)
    splitMD = item.findall(MDtag)
    #ex: will split "0A0C30^A50C0" into ['0', 'A', '0', 'C', '30', '^A', '50', 'C', '0']
    deletions = delItem.findall(MDtag)
    sumDel = 0
    for d in deletions: #need to consider deletions in read as num bp that need to be added to mutArray
        sumDel += len(d)-1
    posInRead = 0
    #add zeros to mutArray and covArray as needed for read's length
    toAdd = max(0,len(qualstring)+sumDel - (len(mutArray)-startAt))
    for i in range(0,toAdd):
        mutArray.append(0)
        covArray.append(0)
    for unit in splitMD:
        if unit.isalpha():
            if ord(qualstring[posInRead])>40 and posInRead>=trim5 and posInRead<trim3:
                mutArray[startAt]+=1
                covArray[startAt]+=1
            posInRead+=1
            startAt+=1
        elif unit.isdigit():
            matches = int(unit)
            for i in range(0,matches):
                #if startAt==len(covArray): #debugging
                #    print(originalStarAt, originalMutLen, toAdd)
                #    print(qualstring,MDtag)
                if posInRead>=trim5 and posInRead<trim3:
                    covArray[startAt]+=1
                startAt+=1
                posInRead+=1
        else: #unit is a deletion of the form ^[A-Z]+
            if len(unit)==2 and ord(qualstring[posInRead])>40 and posInRead>=trim5 and posInRead<trim3: #is ^N, a single deletion
                mutArray[startAt]+=1
                covArray[startAt]+=1
            startAt+=len(unit)-1 #increment genome pos not read pos for deletion
            #not counting deletions as coverage otherwise...
    return mutArray,covArray


def outputMutations(mdtagfile,outputDir): 
    ''' MAIN go through sorted MDtag file
    #MDtagfile format is: chrom(0) strand(1) pos(2) quals(3) MDtag(4)
    '''
    
    strandDict = {'+':0,'-':1}
    file = open(mdtagfile,'r')
    
    #allLines = file.read().splitlines()
    #read through one line at a time is best
    line = file.readline()
    covArray = {'+':[],'-':[]}
    mutArray = {'+':[],'-':[]} #have to keep + and - strand mutation counts separate
    mutStart = {'+':0,'-':0}
    mutChrom = {'+':'none','-':'none'}
    while line:
        #if line[0]=="@":
        #    line = file.readline()
        #    continue
        thisline = line.rstrip('\n').split()
        #print(thisline)
        quals = thisline[3]
        strand = thisline[1] #strandDict[thisline[1]]
        chrom = thisline[0]
        pos = int(thisline[2]) #mdtag genomepos is 0-based
        MDtag = thisline[4]
        
        #split cigar string on 'N'
        thisread = [chrom,strand,pos,quals,MDtag]
        mutStop = mutStart[strand]+len(mutArray[strand]) 
        if pos >= mutStop or chrom!=mutChrom[strand]:
            strand_str = ".pos" if strand == '+' else ".neg"
            #current read falls outside of current mut/cov tracking array.
            #write current mutArray/covArray to appropriate file and create new array.
            outfile1 = open(outputDir+"/"+mutChrom[strand]+strand_str+".cov",'a')
            outfile2 = open(outputDir+"/"+mutChrom[strand]+strand_str+".mut",'a')
            thisarray = mutArray[strand]
            for i in range(0,len(thisarray)):
                if covArray[strand][i] > 0:
                    outfile1.write(mutChrom[strand]+"\t"+str(i+mutStart[strand])+"\t"+str(covArray[strand][i])+"\n")
                if thisarray[i] > 0:
                    outfile2.write(mutChrom[strand]+"\t"+str(i+mutStart[strand])+"\t"+str(thisarray[i])+"\n")
            outfile1.close()
            outfile2.close()
            mutArray[strand] = []
            covArray[strand] = []
            mutStart[strand] = pos
            mutChrom[strand] = chrom
        #add this read's coverage and mutation counts of covArray and mutArray
        #start at position pos - mutStart, appending to arrays as necessary
        startAt = pos - mutStart[strand]
        thismut, thiscov  = countMutationsRead(startAt,quals,MDtag,mutArray[strand],covArray[strand],strand)
        #Dont need to return thismut, thiscov...They get modified in memory
        line = file.readline()
    file.close()

if __name__=='__main__':
    # TODO replace with pileup and count mutation
    mdtag_sorted_file = sys.argv[1]
    coverage_outputs = sys.argv[2]
    outputMutations(mdtag_sorted_file, coverage_outputs)
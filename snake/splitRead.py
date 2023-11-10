import sys
import re
def splitRead(read,cigar): #For reads with N in them
    #read is of format: [chrom,strand,pos,quals,MDtag]
    #Returns split reads in format: [[read1],[read2]]
    #Need: M+I+S+=+X occuring before any N to get qualsPos: qualPosSum
    #Numbers before M+D+N+=+X gives genome pos after N: gPosSum
    #Numbers M+D(-1)+=+X before N is the amount of MD tag I need: MDtagSum
    #example CIGARstring: "32M1026N5M"
    #example MDtag: 0A0G29T5
    #How to split MDtag? Add up int(numbers) and len(letters) (del-1) until
    #When splitting cigar, sum up without I or N separately as info for MDtag
    #May have to "split" number
    splitN = cigar.split('N')
    readItem = re.compile(r"[0-9]+[M|I|S|=|X]")
    refItem = re.compile(r"[0-9]+[M|D|N|=|X]")
    mdItem = re.compile(r"[0-9]+[M|D|=|X]")
    number = re.compile(r"([0-9]+)")
    allreads = []
    remainingMDtag = read[-1]
    qualStart = 0
    refStart = read[2]
    for i in range(0,len(splitN)-1): #for each cigar piece split by N
        thisread = read[0:2]
        thisread.append(refStart)
        spliceWithN = splitN[i]+'N' #ex: 32M1026 -> 32M1026N
        
        readSum = 0
        splitCigarString = readItem.findall(splitN[i])
        cigarList = [number.split(x)[1:] for x in splitCigarString]
        for unit in cigarList: #looks like [['32', 'M']]
            readSum+=int(unit[0])
        thisQuals = read[-2][qualStart:qualStart+readSum]
        thisread.append(thisQuals)
        qualStart +=readSum
        
        refSum = 0
        splitWithN = splitN[i]+'N' #ex: 32M1026 -> 32M1026N
        splitCigarString = refItem.findall(splitWithN)
        cigarList = [number.split(x)[1:] for x in splitCigarString]
        #print(cigarList)
        for unit in cigarList: #looks like [['32', 'M'],['1026','N']]
            refSum+=int(unit[0]) 
        refStart+=refSum
        
        mdSum = 0
        splitCigarString = mdItem.findall(splitN[i])
        cigarList = [number.split(x)[1:] for x in splitCigarString]
        for unit in cigarList: #looks like [['32', 'M']]
            mdSum+=int(unit[0]) #sum of matches/mismatches and deletions occuring before N
        item = re.compile(r"[0-9]+|\^[A-Z]+|[A-Z]") #split MD tag into sections
        splitMD = item.findall(remainingMDtag)
        fulfilledSum = False
        thisMDtag = ""
        counter = 0
        while not fulfilledSum:
        #ex: will split "0A0C30^A50C0" into ['0', 'A', '0', 'C', '30', '^A', '50', 'C', '0']
            if splitMD[counter].isalpha():
                mdSum-=1
                thisMDtag+=splitMD[counter]
                fulfilledSum = (mdSum==0)
                if fulfilledSum:
                    remainingMDtag = ""
            elif splitMD[counter].isdigit():
                numBp = int(splitMD[counter])
                sumUsed = min(mdSum,numBp)
                remainder = max(0,numBp-mdSum)
                mdSum-=sumUsed
                thisMDtag+=str(sumUsed)
                fulfilledSum = (mdSum==0)
                if fulfilledSum:
                    remainingMDtag = str(remainder) #extra 0s in mdtag could be a prob, maybe not for me  
                
            else: #item is a deletion, ex: ^A
                #Im assuming a deletion can't span an Ntype splice..that wouldnt make sense
                mdSum-=len(splitMD[counter])-1
                thisMDtag+=splitMD[counter]
                fulfilledSum = (mdSum==0)
                if fulfilledSum:
                    remainingMDtag = ""
            counter+=1
        thisread.append(thisMDtag)    
        for i in range(counter,len(splitMD)):
            remainingMDtag+= splitMD[i]
            
        allreads.append(thisread)   
    #build read for last splice read
    lastread = read[0:2]
    lastread.append(refStart)
    lastread.append(read[-2][qualStart:])
    lastread.append(remainingMDtag)
    
    allreads.append(lastread)
    return allreads

def splitReads(samfile,outputName):
    '''
    #Go through reads and split up spliced reads into two reads
    #(basically splitting up mdtag and sequence)
    #Need to output to a new file, which needs to be sorted by chrom
    # 
    samfile: .sam without header
    outputName
    '''
    strandDict = {'0':'+','16':'-','4':"unmapped"}
    outfile = open(outputName,'w')
    file = open(samfile,'r')
    #allLines = file.read().splitlines()
    #read through one line at a time is best
    line = file.readline()
    while line:
        #if line[0]=="@":
        #    line = file.readline()
        #    continue
        thisline = line.rstrip('\n').split()
        print(thisline[1])
        #print(thisline)
        quals = thisline[10]
        strand = strandDict[thisline[1]]
        if strand=="unmapped": # filter unmapped reads
            line = file.readline()
            continue
        chrom = thisline[2]
        pos = int(thisline[3])-1 #SAM is 1-based, changing to 0-based
        MDtag = thisline[16].split(':')[-1] #MD tags look like: MD:Z:25G35T1
        cigar = thisline[5] #Need: M+I-D occuring before any N to get quals 
        #ignore I and D to get pos
        item = re.compile(r"[0-9]+[M|I|D|S|N|=|X]")
        number = re.compile(r"([0-9]+)")
        #split cigar string on 'N'
        thisread = [chrom,strand,pos,quals,MDtag]
        if 'N' in cigar: #if cigar has no N in it, this will assess to 0
            splitread = splitRead(thisread,cigar) #returns [[read1],[read2],[read3],etc]
            for read in splitread:
                outfile.write(read[0]+"\t"+read[1]+"\t"+str(read[2])+"\t"+read[3]+"\t"+read[4]+"\n")
        else:
            outfile.write(chrom+"\t"+strand+"\t"+str(pos)+"\t"+quals+"\t"+MDtag+"\n")
        line = file.readline()
    file.close()
    outfile.close()

if __name__=='__main__':
    sam_file = sys.argv[1]
    mdtag_file=sys.argv[2]
    splitReads(sam_file, mdtag_file)
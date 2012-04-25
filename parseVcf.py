#!/usr/bin/env python

#options [compress_nocall, compress_reference]

import editdist
import math
import operator
import re
import dxpy

compressNoCall = True
compressReference = True

##Options for prior are ref, variant, no-call

prior = "None"
start = 0

#inputFile = open("SRR0308Combine._ALL_VARIANTS.norg.sort.txt", 'r')

inputFile = job['input']['vcf']


header = ''
count = 0


def getInfoField(fieldName, infoColumns, infoContents):
    if infoColumns.count(fieldName) > 0:
        entrySplitColumn = infoColumns.split(":")
        position = -1
        for i in range(len(entrySplitColumn)):
            if entrySplitColumn[i] == fieldName:
                position = i
                entrySplitInfo = infoContents.split(":")
                if len(entrySplitInfo) == len(entrySplitColumn):
                    return entrySplitInfo[position]
    return False


def convertPhredToProbability(phred):
    return 10.0**(-float(phred)/10.0)

def convertProbabilityToPhred(probability):
    return -10.0*math.log10(probability)


mappings_schema = [
        {"name": "chr", "type": "string"}, 
        {"name": "lo", "type": "int32"},
        {"name": "hi", "type": "int32"},
        {"name": "type", "type": "string"},     #change this type to uint once there is an abstraction method for enum
        {"name": "ref", "type": "string"},
        {"name": "alt", "type": "string"},
        {"name": "qual", "type": "int32"},
        {"name": "coverage", "type": "int32"},
        {"name": "genotypeQuality", "type": "int32"},
        
    ]

gtable = dxpy.new_dxgtable(mappings_schema, indices=[dxpy.DXGTable.genomic_range_index("chr","lo","hi",'gri')])
tableId = gtable.get_id()
gtable = dxpy.open_dxgtable(tableId)


print tableId

while count < 100:
    if count%100000 == 0:
        print count
    count += 1

    input = inputFile.readline()
    if input == '':
        break
    if input[0] == "#":
        header += input
    else:
        tabSplit = input.split("\t")
        chr = tabSplit[0]
        lo = int(tabSplit[1])
        hi = lo + len(tabSplit[3])
        ref = tabSplit[3]
        altOptions = [ref]
        altOptions.extend(tabSplit[4].split(","))
        qual = int(float(tabSplit[5]))
        formatColumn = tabSplit[7]
        infoColumns = tabSplit[8]
        genotypeQuality = 0

        coverage = re.findall("DP=(\d+);", formatColumn)
        coverage = int(coverage[0])
        
        
        if altOptions == [ref, '.']:
            #Print a reference row
            blah = 1
        else:
            genotypePossibilities = {}
            
            for i in range(9, len(tabSplit)):
                infoContents = tabSplit[9]
                genotype = getInfoField("GT", infoColumns, infoContents)
                genotypeQuality = float(getInfoField("GQ", infoColumns, infoContents))         
                if genotype != False and genotypeQuality != False:
                    if genotypePossibilities.get(genotype) == None:
                        genotypePossibilities[genotype] = float(genotypeQuality)
                    else:
                        genotypePossibilities[genotype] += float(genotypeQuality)
                else:
                    genotypeQuality = 0
            genotypePossibilities = sorted(genotypePossibilities.iteritems(), key=operator.itemgetter(1))
                
            #print genotypePossibilities[0]
            if len(genotypePossibilities) > 1:
                print genotypePossibilities[0][1]/genotypePossibilities[1][1]
             
            alt = ""
            if genotype == "0/0" or genotype == "0|0" or genotype == False:
                type  = "Ref"

            else:
                genotypeSplit = re.split("[\|\/]", genotype)
                for i in range(len(genotypeSplit)):
                #This is done to ensure the convention of placing the ref allele first
                #   in practice, it seems that all VCFs already place the ref first
                    genotypeSplit[i] = int(genotypeSplit[i])
                genotypeSplit.sort()
                for x in genotypeSplit:
                    if len(alt) > 0:
                        alt += ","
                    alt += altOptions[x]
    
                typeList = []
            
                for x in altOptions:        
                    #This categorizes the type of variant that has occurred, based on rules that describe SNPs
                    #Done per allele, if two non-ref alleles don't match type, the type is mixed.
                    type = "Unknown"
                    if x == ref:
                        type = "Ref"
                    elif editdist.distance(ref, x) == len(ref) - len(x):
                        type = "Ins"
                    elif editdist.distance(ref, x) == len(x) - len(ref):
                        type = "Del"
                    elif editdist.distance(ref, x) == 1:
                        type = "SNP"
                    elif len(ref) != len(x):
                        type = "Complex"
                    else:
                        type = "MNP"
                    if type != "Ref" and type != "Unknown":
                        typeList.append(type)
                for x in typeList[1::]:
                    if typeList[0] != x:
                        type = "Mixed"

                if type == "Mixed":
                    print "Mixed"
                    print chromosome + "\t" + str(lo) + "\t" + str(hi) + "\t" + ref + "\t" + alt + "\t" + str(qual) + "\t" + str(genotypeQuality) + "\t"
        
                
                #gtable.add_rows([chr])
        gtable.add_rows([[chr, lo, hi, type, ref, alt, qual, coverage, int(genotypeQuality)]])
        print [chr, lo, hi, type, ref, alt, qual, coverage, int(genotypeQuality)]
        
        
gtable.close(True)
q = gtable.genomic_range_query('chr1', 10000, 11000, 'enclose', 'gri')
result = gtable.get_rows(query=q)
print result
    
    
    
    
                    

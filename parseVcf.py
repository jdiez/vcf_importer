#!/usr/bin/env python

#options [compress_nocall, compress_reference]

import math
import operator
import re
import dxpy


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


def main():
    
    
    
    compressNoCall = True
    compressReference = True
    
    ##Options for prior are ref, variant, no-call

    print "1"
    
    prior = "None"
    start = 0
    
    #inputFile = open("SRR0308Combine._ALL_VARIANTS.norg.sort.txt", 'r')
    
    print job['input']['vcf']
    
    
    
    #inputFile = dxpy.open_dxfile()
    
    
    header = ''
    count = 0
    
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
    

    print "2"

    inputFile = dxpy.open_dxfile(job['input']['vcf'])
    fileIter = inputFile.__iter__()
    
    #print tableId
    #print dxpy.open_dxfile(job['input']['vcf']).describe()
    while 1:
        try:
            input = fileIter.next()
            print "3"
            #if count%100000 == 0:
                #print count
            #count += 1

            print input
            
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
                altOptions = [ref.upper()]
                altOptions.extend(tabSplit[4].upper().split(","))
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
                    #if len(genotypePossibilities) > 1:
                        #print genotypePossibilities[0][1]/genotypePossibilities[1][1]
                     
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
                        
                        #These rules determine how to characterize the type of change that has occurred
                        for x in altOptions:
                            if x == ref:
                                type = "Ref"
                            elif len(x) == len(ref) and len(ref) == 1:
                                type = "SNP"
                            elif ref in x:
                                type = "Ins"
                            elif x in ref:
                                type = "Del"
                            else:
                                type = "Complex"
                        for x in typeList[1::]:
                            if typeList[0] != x:
                                type = "Mixed"
                        print type
                        
                        #gtable.add_rows([chr])
                gtable.add_rows([[chr, lo, hi, type, ref, alt, qual, coverage, int(genotypeQuality)]])
                #print [chr, lo, hi, type, ref, alt, qual, coverage, int(genotypeQuality)]
        except StopIteration:
            break
        print "3"
            
    gtable.close(block=True)
    
    print json.dumps({'table_id':gtable.get_id()})
    
    #job['output'] = {'simplevar': {'job': reduceJobId, 'field': 'mappings'}}
    #q = gtable.genomic_range_query('chr1', 10000, 11000, 'enclose', 'gri')
    #result = gtable.get_rows(query=q)
    #print result
    
    #logging.info("Waiting for table %s to close..." % readsTable.get_id())

    job['output']['simplevar'] = dxpy.dxlink(gtable.get_id())

    
    
                    

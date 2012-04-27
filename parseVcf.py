#!/usr/bin/env python

#options [compress_nocall, compress_reference]

import math
import operator
import re
import dxpy

def main():

    compressNoCall = True
    compressReference = True
    
    print job['input']['vcf']
    header = ''

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

    inputFile = dxpy.open_dxfile(job['input']['vcf'])
    fileIter = inputFile.__iter__()
    count = 1
    
    additionalData = []
    
    while 1:
        try:
            input = fileIter.next()
            if count%100000 == 0:
                print "Processed count %i variants " % count
            count += 1
            
            if input[0] == "#":
                header += input
                #extract additional column header data
                if(input[1] != "#"):
                    tabSplit = input.split("\t")
                    additionalColumns = tabSplit[9:]
                    for i in range(len(additionalColumns)):
                        additionalColumns[i] += ":string"
                    
                        
            else:
                tabSplit = input.split("\t")
                chr = tabSplit[0]
                lo = int(tabSplit[1])
                hi = lo + len(tabSplit[3])
                ref = tabSplit[3]
                altOptions = [ref.upper()]
                altOptions.extend(tabSplit[4].upper().split(","))
                qual = tabSplit[5]
                if qual == ".":
                    type = "No-call"
                else:
                    qual = int(float(tabSplit[5]))
                type = "Unknown"
                
                formatColumn = tabSplit[7]
                infoColumns = tabSplit[8]
                
                additionalData.append(tabSplit[9:])
                genotypeQuality = 0

                coverage = re.findall("DP=(\d+);", formatColumn)
                if(len(coverage) > 0):
                    coverage = int(coverage[0])
                else:
                    coverage = 0
 
                
                if altOptions == [ref, '.']:
                    if type != "No-call":
                        type = "Ref"
                        gtable.add_rows([[chr, lo, hi, type, "", "", 0, 0, 0]])
                        print [chr, lo, hi, type, "", "", 0, 0, 0]
                else:
                    genotypePossibilities = {}
                    
                    #Find all of the genotypes 
                    for x in tabSplit[9:]:
                        genotype = getInfoField("GT", infoColumns, x)
                        genotypeQuality = float(getInfoField("GQ", infoColumns, x))
                        if genotype != False and genotypeQuality != False:
                            if genotypePossibilities.get(genotype) == None:
                                genotypePossibilities[genotype] = float(genotypeQuality)
                            else:
                                genotypePossibilities[genotype] += float(genotypeQuality)
                        else:
                            genotypeQuality = 0
                    genotypePossibilities = sorted(genotypePossibilities.iteritems(), key=operator.itemgetter(1), reverse=True)
                    genotype = genotypePossibilities[0][0]
                    genotypeQuality = genotypePossibilities[0][1]
                    if len(genotypePossibilities) > 1:
                        genotypeQuality -= genotypePossibilities[1][1]
                    
                    alt = ""
                    if genotype == "0/0" or genotype == "0|0" or genotype == False:
                        if(len(genotypePossibilities) > 1):
                            type  = "Weak"
                            genotype = genotypePossibilities[1][0]
                            genotypeQuality = 0
                            
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
                        if type != "Weak":
                            typeList = []                
                            #These rules determine how to characterize the type of change that has occurred
                            for x in altOptions:
                                if len(x) == len(ref) and len(ref) == 1:
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
                    print [chr, lo, hi, type, ref, alt, qual, coverage, int(genotypeQuality)]
                    gtable.add_rows([[chr, lo, hi, type, ref, alt, qual, coverage, int(genotypeQuality)]])
        except StopIteration:
            break
    gtable.set_details({"header":header})
    print additionalColumns
    
    gtable.close(block=True)
    print json.dumps({'table_id':gtable.get_id()})
    
    #extendedTable = gtable.extend(additionalColumns)
    #extendedTable = dxpy.extend_dxgtable(gtable.get_id(), columns=additionalColumns, indices=[dxpy.DXGTable.genomic_range_index("chr","lo","hi",'gri')])
    #extendTable.add_rows(additionalData)
    
    #print json.dumps({'table_id':extendTable.get_id()})
    
    job['output']['simplevar'] = dxpy.dxlink(gtable.get_id())
    #job['output']['extendedvar'] = dxpy.dxlink(extendTable.get_id())

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
    
                    

#!/usr/bin/env python

#options [compress_nocall, compress_reference]

import math
import operator
import re
import dxpy




def main():

    print "Running VCF to SimpleVar"

    compressNoCall = job['input']['compressNoCall']
    compressReference = job['input']['compressReference']
    storeFullVcf = job['input']['storeFullVcf']
    
    print job['input']['vcf']
    header = ''


    #These prior variables are used for keeping track of contiguous reference/no-call
    #   in the event that compressReference or compressNoCall is True
    priorType = "None"
    priorPosition = -1

    inputFile = dxpy.open_dxfile(job['input']['vcf'])
    fileIter = inputFile.__iter__()
    count = 1

    #Additional data will contain the extra format and info columns that are optional in VCF and may not be
    #   present in the VCF file. These are stored in an extended table 
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
                    tabSplit = input.strip().split("\t")
                    additionalColumns = tabSplit[7:]
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
                    if storeFullVcf:
                        mappings_schema.extend([{"name": "vcf_alt", "type": "string"}, {"name": "vcf_additional_data", "type": "string"}])
    
                    
                    #This line commented until substring index has been implemented
                    #simpleVar = dxpy.new_dxgtable(mappings_schema, indices=[dxpy.DXGTable.genomic_range_index("chr","lo","hi", 'gri'), dxpy.DXGTable.substring_index("type", "typeIndex")])
                    simpleVar = dxpy.new_dxgtable(mappings_schema, indices=[dxpy.DXGTable.genomic_range_index("chr","lo","hi", 'gri')])
                    tableId = simpleVar.get_id()
                    simpleVar = dxpy.open_dxgtable(tableId)
                    simpleVar.set_details({'header':header})
            else:
                tabSplit = input.split("\t")
                chr = tabSplit[0]
                lo = int(tabSplit[1])
                hi = lo + len(tabSplit[3])
                ref = tabSplit[3]
                
                #In VCF format, the ALT column holds possible candidate alleles. The actual call as to the
                #   variant and its zygosity is a combination of ALT and the genotype specified in the info field.
                #   We store all of the options (including ref) and calculated the actual calls later
                altOptions = [ref.upper()]
                altOptions.extend(tabSplit[4].upper().split(","))
                qual = tabSplit[5]
                type = "Unknown"
                if qual == ".":
                    type = "No-call"
                else:
                    qual = int(float(tabSplit[5]))

                formatColumn = tabSplit[7]
                infoColumn = tabSplit[8]
                
                genotypeQuality = 0
                
                coverage = re.findall("DP=(\d+);", formatColumn)
                if(len(coverage) > 0):
                    coverage = int(coverage[0])
                else:
                    coverage = 0
                    
                if altOptions == [ref, '.']:
                    if type == "No-call":
                        if compressNoCall == False:
                            entry = [chr, lo, hi, type, "", "", 0, 0, 0]
                            entry.append(tabSplit[4])
                            vcfSpecificData = ''
                            for x in tabSplit[7:]:
                                vcfSpecificData += x+"\t"
                            entry.append(vcfSpecificData.strip())
                    else:
                        type = "Ref"
                        if compressReference == False:
                            entry = [chr, lo, hi, type, "", "", 0, 0, 0]
                            entry.append(tabSplit[4])
                            vcfSpecificData = ''
                            for x in tabSplit[7:]:
                                vcfSpecificData += x+"\t"
                            entry.append(vcfSpecificData.strip())
                else:
                    #Find all of the genotypes 
                    genotypePossibilities = {}
                    for x in tabSplit[9:]:
                        genotype = getInfoField("GT", infoColumn, x)
                        genotypeQuality = float(getInfoField("GQ", infoColumn, x))
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
                            genotype = genotypePossibilities[1][0]
                            genotypeQuality = 0
                            
                    genotypeSplit = re.split("[\|\/]", genotype)
                    for i in range(len(genotypeSplit)):
                        
                    #This is done to ensure the convention of placing the ref allele first
                    #   in practice, it seems that all VCFs already place the ref first
                        genotypeSplit[i] = int(genotypeSplit[i])
                    genotypeSplit.sort()

                    #In VCF format, the prior character to a sequence change is given in some cases (Ins, Del)
                    #   we are removing this in our format, and so need to figure out which characters to filter   
                    overlap = findMatchingSequence(ref, altOptions)

                    for x in genotypeSplit:
                        if len(alt) > 0:
                            alt += "/"
                        alt += altOptions[x][overlap:]
                        if len(altOptions[x][overlap:]) == 0:
                            alt += "-"
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
                        
                    ref = ref[overlap:]
                    if len(ref) == 0:
                        ref = "-"
                    entry = [chr, lo-overlap, lo+len(ref), type, ref, alt, qual, coverage, int(genotypeQuality)]
                    if storeFullVcf:
                        entry.append(tabSplit[4])
                        vcfSpecificData = ''
                        for x in tabSplit[7:]:
                            vcfSpecificData += x+"\t"
                        entry.append(vcfSpecificData.strip())
                    simpleVar.add_rows([entry])
                if compressReference:
                    if priorType == "Ref" and type != priorType:
                        entry = [chr, priorPosition, hi, type, "", "", 0, 0, 0]
                        if storeFullVcf:
                            entry.extend([".", ""])
                        simpleVar.add_rows([entry])                        
                if compressNoCall:
                    if priorType == "No-call" and type != priorType:
                        entry = [chr, priorPosition, hi, type, "", "", 0, 0, 0]
                        if storeFullVcf:
                            entry.extend([".",""])
                        simpleVar.add_rows([entry])
                if type != priorType:
                    priorType = type
                    priorPosition = lo  
        except StopIteration:
            break
    
    simpleVar.close(block=True)
    print "SimpleVar table" + json.dumps({'table_id':simpleVar.get_id()})    
    job['output']['simplevar'] = dxpy.dxlink(simpleVar.get_id())

def findMatchingSequence(ref, altOptions):
    position = 0
    minLength = len(ref)
    for x in altOptions:
        if len(x) < minLength:
            minLength = len(x)
    for i in range(minLength):
        for x in altOptions:
            if ref[i] != x[i]:
                return i
    return minLength

def getInfoField(fieldName, infoColumn, infoContents):
    if infoColumn.count(fieldName) > 0:
        entrySplitColumn = infoColumn.split(":")
        position = -1
        for i in range(len(entrySplitColumn)):
            if entrySplitColumn[i] == fieldName:
                position = i
                entrySplitInfo = infoContents.split(":")
                if len(entrySplitInfo) == len(entrySplitColumn):
                    return entrySplitInfo[position]
    return False
    
def generateEmptyList(columns):
    result = []
    for i in range(columns):
        result.append('')
    return result
                    

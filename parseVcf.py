#!/usr/bin/env python

#options [compress_nocall, compress_reference]

import math
import re
import dxpy
import subprocess
import time




def main():

    print "Running VCF to SimpleVar"
    print job['input']['vcf']
    header = ''
    
    
    inputFile = dxpy.download_dxfile(job['input']['vcf'], 'output.vcf')
    variants_schema = [
      {"name": "chr", "type": "string"}, 
      {"name": "lo", "type": "int32"},
      {"name": "hi", "type": "int32"},
      {"name": "type", "type": "string"},
      {"name": "ref", "type": "string"},
      {"name": "alt", "type": "string"},
      {"name": "qual", "type": "int32"},
      {"name": "coverage", "type": "int32"},
      {"name": "genotypeQuality", "type": "int32"},    
         ]
    
    if 'output name' in job['input']:
        name =  job['input']['output name']
    else:
        fileName = dxpy.DXFile(job['input']['vcf']['$dnanexus_link']).describe()['name']
        name = fileName.split(".")[0]
        for x in fileName.split(".")[1:-1]:
            name += "."+x
    
    simpleVarArray = []
    if job['input']['store_full_vcf']:
        variants_schema.extend([{"name": "vcf_alt", "type": "string"}, {"name": "vcf_additional_data", "type": "string"}])

    if job['input']['store_samples_individually']:
        header = extractHeader(open('output.vcf', 'r')).strip()
        tabSplit = header.strip().split("\t")
        for x in tabSplit[9:]:
            table = dxpy.new_dxgtable(variants_schema, indices=[dxpy.DXGTable.genomic_range_index("chr","lo","hi", 'gri')])
            details = {'sample':x, 'original_contigset':job['input']['reference']}
            table.set_details(details)
            table.add_types(["SimpleVar", "gri"])
            table.rename(name + "sample " + x)
            print table.get_details()
            simpleVarArray.append(table)        
    else:
        table = dxpy.new_dxgtable(variants_schema, indices=[dxpy.DXGTable.genomic_range_index("chr","lo","hi", 'gri')])
        details = {'original_contigset':job['input']['reference']}
        setTableDetails(table, details)
        table.rename(name+".vcf")
        print table.get_details()
        simpleVarArray.append(table)
        
        
    
    command = "dx_vcfToSimplevar2"
    for x in simpleVarArray:
        command += " --table_id " + str(x.get_id())
    command += " --vcf_file output.vcf"
    if job['input']['compress_reference']:
        command += " --compress_reference"
    if job['input']['compress_no_call']:
        command += " --compress_no_call"
    if job['input']['store_full_vcf']:
        command += " --store_full_vcf"
    command += " --extract_header"
    if job['input']['store_samples_individually']:
        command += " --store_samples_individually"

    print command
    subprocess.check_call(command, shell=True)
    
    for simpleVar in simpleVarArray:
        simpleVar.close()
    
    completed = False
    while completed == False:
        completed = True
        for x in simpleVarArray:
            print x.describe()['state']
            if x.describe()['state'] != 'closed':
                completed = False
            time.sleep(2)
    result = []
    for x in simpleVarArray:    
        print "SimpleVar table" + json.dumps({'table_id':x.get_id()})
        result.append(dxpy.dxlink(x.get_id()))

    job['output']['simplevar'] = result

def setTableDetails(table, details):
    tableId = table.get_id()
    simpleVar = dxpy.open_dxgtable(tableId)
    simpleVar.set_details(details)
    

def extractHeader(vcfFile):
    header = ''
    fileIter = vcfFile.__iter__()
    while 1:
        try:
            input = fileIter.next()
            if(input[1] != "#"):
                return input
        except:
            break

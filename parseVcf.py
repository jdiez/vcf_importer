#!/usr/bin/env python

#options [compress_nocall, compress_reference]

import math
import re
import dxpy
import subprocess



def main():

    print "Running VCF to SimpleVar"
    print job['input']['vcf']
    header = ''
        
    inputFile = dxpy.download_dxfile(job['input']['vcf'], 'output.vcf')
    mappings_schema = [
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
    if job['input']['store_full_vcf']:
        mappings_schema.extend([{"name": "vcf_alt", "type": "string"}, {"name": "vcf_additional_data", "type": "string"}])
    simpleVar = dxpy.new_dxgtable(mappings_schema, indices=[dxpy.DXGTable.genomic_range_index("chr","lo","hi", 'gri')])
    tableId = simpleVar.get_id()
    simpleVar = dxpy.open_dxgtable(tableId)
    simpleVar.set_details({'original_contigset':job['input']['reference']})
        
    
    command = "dx_vcfToSimplevar --table_id %s --vcf_file output.vcf" % (tableId)
    if job['input']['compress_reference']:
        command += " --compress_reference"
    if job['input']['compress_no_call']:
        command += " --compress_no_call"
    if job['input']['store_full_vcf']:
        command += " --store_full_vcf"
    command += " --extract_header"
    print command    
    subprocess.call(command, shell=True)
    
    
    simpleVar.close(block=True)
    print "SimpleVar table" + json.dumps({'table_id':simpleVar.get_id()})    
    job['output']['simplevar'] = dxpy.dxlink(simpleVar.get_id())
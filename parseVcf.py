#!/usr/bin/env python

#options [compress_nocall, compress_reference]

import math
import re
import dxpy
import subprocess
import time
import logging

sys.path.append('/usr/local/lib/')
import magic

logging.basicConfig(level=logging.DEBUG)



def main():

    print "Running VCF to SimpleVar"
    print job['input']['vcf']
    header = ''

    inputFile = dxpy.download_dxfile(job['input']['vcf'], 'output.vcf')
    
    header = extractHeader('output.vcf').strip()
    
    variants_schema = [
      {"name": "chr", "type": "string"}, 
      {"name": "lo", "type": "int32"},
      {"name": "hi", "type": "int32"},
      {"name": "type", "type": "string"},
      {"name": "ref", "type": "string"},
      {"name": "alt", "type": "string"},
      {"name": "qual", "type": "int32"},
      {"name": "coverage", "type": "string"},
      {"name": "total_coverage", "type": "int32"},
      {"name": "genotype_quality", "type": "int32"}
         ]
    
    #if re.search("ID=AF,Number=A,Type=Float", header) != None:
    #    print "Affirmative"
    #    variants_schema.append({"name": "coverage", "type": "string"})
    #print "Negative"
    #variants_schema.extend([{"name": "total_coverage", "type": "int32"}, {"name": "genotype_quality", "type": "int32"} ])
    
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
        tabSplit = header.strip().split("\t")
        print tabSplit
        for x in tabSplit[9:]:
            table = dxpy.new_dxgtable(variants_schema, indices=[dxpy.DXGTable.genomic_range_index("chr","lo","hi", 'gri')])
            details = {'sample':x, 'original_contigset':job['input']['reference'], 'original_file':job['input']['vcf']}
            table.set_details(details)
            table.add_types(["SimpleVar", "gri"])
            table.rename(name + "sample " + x)
            print table.get_details()
            simpleVarArray.append(table)        
    else:
        table = dxpy.new_dxgtable(variants_schema, indices=[dxpy.DXGTable.genomic_range_index("chr","lo","hi", 'gri')])
        details = {'sample':name, 'original_contigset':job['input']['reference'], 'original_file':job['input']['vcf']}
        setTableDetails(table, details)
        table.add_types(["SimpleVar", "gri"])
        table.rename(name+".vcf")
        print table.get_details()
        simpleVarArray.append(table)
        
        
    
    command = "dx_vcfToSimplevar2"
    for x in simpleVarArray:
        command += " --table_id " + str(x.get_id())
    command += " --vcf_file output.vcf"
    if job['input']['compress_reference']:
        command += " --compress_reference"
    if job['input']['infer_no_call']:
        command += " --infer_no_call"
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

    job['output']['variants'] = result

def setTableDetails(table, details):
    tableId = table.get_id()
    simpleVar = dxpy.open_dxgtable(tableId)
    simpleVar.set_details(details)
    

def extractHeader(vcfFileName):
    header = ''
    with unpack_and_open(vcfFileName) as vcfFile:
        while 1:
            line = vcfFile.next()
            if line == '':
                break
            if line[1] != "#":
                return line
            
def unpack_and_open(input):
    m = magic.Magic()

    # determine compression format
    try:
        file_type = m.from_file(input)
    except:
        raise dxpy.ProgramError("Unable to identify compression format")

    print file_type

    # if uncompressed open the file and return a handle to it
    try:
        if file_type == 'ASCII text' or file_type == 'ASCII English text, with very long lines':
            return open(input)
    except:
        raise dxpy.ProgramError("Detected uncompressed input but unable to open file. File may be corrupted.")

    # if we find a tar file throw a program error telling the user to unpack it
    if file_type == 'application/x-tar':
        raise dxpy.ProgramError("Program does not support tar files.  Please unpack.")

    # since we haven't returned, the file is compressed.  Determine what program to use to uncompress
    uncomp_util = None
    if file_type == 'XZ compressed data':
        uncomp_util = 'xzcat'
    elif file_type[:21] == 'bzip2 compressed data':
        uncomp_util = 'bzcat'
    elif file_type[:20] == 'gzip compressed data':
        uncomp_util = 'zcat'
    elif file_type == 'POSIX tar archive (GNU)' or 'tar' in file_type:
        raise dxpy.ProgramError("Found a tar archive.  Please untar your sequences before importing")
    else:
        raise dxpy.ProgramError("Unsupported compression type.  Supported formats are xz, gzip, bzip, and uncompressed")

    # with that in hand, open file for reading.  If we find a tar archive then exit with error.
    try:
        with subprocess.Popen([uncomp_util, input], stdout=subprocess.PIPE).stdout as pipe:
            line = pipe.next()
        uncomp_type = m.from_buffer(line)
        print uncomp_type
    except:
        raise dxpy.ProgramError("Error detecting file format after decompression")

    if uncomp_type == 'POSIX tar archive (GNU)' or 'tar' in uncomp_type:
        raise dxpy.ProgramError("Found a tar archive after decompression.  Please untar your sequences before importing")
    elif uncomp_type != 'ASCII text':
        raise dxpy.ProgramError("After decompression found file type other than plain text")
    
    try:
        return subprocess.Popen([uncomp_util, input], stdout=subprocess.PIPE).stdout
    except:
        raise dxpy.ProgramError("Unable to open compressed input for reading")
            

            
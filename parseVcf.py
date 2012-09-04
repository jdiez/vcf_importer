#!/usr/bin/env python

#options [compress_nocall, compress_reference]

import math
import re
import dxpy
import subprocess
import time
import logging
import os
import dxpy

sys.path.append('/usr/local/lib/')
import magic

logging.basicConfig(level=logging.DEBUG)

@dxpy.entry_point('main')
def main(**job_inputs):
  
    print "Running VCF to Variants"
    print job_inputs['vcf']
    job_outputs = {}
    header = ''
    
    inputFile = dxpy.download_dxfile(job_inputs['vcf'], 'output.file')
    

    
    decompressFile('output.file')
    headerInfo = extractHeader('output.vcf')
    
    variants_schema = [
      {"name": "chr", "type": "string"}, 
      {"name": "lo", "type": "int32"},
      {"name": "hi", "type": "int32"},
      {"name": "ref", "type": "string"},
      {"name": "alt", "type": "string"},
      {"name": "qual", "type": "double"},
      {"name": "filter", "type": "string"},
      {"name": "ids", "type": "string"}
         ]
    
    description = {}
    samples = []

    print headerInfo
    print headerInfo['tags']['format']

    elevatedTags = ['format_GT', 'format_DP', 'format_AD']
    indices = [dxpy.DXGTable.genomic_range_index("chr","lo","hi", 'gri')]
    
    
    formats = {}
    infos = {}
    filters = {}
    
    for k, v in headerInfo['tags']['info'].iteritems():
        variants_schema.append({"name": "info_"+k, "type":translateTagTypeToColumnType(v)})
        description[k] = {'name' : k, 'description' : v['description'], 'type' : v['type'], 'number' : v['number']}
    
    numSamples = len(headerInfo['columns'].strip().split("\t"))
    if numSamples > 10:
      raise dxpy.AppError("The VCF file contained too many samples, can't import a VCF containing more than 10 samples")
    #For each sample, write the sample-specific columns
    for i in range(len(headerInfo['columns'].strip().split("\t")[9:])):
      #This prevents name collision in columns
      variants_schema.extend([
        {"name": "genotype_"+str(i), "type": "string"},
        {"name": "phasing_"+str(i), "type": "string"},
        {"name": "type_"+str(i), "type": "string"},
        {"name": "variation_qual_"+str(i), "type": "double"},
        {"name": "genotype_qual_"+str(i), "type": "double"},
        {"name": "coverage_"+str(i), "type": "string"},
        {"name": "total_coverage_"+str(i), "type": "int32"}
      ])
      #indices.append(dxpy.DXGTable.lexicographic_index([["type_"+str(i), "ASC"]], 'type_'+str(i)))
      samples.append(headerInfo['columns'].strip().split("\t")[9:][i])
      for k, v in headerInfo['tags']['format'].iteritems():
        if "format_"+k not in elevatedTags:
          variants_schema.append({"name": "format_"+k+"_"+str(i), "type":translateTagTypeToColumnType(v)})
        
    #for x in variants_schema:
    #  print x['name']
    
    if 'output name' in job_inputs:
        name =  job_inputs['output name']
    else:
        fileName = dxpy.DXFile(job_inputs['vcf']['$dnanexus_link']).describe()['name']
        name = fileName.split(".")[0]
        for x in fileName.split(".")[1:-1]:
            name += "."+x      
    
    details = {'samples':samples, 'original_contigset':job_inputs['reference'], 'original_file':job_inputs['vcf'], 'formats':headerInfo['tags']['format'], 'infos':headerInfo['tags']['info']}
    if headerInfo.get('filters') != {}:
      details['filters'] = headerInfo['filters']
      
    table = dxpy.new_dxgtable(variants_schema, indices=indices)
    table.set_details(details)
    table.add_types(["Variants", "gri"])
    table.rename(name)

    command = "dx_vcfToVariants2"
    command += " --table_id " + str(table.get_id())
    command += " --vcf_file output.vcf"
    if job_inputs['compress_reference']:
        command += " --compress_reference"
    if job_inputs['infer_no_call']:
        command += " --infer_no_call"
    if job_inputs['compress_no_call']:
      command += " --compress_no_call"
    
    print command
    subprocess.check_call(command, shell=True)
    

    table.close()
    result = dxpy.dxlink(table.get_id())
    
    job_outputs['variants'] = result
    return job_outputs

def setTableDetails(table, details):
    tableId = table.get_id()
    table = dxpy.open_dxgtable(tableId)
    table.set_details(details)
    

def extractHeader(vcfFileName):
  result = {'columns': '', 'tags' : {'format' : {}, 'info' : {} }, 'filters' : {}}
  vcfFile = open(vcfFileName, 'r')
  while 1:
    line = vcfFile.next().strip()
    tag = re.findall("ID=(\w+),", line)
    if len(tag) > 0:
      tagType = ''
      if line.count("##FORMAT") > 0:
        tagType = 'format'
      elif line.count("##INFO") > 0:
        tagType = 'info'
      elif line.count("##FILTER") > 0:
        result['filters'][re.findall("ID=(\w+),", line)[0]] = re.findall('Description="(.*)"', line)[0]

      typ = re.findall("Type=(\w+),", line)
      if tagType != '':
        number = re.findall("Number=(\w+)", line)
        description = re.findall('Description="(.*)"', line)
        if len(number) == 0:
          number = ['.']
        if len(description) == 0:
          description = ['']
        result['tags'][tagType][tag[0]] = {'type':typ[0], 'description' : description[0], 'number' : number[0]}
    if line[0] == "#" and line[1] != "#":
      result['columns'] = line.strip()        
    if line == '' or line[0] != "#":
        break
  return result


def checkUncompressedFile(input):
    m = magic.Magic()

    # determine compression format
    try:
        file_type = m.from_file(input)
    except:
        raise dxpy.ProgramError("Unable to identify compression format")

    print file_type

    # if uncompressed open the file and return a handle to it
    try:
        if 'ASCII' in file_type:
            return True
        else:
          return False
    except:
        raise dxpy.ProgramError("Detected uncompressed input but unable to open file. File may be corrupted.")
  
            
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
        if 'ASCII' in file_type:
            return open(input)
    except:
        raise dxpy.ProgramError("Detected uncompressed input but unable to open file. File may be corrupted.")

    # if we find a tar file throw a program error telling the user to unpack it
    print file_type
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

def decompressFile(inputFile):
    m = magic.Magic()
    # determine compression format
    try:
        file_type = m.from_file(inputFile)
    except:
        raise dxpy.ProgramError("Unable to identify compression format")

    print file_type
    # if uncompressed open the file and return a handle to it
    if 'ASCII' in file_type:
        subprocess.call("mv output.file output.vcf", shell=True)
    else:
      # if we find a tar file throw a program error telling the user to unpack it
      print file_type
      if file_type == 'application/x-tar':
          raise dxpy.ProgramError("Program does not support tar files.  Please unpack.")  
      uncomp_util = None
      if file_type == 'XZ compressed data':
          subprocess.call("mv output.file output.vcf.xz", shell=True)
          subprocess.call('xz -d output.xz', shell=True)
      elif file_type[:21] == 'bzip2 compressed data':
          subprocess.call("mv output.file output.vcf.bz2", shell=True)
          subprocess.call('bzip2 -d output.vcf.bz2', shell=True)
      elif file_type[:20] == 'gzip compressed data':
          subprocess.call('mv output.file output.vcf.gz', shell=True)
          subprocess.call('gzip -d output.vcf.gz', shell=True)
      elif file_type == 'POSIX tar archive (GNU)' or 'tar' in file_type:
          raise dxpy.ProgramError("Found a tar archive.  Please untar your sequences before importing")
      else:
          raise dxpy.ProgramError("Unsupported compression type.  Supported formats are xz, gzip, bzip, and uncompressed")


def translateTagTypeToColumnType(tag):
  if tag['type'] == "Flag":
    return "boolean"
  if tag['number'] != '1':
    return 'string'
  if tag['type'] == "Integer":
    return 'int32'
  if tag['type'] == "Float":
    return "double"
  return "string"


dxpy.run()
            

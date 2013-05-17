#!/usr/bin/env python
#
# Copyright (C) 2013 DNAnexus, Inc.
#
# This file is part of vcf_importer.
#
#   Licensed under the Apache License, Version 2.0 (the "License"); you may not
#   use this file except in compliance with the License. You may obtain a copy
#   of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
#   WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
#   License for the specific language governing permissions and limitations
#   under the License.

#options [compress_nocall, compress_reference]

import math
import re
import dxpy
import subprocess
import time
import logging
import os
import dxpy

import magic

@dxpy.entry_point('main')
def main(**job_inputs):

    job_outputs = {}
    header = ''

    print "Downloading input VCF file"
    inputFile = dxpy.download_dxfile(job_inputs['vcf'], 'output.file')

    decompressFile('output.file')

    print "Constructing table schema"
    elevatedTags = ['format_GT', 'format_DP', 'format_AD']
    headerInfo = extractHeader('output.vcf', elevatedTags)

    variants_schema = [
      {"name": "chr", "type": "string"},
      {"name": "lo", "type": "int32"},
      {"name": "hi", "type": "int32"},
      {"name": "ref", "type": "string"},
      {"name": "alt", "type": "string"},
      {"name": "qual", "type": "double"},
      {"name": "ids", "type": "string"}
         ]

    description = {}
    samples = []

    if headerInfo.get('filters') != {}:
      variants_schema.append({"name": "filter", "type": "string"})
    indices = [dxpy.DXGTable.genomic_range_index("chr","lo","hi", 'gri')]

    formats = {}
    infos = {}
    filters = {}

    for k, v in headerInfo['tags']['info'].iteritems():
        variants_schema.append({"name": "info_"+k, "type":translateTagTypeToColumnType(v)})
        description[k] = {'name' : k, 'description' : v['description'], 'type' : v['type'], 'number' : v['number']}

    numSamples = len(headerInfo['columns'].strip().split("\t")[9:])
    if numSamples > 10:
      raise dxpy.AppError("The VCF file contained too many samples, can't import a VCF containing more than 10 samples")
    if job_inputs['searchable_ids']:
      indices.append(dxpy.DXGTable.lexicographic_index([
        dxpy.DXGTable.lexicographic_index_column("ids", True, False),
        dxpy.DXGTable.lexicographic_index_column("chr"),
        dxpy.DXGTable.lexicographic_index_column("lo"),
        dxpy.DXGTable.lexicographic_index_column("hi")], "search"))
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


    if 'output_name' in job_inputs:
        name =  job_inputs['output_name']
    else:
        fileName = dxpy.DXFile(job_inputs['vcf']['$dnanexus_link']).describe()['name']
        name = fileName.split(".")[0]
        for x in fileName.split(".")[1:-1]:
            name += "."+x

    details = {'samples':samples, 'original_contigset':job_inputs['reference'], 'original_file':job_inputs['vcf'], 'formats':headerInfo['tags']['format'], 'infos':headerInfo['tags']['info'], 'alts':headerInfo['tags']['alt']}
    if headerInfo.get('filters') != {}:
      details['filters'] = headerInfo['filters']

    table = dxpy.new_dxgtable(variants_schema, indices=indices)
    table.set_details(details)
    types = ["Variants", "gri"]
    if 'additional_types' in job_inputs:
      for x in job_inputs['additional_types'].split(","):
        if x != '':
          types.append(x)
    table.add_types(types)

    if 'tags' in job_inputs:
        table.add_tags(job_inputs['tags'])
    if 'properties' in job_inputs:
        table.set_properties(job_inputs['properties'])


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
    command += " --encoding "+job_inputs["file_encoding"]

    print "Importing variants by running:", command
    try:
      subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
      try:
        errorData = open("AppError.txt", 'r').read()
        raise dxpy.AppError(errorData)
      except IOError:
        raise dxpy.AppError("An unknown error occurred. Please check the log file")

    attach_empty_trackspec(table)
    table.close()
    result = dxpy.dxlink(table.get_id())

    job_outputs['variants'] = result
    return job_outputs

def extractHeader(vcfFileName, elevatedTags):
  result = {'columns': '', 'tags' : {'format' : {}, 'info' : {}, 'alt': {} }, 'filters' : {}}
  vcfFile = open(vcfFileName, 'r')
  while 1:
      line = vcfFile.next().strip()
      tag = re.findall("ID=([^,]*),", line)
      tagType = ''
      if len(tag) > 0:
          tagType = ''
      if line.count("##FORMAT") > 0:
          tagType = 'format'
      elif line.count("##INFO") > 0:
          tagType = 'info'
      elif line.count("##ALT") > 0:
          tagType = 'alt'
      elif line.count("##FILTER") > 0:
          result['filters'][re.findall("ID=([^,]*),", line)[0]] = re.findall('Description="(.*)"', line)[0]

      typ = re.findall("Type=([^,]*),", line)
      if tagType == 'alt':
          if tagType == 'alt':
              description = re.findall('Description="(.*)"', line)
              result['tags'][tagType][tag[0]] = {'description' : description[0]}
      elif tagType != '':
          number = re.findall("Number=([^,]*),", line)
          description = re.findall('Description="(.*)"', line)
          if len(number) == 0:
              number = ['.']
          if len(description) == 0:
              description = ['']
          if tagType+"_"+tag[0] not in elevatedTags:
              result['tags'][tagType][tag[0]] = {'type':typ[0], 'description' : description[0], 'number' : number[0]}
      if line == '' or line[0] != "#":
          break

      if line[0] == "#":
          try:
              if line[1] != "#":
                  if result['columns'] == '':
                      result['columns'] = line.strip()
                  else:
                      raise dxpy.AppError("Could not identify a unique header line signifying the VCF samples.")
          except:
              print "ignoring single comment line"

  if len(result['tags']['format']) + len(result['tags']['info']) == 0:
      raise dxpy.AppError("Could not find any VCF-specific format or info tags. Is this a valid VCF?")

  if result['columns'] == '':
      raise dxpy.AppError("Couldn't find a header line singifying the VCF samples. Is this a valid VCF?")

  return result


def decompressFile(inputFile):
    print "Detecting input file type"
    m = magic.Magic()
    # determine compression format
    try:
        file_type = m.from_file(inputFile)
    except:
        raise dxpy.AppError("Unable to identify input file format")

    print "Input file type is:", file_type
    # if uncompressed open the file and return a handle to it
    if 'ASCII' in file_type:
        subprocess.call("mv output.file output.vcf", shell=True)
    else:
      # if we find a tar file throw a program error telling the user to unpack it
      if file_type == 'application/x-tar':
          raise dxpy.AppError("Unsupported compression type (tar). Supported compression formats are xz, gzip, and bzip2")
      if 'Zip archive data' in file_type:
          raise dxpy.AppError("Unsupported compression type (zip). Supported compression formats are xz, gzip, and bzip2")
      if file_type == 'XZ compressed data':
          print "Uncompressing data with xz"
          subprocess.call("mv output.file output.vcf.xz", shell=True)
          subprocess.call('xz -d output.vcf.xz', shell=True)
      elif file_type[:21] == 'bzip2 compressed data':
          print "Uncompressing data with bzip2"
          subprocess.call("mv output.file output.vcf.bz2", shell=True)
          subprocess.call('bzip2 -d output.vcf.bz2', shell=True)
      elif file_type[:20] == 'gzip compressed data':
          print "Uncompressing data with gzip"
          subprocess.call('mv output.file output.vcf.gz', shell=True)
          subprocess.call('gzip -d output.vcf.gz', shell=True)
      elif file_type == 'POSIX tar archive (GNU)' or 'tar' in file_type:
          raise dxpy.AppError("Unsupported compression type (tar). Supported compression formats are xz, gzip, and bzip2")
      else:
          raise dxpy.AppError("Unknown file format.  Supported formats are vcf (text, uncompressed), or compressed with xz, gzip, or bzip2")

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

def attach_empty_trackspec(table):
  description = table.describe(False, True)

  # Extract row count from unclosed table
  rows = 0
  for part_id, part_description in description['parts'].items():
    rows += part_description['length']

  # Calculate genome length
  genome_length = 0
  try:
    genome_length = sum(dxpy.DXRecord(description['details']['original_contigset']).describe(False, True)['details']['contigs']['sizes'])
  except:
    return

  # Skip extreme cases
  if genome_length < 1 or rows < 1:
    return

  fetch_limit = 10000
  typical_browser_width = 1200
  bpp = int((genome_length * fetch_limit) / (rows * typical_browser_width))

  # Skip extremely high numbers
  if bpp >= genome_length:
    return

  print "Attaching trackspec so that no data are displayed beyond",bpp,"bases per pixel"
  table.add_types(["TrackSpec"])
  details = description['details']
  details['representations'] = [[0, {"type": "variants", "source": dxpy.dxlink(table.get_id())}], [bpp, {"type": "empty"}]]
  table.set_details(details)

dxpy.run()

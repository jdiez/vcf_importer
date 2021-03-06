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

import os, sys, unittest, json, subprocess

import dxpy, dxpy.app_builder
from dxpy.exceptions import *

src_dir = os.path.join(os.path.dirname(__file__), "..")
test_resources_dir = os.path.join(src_dir, "test", "resources")



def makeInputs():
    #try:
    #    contigset_importer = dxpy.DXApplet(dxpy.find_data_objects(classname="applet", properties={"name": "fasta_contigset_importer"}).next()['id'])
    #except StopIteration:
    #    raise Exception("fasta_contigset_importer not found, please upload them")
    #
    #genome_archive = dxpy.upload_local_file(os.path.join(test_resources_dir, "hg19_chr22.fa.xz"), wait_on_close=True)
    #contigset_importer_input = {"name": "hg19_chr22", "sequence_file": dxpy.dxlink(genome_archive)}
    #print "Running fasta_contigset_importer with", contigset_importer_input
    #job = contigset_importer.run(contigset_importer_input)
    #job.wait_on_done()
    #contig_set = job.describe()["output"]["contig_set"]

    vcf = dxpy.upload_local_file(os.path.join(test_resources_dir, "SRR10022_GATK_Subsample.vcf"), wait_on_close=True)
    program_input = {"vcf": dxpy.dxlink(vcf), "compress_reference":False, "compress_no_call":True, "infer_no_call": False, "reference": {"$dnanexus_link":"record-9ykz7KQ00006B3PXk1b00005"}}
    #program_input = {"vcf": dxpy.dxlink(vcf), "compress_reference":False, "compress_no_call":True, "infer_no_call": False, "reference": {"$dnanexus_link":"record-9zV2FBQ0000293088JZ00005"}}
    print program_input
    return program_input

class TestMyApp(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.base_input = makeInputs()
        bundled_resources = dxpy.app_builder.upload_resources(src_dir)
        cls.program_id = dxpy.app_builder.upload_applet(src_dir, bundled_resources, overwrite=True)

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_default_vcf(self):

        job = dxpy.DXApplet(self.program_id).run(self.base_input)
        print "Waiting for job to complete"
        job.wait_on_done()
        print json.dumps(job.describe()["output"])


if __name__ == '__main__':
    unittest.main()

#!/usr/bin/env python
import os, sys, unittest, json, subprocess

import dxpy, dxpy.app_builder
from dxpy.exceptions import *

src_dir = os.path.join(os.path.dirname(__file__), "..")
test_resources_dir = os.path.join(src_dir, "test", "resources")



def makeInputs():
    try:
        contigset_importer = dxpy.DXApplet(dxpy.find_data_objects(classname="applet", properties={"name": "fasta_contigset_importer"}).next()['id'])
    except StopIteration:
        raise Exception("fasta_contigset_importer not found, please upload them")
    
    genome_archive = dxpy.upload_local_file(os.path.join(test_resources_dir, "hg19_chr22.fa.xz"), wait_on_close=True)
    contigset_importer_input = {"name": "hg19_chr22", "sequence_file": dxpy.dxlink(genome_archive)}
    print "Running fasta_contigset_importer with", contigset_importer_input
    job = contigset_importer.run(contigset_importer_input)
    job.wait_on_done()
    contig_set = job.describe()["output"]["contig_set"]

    vcf = dxpy.upload_local_file(os.path.join(test_resources_dir, "variants.vcf"), wait_on_close=True)
    program_input = {"vcf": dxpy.dxlink(vcf), "compress_reference":True, "compress_no_call":True, "reference": contig_set, "store_samples_individually":True}        
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

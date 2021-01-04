#!/usr/bin/env python3
import pysam
from pysam import VariantFile
import argparse
import os
import xml.etree.ElementTree as ET
import subprocess
vcf='DBA_2J.mgp.v5.snps.dbSNP142.chr.vcf'
myvcf=VariantFile(vcf,"r")
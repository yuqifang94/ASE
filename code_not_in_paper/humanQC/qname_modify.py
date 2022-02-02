import pysam
from itertools import izip
bam1=pysam.AlignmentFile('H1_WGBS_1_Arioc_WGBS_st.bam','rb')
outbam = pysam.AlignmentFile("H1_WGBS_1_Arioc_WGBS_st_rq.bam", "wb", template = bam1)
for read in bam1.fetch(until_eof=True):
	read.query_name=read.query_name.split(" ")[0]
	outbam.write(read)
outbam.close()

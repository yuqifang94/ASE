import glob
import pandas as pd
def dupExtract(fn):
    f= open(fn,'r')
    datIn=f.readlines()
    f.close()
    return pd.DataFrame(data={'Sample':[fn.split('/')[-1].split('.')[0]],"duplication":[float(datIn[-3].split('\t')[-2])]})
input_dir="../downstream/input/mouse_analysis/dedup_report/"
fn_all=glob.glob(input_dir+"*.sort.dup_report.txt", recursive=False)
dupDt = pd.DataFrame()
for fn  in fn_all:
    dupDt=pd.concat([dupDt,dupExtract(fn)])
dupDt.to_csv('../downstream/output/mouse_analysis/duplication_rate_fn.csv',index=False)
#Example file
# ## htsjdk.samtools.metrics.StringHeader
# # MarkDuplicates INPUT=[mm10_heart_day10_5_merged1.sort.bam] OUTPUT=mm10_heart_day10_5_merged1.sort.dup.bam METRICS_FILE=mm10_heart_day10_5_merged1.sort.dup_report.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true READ_NAME_REGEX=null TMP_DIR=[/home-4/yfang27@jhu.edu/scratch_feinberg/yfang/allele_agnostic_mouse_all/CPEL/CPEL_agnostic/cpelasm/bam/picard_tmp] QUIET=true VALIDATION_STRINGENCY=LENIENT    MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 SORTING_COLLECTION_SIZE_RATIO=0.25 TAG_DUPLICATE_SET_MEMBERS=false REMOVE_SEQUENCING_DUPLICATES=false TAGGING_POLICY=DontTag CLEAR_DT=true DUPLEX_UMI=false ADD_PG_TAG_TO_READS=true DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES PROGRAM_RECORD_ID=MarkDuplicates PROGRAM_GROUP_NAME=MarkDuplicates OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 MAX_OPTICAL_DUPLICATE_SET_SIZE=300000 VERBOSITY=INFO COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json USE_JDK_DEFLATER=false USE_JDK_INFLATER=false
# ## htsjdk.samtools.metrics.StringHeader
# # Started on: Wed Aug 19 17:01:33 EDT 2020

# ## METRICS CLASS        picard.sam.DuplicationMetrics
# LIBRARY UNPAIRED_READS_EXAMINED READ_PAIRS_EXAMINED     SECONDARY_OR_SUPPLEMENTARY_RDS  UNMAPPED_READS  UNPAIRED_READ_DUPLICATES        READ_PAIR_DUPLICATES    READ_PAIR_OPTICAL_DUPLICATES    PERCENT_DUPLICATION     ESTIMATED_LIBRARY_SIZE
# Unknown Library 832583476       0       0       0       123018435       0       0       0.147755
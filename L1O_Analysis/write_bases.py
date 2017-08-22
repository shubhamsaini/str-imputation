### python write_bases.py imputed.str SS0012978 SS0012978.imputeResult.txt

from cyvcf2 import VCF
import sys
import re

vcf = VCF(sys.argv[1], samples=sys.argv[2])
file = open(sys.argv[3],'a')
for v in vcf:
	pos = ":".join([str(v.CHROM), str(v.POS)])
	gt_bases = v.gt_bases[0]
	gt_bases = re.split('/|\|',gt_bases)
	gt_bases_len = [str(len(i)) for i in gt_bases]
	final_list = [pos]+gt_bases_len+["\n"]
	file.write("\t".join(final_list))
file.close()
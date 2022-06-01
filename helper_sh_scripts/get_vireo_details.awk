#!/usr/bin/awk -f

# The input is GT_donors.vireo.vcf.gz

BEGIN {
	OFS=" "
	FS="\t"
}

NR==1 {
	print "CHR", "POS", "REF", "ALT", "AD_total", "DP_total", "AD_1", "DP_1", "AD_2", "DP_2", "AD_3", "DP_3", "AD_4", "DP_4", "AD_5", "DP_5", "AD_6", "DP_6"
}
$0 !~ /^#/ { 
	l=""
	ad_id=0
	dp_id=0
	for (i=1; i<=NF;i++) {
		if (i == 8) {
			a=gensub(".*AD=([0-9]+).*", "\\1", 1, $i)
			b=gensub(".*DP=([0-9]+).*", "\\1", 1, $i)
			l = l OFS a OFS b
		}
		# identify where is AD and DP column
		else if (i == 9) {
			split($i, id, ":")
			for (j in id) {
				if (id[j] == "AD") {
					ad_id=j
				}
				else if (id[j] == "DP") {
					dp_id=j
				}
			}
		}
		else if (i >= 10) {
			split($i, x, ":")
			l = l OFS x[ad_id] OFS x[dp_id]
		}
		else if (i == 3 || i == 6 || i == 7) {
			continue
		}		
		else {
			l = l OFS $i
		}
	}
	print l
}

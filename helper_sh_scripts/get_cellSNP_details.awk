#!/usr/bin/awk -f

# Input is cellSNP.cells.vcf.gz
BEGIN {
	OFS=" "
	FS="\t"
}

NR==1 {
	print "CHR", "POS", "REF", "ALT", "Total_counts", "homo_ref_cells", "homo_alt_cells", "het_cells", "non_empty_cells"
}
NR>38 { 
	l=""
	homo_ref=0
	homo_alt=0
	het=0
	gt_id=0
	ad_id=0
	dp_id=0
	for (i=1; i<=NF;i++) {
		if (i == 8) {
			a=gensub(".*AD=([0-9]+).*", "\\1", 1, $i)
			b=gensub(".*DP=([0-9]+).*", "\\1", 1, $i)
			c=gensub(".*OTH=([0-9]+).*", "\\1", 1, $i)
			s=b+c
		}
		# identify where is GT, AD and DP column
		else if (i == 9) {
			split($i, id, ":")
			for (j in id) {
				if (id[j] == "AD") {
					ad_id=j
				}
				else if (id[j] == "DP") {
					dp_id=j
				}
				else if (id[j] == "GT") {
					gt_id=j
				}
			}
		}
		else if (i>=10 && $i ~ /0\/0/) {
			homo_ref++
		}
		else if (i>=10 && $i ~ /1\/0/) {
			het++
		}
		else if (i>=10 && $i ~ /1\/1/) {
			homo_alt++
		}
		else if (i==1 || i==2 || i==4 || i==5) {
			l = l OFS $i
		}	
		else {
			continue
		}
	}
	t_c = homo_ref + homo_alt + het
	printf "%s %d %d %d %d %d\n", l, s, homo_ref, homo_alt, het, t_c
}

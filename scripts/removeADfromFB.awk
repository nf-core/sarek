# this script removes multiple AD fields from FreeBayes genotypes
# the issues is that for multi-allelic calls vcfflatten writes out the most likely allele, but leaves the genotype as multi-allelic
# so we will have calls with genotypes like 0/3 or 5/5 instead of 0/1 and 1/1 
# also the AD field is separated to multiple fields
# 
# use like 
# awk -f removeADfromFB.awk infile.vcf > outfile.vcf 


function fixGT(gtStr) {
# this makes a multi-allele genotype string to biallelic (like 0/3 or 5/5 changed to 0/1 or 1/1)
		split(gtStr,GT,"/") 
		a1 = GT[1] ~ /0/ ? "0" : "1"
		a2 = GT[2] ~ /0/ ? "0" : "1"
	
		return a1 "/" a2
}

function exciseAD(AD_index,gt) {
		newGT = ""
		for(i=2;i<=length(gt); i++) {
				if(i != AD_index)
					newGT = newGT gt[i]
				if(i != AD_index && i != length(gt))
					newGT = newGT ":"
		}

		return newGT
}

/^#/{print}
!/^#/{
		for(i=1;i<=8;i++) printf("%s\t",$i);
		# 9th field is the genotype format: remove AD
		split($9, format, ":")
		formatString = ""
		for(i=1;i<=length(format);i++) {
				# print out fields, but not AD - also remember the index of AD
				if(format[i] !~ /AD/)
					formatString = formatString format[i]
				else 
					ADidx = i
				# add ":" only if it is not the last field (and not AD of course)
				if(i != length(format) && format[i] !~ /AD/)
						formatString = formatString ":"
		}
		printf("%s\t",formatString)

		# excise AD from the next genotype, and make a [01]/[01] from [0-9]/[0-9] genotypes
		# we just do not know which one is the normal and the tumour genotype, but we have to do the
		# same for both nevertheless

		split($10,gt,":")
		genotype = fixGT(gt[1]) ":" exciseAD(ADidx,gt)
		printf("%s\t",genotype)

		split($11,gt,":")
		genotype = fixGT(gt[1]) ":" exciseAD(ADidx,gt)
		print genotype

}

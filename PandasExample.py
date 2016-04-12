#!/usr/bin/python
import os
import re
import pandas as pd
import numpy as np

#GENE FORMAT
#PIK3CA	3	178866311	178952497


#BIO-VCF
#


biovcf = "bio-vcf -i --num-threads=2"
vocabulary = {"None":-1, "clean":0, "sequence_feature":0, "synonymous_variant":0, "intron_variant":0, "5_prime_UTR_variant":0.5, "non_coding_exon_variant":0.5, "missense_variant":1, "splice_region_variant":2, "inframe_deletion":3, "stop_gained":4, "nonsense_mediated_decay":4, "frameshift_variant":5}
toselect =["missense_variant", "splice_region_variant", "inframe_deletion", "stop_gained", "nonsense_mediated_decay", "frameshift_variant"]


variantfile="PIK3CA.vcf"
#"iap_all.filtered_variants_snpEff_snpSift_Cosmicv72_GoNLv5.vcf"

genelist=open('genelist.txt','r')


min_cov = 10
min_pnr = 0.15
max_pop = 0.05




for gene in genelist:

	thisgene = gene.strip().split('\t')
	subfile = thisgene[0]+".vcf"

	print(thisgene)

	# make gene subfile
	if not os.path.exists(subfile):
		print("Making Gene file")
		os.system(biovcf+" --filter \'rec.chrom==\"%s\" and rec.pos>=%s and rec.pos<=%s\'"%(thisgene[1],thisgene[2],thisgene[3])+ " < "+variantfile+" > "+subfile)

	# make INFO, GT, DP and PNR files
	if not os.path.exists(subfile+".gt"):
		print("Making GT/DP/PNR files")
		formatting = " --set-header \"chr,pos,#samples\" < "+subfile
		os.system(biovcf+" --seval \'s.gt\'"+formatting+" > "+subfile+".gt")
		os.system(biovcf+" --seval \'s.dp\'"+formatting+" > "+subfile+".dp")
		os.system(biovcf+" --seval \'tot=s.ad.reduce(:+) ; (tot-s.ad[0].to_f)/tot\'"+formatting+" > "+subfile+".pnr")
	
	#print("DONE")
	#break
	gtdf = pd.read_table(subfile+".gt", sep='\t').iloc[0:,2:]
	dpdf = pd.read_table(subfile+".dp", sep='\t').iloc[0:,2:]
	pnrdf = pd.read_table(subfile+".pnr", sep='\t').iloc[0:,2:]
	
	#print(gtdf.iloc[0:,133:140])
	gtdf[dpdf  < min_cov] = np.nan
	pnrdf[dpdf < min_cov] = np.nan
	pnrdf[gtdf=="./."] = np.nan
	#print(gtdf.iloc[0:,133:140])
	gtdf[pnrdf < min_pnr] = "0/0"
	
	#print(gtdf.iloc[0:,133:140])
	#break

	#samcov = dpdf.mean(0)
	#poscov = dpdf.mean(1)

	posdict = []
	varcounter = 0
	header = None

	subreader = open(subfile,'r')
	for variant in subreader:
		if variant.startswith("##"):
			continue
		if variant.startswith("#"):
			header=variant.strip().split("\t")
			continue
		
		data = dict(zip(header, variant.strip().split("\t")))



		posdict.append({})

		posdict[varcounter]["POS"] = data["POS"]
		posdict[varcounter]["REF"] = data["REF"]
		posdict[varcounter]["ALT"] = data["ALT"]

		alt = data["ALT"].split(",")
		popfreq=[0.0]*len(alt)
		
		ExAC_search = re.search('dbNSFP_ExAC_AF=(.*?);', data["INFO"])
		if ExAC_search:
			#print(ExAC_search.group(1))
			popfreq = [float(x) for x in ExAC_search.group(1).split(",")]
		
		ExAC_search = re.search('dbNSFP_ExAC_Adj_AF=(.*?);', data["INFO"])
		if ExAC_search:
			#print(ExAC_search.group(1))
			popfreq = [float(x) for x in ExAC_search.group(1).split(",")]
			
		GoNL_search = re.search('GoNLv5_Freq=(.*?);', data["INFO"])
		if GoNL_search:
			#print(GoNL_search.group(1))
			popfreq = [float(x) for x in GoNL_search.group(1).split(",")]

		posdict[varcounter]["PopFreq"] = popfreq
		posdict[varcounter]["Effects"] = ["None"]*len(alt)


		ann_search = re.search('ANN=(.*?);', data["INFO"])
		ann=None
		if ann_search:
			ann = ann_search.group(1).split(",")
			for pred in ann:
				items = pred.split("|")
				allele = items[0]
				effects = items[1].split("&")
				for effect in effects:
					if effect not in vocabulary:
						print effect
						print ann
					else:
						if vocabulary[effect] > vocabulary[posdict[varcounter]["Effects"][alt.index(allele)]]:
							posdict[varcounter]["Effects"][alt.index(allele)] = effect

		varcounter+=1
			
	#print(posdict)
	subreader.close()
	
	#os.system("rm "+' '.join(subfile+".gt",subfile+".dp",subfile+".pnr"))
	
	scount = 0
	varcount = ["None"]*len(gtdf.iloc[0,])

	# print(varcount)
	# Determine per gene per sample load and per gene popultaion frequency
	for i in range(0,len(posdict)):
		# Check if minumum coverage is good
		if (np.median(dpdf.iloc[i,]) > min_cov):
			# Check if population frequency within bounds
			if (np.max(posdict[i]["PopFreq"]) < max_pop):
				# Compute number of samples
				this_scount = np.sum(dpdf.iloc[i,] >= min_cov)
				
				scount = max(this_scount, scount)

				

				gts = gtdf.iloc[i,]
				#print(gts)
				for j in range(0, len(gts)):
					effect = "None"
					if gts[j] == "0/0":
						effect = "clean"

					if gts[j] == "0/1" or gts[j] == "1/1":
						effect = posdict[i]["Effects"][0]

					if gts[j] == "0/2" or gts[j] == "2/2":
						effect = posdict[i]["Effects"][1]

					if gts[j] == "1/2":
						print("1/2 detected: ",i,j)
						effect = posdict[i]["Effects"][1]
						# FIXME

					if vocabulary[effect] > vocabulary[varcount[j]]:
						varcount[j] = effect

				
	#print(varcount)
	#print(scount)
	measured = [x!="None" for x in varcount]
	print(np.sum(measured))
	#print(np.sum([x=="missense_variant" for x in varcount]))
	affected = [x in toselect for x in varcount ]
	print(np.sum(affected))

	samplenames = list(gtdf.columns.values)
	tumors = ["R" not in x for x in samplenames]
	#print(dict(zip(samplenames,tumors)))
	print(dict(zip(samplenames,affected)))
	print(dict(zip(samplenames,varcount)))

	print(np.sum([a and b for a, b in zip(tumors, measured)]))
	print(np.sum([a and b for a, b in zip(tumors, affected)]))

	#hot = [a or b for a, b in zip(tumors, affected)]
	#hotpos = [i for i, x in enumerate(hot) if x]
	#for h in hotpos:
	#	print(h)
	#	print(samplenames[h])
	#print(dict(zip(samplenames,varcount)))
	#print(affected/(measured*1.0))

genelist.close()




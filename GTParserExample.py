infile = open("PIK3CA.vcf",'r')

for line in infile:
	header=""
	if line.startswith("##"):
		continue
	if line.startswith("#"):
		header = lines.strip().split('\t')
		
	items = line.strip().split('\t')


	gts = []
	for sample in items[8:]:
		gt=sample.split(":")[0]
		gts.append(gt)

	print(dict(zip(header[8:], gts))

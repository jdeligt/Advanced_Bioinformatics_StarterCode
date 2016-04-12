# AdvancedBioinformatics_2016_1
Repo for teaching purposes in the Advanced Bioinformatics course, UMC Utrecht, 2016


For the Pandas Python examples please see the seperate file in the repository

Example of how to use bio-vcf to extract genotypes:
```bash
bio-vcf --seval 's.gt' < PIK3CA.vcf
```

Example of how to use python to extarct the Genotypes:
```python
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
```

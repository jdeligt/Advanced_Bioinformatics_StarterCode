# AdvancedBioinformatics_2016_1
Repo for teaching purposes in the Advanced Bioinformatics course, UMC Utrecht, 2016


For the Pandas Python examples please see the seperate file in the repository

###Example of how to use bio-vcf to extract genotypes:
```bash
bio-vcf --seval 's.gt' < PIK3CA.vcf
```

###Example of how to use python to extarct the Genotypes:
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

# Selected from favourite examples made by Thomaz and Roel
### A 4 way JOIN with a WHERE
```sql
SELECT Employees.FirstName, SUM(OrderDetails.Quantity*Products.Price) AS CASH FROM OrderDetails
JOIN Orders
ON OrderDetails.OrderID = Orders.OrderID
JOIN Employees
ON Orders.EmployeeID = Employees.EmployeeID
JOIN Products
ON OrderDetails.ProductID = Products.ProductID
GROUP BY Employees.FirstName
```

###CREATE TABLE
```sql
---
--- Table 'Samples'
---
CREATE TABLE samples (
  id              INTEGER PRIMARY KEY AUTO_INCREMENT NOT NULL,
  sample_id       INTEGER NOT NULL,
  --- ... All fields you're interested in ...
);

---
--- Table 'Mutations'
---
CREATE TABLE mutations (
  id               INTEGER PRIMARY KEY AUTO_INCREMENT NOT NULL,
  sample_id        INTEGER NOT NULL,
  chromosome       INTEGER NOT NULL,
  position         INTEGER(12) NOT NULL,
  known_identifier VARCHAR(255) NOT NULL,
  reference_base   VARCHAR(1) NOT NULL,
  alternative_base VARCHAR(1) NOT NULL,
  quality          FLOAT NOT NULL,
  filter           TEXT,
  info             TEXT,
  format           TEXT
);
```

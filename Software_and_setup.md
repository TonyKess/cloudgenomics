# Setup
The first step in working in the cloud is selecting a virtual machine. There are many permutations of choices, but the main things we require for genomics are:
1. A unix environment
2. A flexible software installation environment
3. Large storage space - (e.g. 3TB/Illumina sequencing lane to handle intermediate files in whole genome sequencing projects)
4. Large number of cores (100+) for processes that require parallelization (e.g. adapter trimming), large amount of memory (e.g. 250 GB) for processes such as alignment. 

Currently we are using Azure as our cloud provider. We can view VM options and additions in the (Azure Price Calculator)[https://azure.microsoft.com/en-ca/pricing/calculator/].
# Software
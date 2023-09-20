# Setup
The first step in working in the cloud is selecting a virtual machine. There are many permutations of choices, but the main things we require for genomics are:
1. A unix environment
2. A flexible software installation environment
3. Large storage space - (e.g. 3TB/Illumina sequencing lane to handle intermediate files in whole genome sequencing projects)
4. Large number of cores (100+) for processes that require parallelization (e.g. adapter trimming), large amount of memory (e.g. 250 GB) for processes such as alignment. 

Currently we are using Azure as our cloud provider. We can view VM options and additions in the [Azure Price Calculator.] (https://azure.microsoft.com/en-ca/pricing/calculator/). We have set up two separate virtual machines for these analyses: an HB and M series, both with 120 + cores. We have also partitioned an 8TB and 32 TB storage disk, and a storage account. All of these are provisioned under the same [resource group.](https://learn.microsoft.com/en-us/azure/azure-resource-manager/management/manage-resource-groups-portal). Info on connecting to a selected VM is provided in internal documentation.

# Software
We also require software for interacting with cloud storage, and working in the cloud environment with code that we are developing. We use three programs for this:
[Visual Studio Code](https://code.visualstudio.com/download) for connecting via the terminal and coding.
[Microsoft Azure Storage Explorer](https://azure.microsoft.com/en-us/products/storage/storage-explorer) for managing large stored genomic datasets.
[azcopy](https://learn.microsoft.com/en-us/azure/storage/common/storage-use-azcopy-v10) for moving data between storage and virtual machines.

With these tools in place and a VM selected, we can begin configuring an enviroment for working with genomic data on a virtual machine.

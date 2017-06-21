## Synopsis
gcn.mops is forked from package cn.mops version 1.18.0. In the original package (ie cn.mops), the program executes on CPU while in gcn.mops, portions of the code was re-written and new code was introduced in order to make the core of cn.mops run on nVidia GPU, using CUDA, to speed up processing.

## Motivation (Abstract)
"cn.MOPS is a model-based algorithm used to quantitatively detect copy-number variations in next-generation, DNA-sequencing data. The algorithm is implemented as an R package and can speed up processing with multi-CPU parallelism. However, the maximum achievable speedup is limited by the overhead of multi-CPU parallelism, which increases with the number of CPU cores used. In this thesis, an alternative mechanism of process acceleration is proposed. Using one CPU core and a GPU device, the proposed solution, gcn.MOPS, achieved a speedup factor of 159× and decreased memory usage by more than half. This speedup was substantially higher than the maximum achievable speedup in cn.MOPS, which was ∼20×."

For more information, refer to [my thesis](https://dspace.library.uvic.ca//handle/1828/8286).

## Plaform
gcn.mops was tested on the following platform:
* OS: Ubuntu 16.04 - x64
* Workstation: Dell T7500
* CPU: 2x Hexa-core Intel Xeon X5650
* RAM: ECC 46 GB
* GPU0: nVidia Tesla C2050 (Fermi)
* GPU1: nVidia Quadro 4000 (Fermi)
* GPU Driver: 367.57
* CUDA: 8.0
* NVCC: 8.0, V8.0.44
* GCC: (Ubuntu 5.4.0-6ubuntu1~16.04.4) 5.4.0 20160609
* R: version 3.2.3 (Wooden Christmas-Tree)

## Installation
All original R function interfaces from cn.mops are available in package gcn.mops besides the ones that were newly introduced. This is done intentionally to ease any side-by-side comparison between cn.mops and gcn.mops. calling functions in gcn.mops package should always be preceeded by "gcn.mops::" in case package cn.mops is also installed.

To install, simply execute the script 'recompile_install.sh' from the terminal.

## Test
Installation script will prompt the user to run a test if desired. Aternatively, you may run the test script ('test_script.R') directly after initiating a new R session in terminal.

## Notes
* Modifications to "./gcn.mops/src/Makevars" might be required depending on your platform configurations.
* The package was only tested and optimized for devices of compute capability 2.0 (Fermi).
* The program partitions the input data according to free GPU memory. If other processes are using the GPU such as X11 and web browsers, gcn.mops might fail to allocate memory due to race condition.  gcn.mops will first prob the device for free memory and then it will attempt to reserve all the available memory. In the case such that another process allocates GPU memory after gcnmops probs the device and before gcnmops allocates memory on the device, memory allocation failure might occur. However, informal testing shows that this never occurred. Once the program succeeds in allocating memory, the memory won't be freed until all processing is finalized.
* directory "./R\_scripts" contains two R scripts to test-run cn.mops and gcn.mops for a small dataset in directory "./R\_scripts/data".

## License
LGPL (>= 2.0)

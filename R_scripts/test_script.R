library(gcn.mops)

if(!exists("bamDataRanges"))
	bamDataRanges <- readRDS("data/bamDataRanges21.rds")

message("**********************************************")
message("****** running gcn.mops::gcn.mops (GPU) ******")
message("**********************************************")
gpu_res <- gcn.mops::gcn.mops(bamDataRanges, I = c(0.025, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4),
                        classes = c("CN0", "CN1", "CN2", "CN3", "CN4", "CN5", "CN6", "CN7", "CN8"),
                         priorImpact = 10, cyc = 20, parallel = 0, gpu=0, norm = 1,
                         normType = "poisson", sizeFactor = "mean", normQu = 0.25,
                         quSizeFactor = 0.75, upperThreshold = 0.5, lowerThreshold = -0.8,
                         minWidth = 5, segAlgorithm = "fast", minReadCount = 1,
                         useMedian = FALSE, returnPosterior = FALSE)
 
gpu_res <- gcn.mops::calcIntegerCopyNumbers(gpu_res)
print(gpu_res)

message("*********************************************")
message("****** running gcn.mops::cn.mops (CPU) ******")
message("*********************************************")
cpu_res <- gcn.mops::cn.mops(bamDataRanges, I = c(0.025, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4),
                        classes = c("CN0", "CN1", "CN2", "CN3", "CN4", "CN5", "CN6", "CN7", "CN8"),
                         priorImpact = 10, cyc = 20, parallel = 0, norm = 1,
                         normType = "poisson", sizeFactor = "mean", normQu = 0.25,
                         quSizeFactor = 0.75, upperThreshold = 0.5, lowerThreshold = -0.8,
                         minWidth = 5, segAlgorithm = "fast", minReadCount = 1,
                         useMedian = FALSE, returnPosterior = FALSE)

cpu_res <- gcn.mops::calcIntegerCopyNumbers(cpu_res)
print(cpu_res)

message("***************************************************")
message("****** Execution Time for local modeling (s) ******")
message("***************************************************")

message(sprintf("gcn.mops: %.2f sec", gpu.time.modeling))
message(sprintf("cn.mops:  %.2f sec", cpu.time.modeling))

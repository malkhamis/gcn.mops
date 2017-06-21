library(gcn.mops)
bamDataRanges <- readRDS("data/bamDataRangesMT.rds")

library(future)
maxObjSize <- 2 #GB
options(future.globals.maxSize= maxObjSize*1024^3)

#plan(multicore)
plan(cluster, workers = c("localhost", "localhost"))
experiment1  %<-% {   
  library(gcn.mops)
  res <- gcn.mops::gcn.mops(bamDataRanges, I = c(0.025, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4),
                classes = c("CN0", "CN1", "CN2", "CN3", "CN4", "CN5", "CN6", "CN7", "CN8"),
                priorImpact = 10, cyc = 20, parallel = 0, gpu=0, norm = 0,
                normType = "poisson", sizeFactor = "mean", normQu = 0.25,
                quSizeFactor = 0.75, upperThreshold = 0.5, lowerThreshold = -0.8,
                minWidth = 5, segAlgorithm = "fast", minReadCount = 1,
                useMedian = FALSE, returnPosterior = FALSE)
	res <- gcn.mops::calcIntegerCopyNumbers(res)
}

experiment2 %<-% {   
  library(gcn.mops)
  res <- gcn.mops::gcn.mops(bamDataRanges, I = c(0.025, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4),
                classes = c("CN0", "CN1", "CN2", "CN3", "CN4", "CN5", "CN6", "CN7", "CN8"),
                priorImpact = 10, cyc = 20, parallel = 0, gpu=1, norm = 0,
                normType = "poisson", sizeFactor = "mean", normQu = 0.25,
                quSizeFactor = 0.75, upperThreshold = 0.5, lowerThreshold = -0.8,
                minWidth = 5, segAlgorithm = "fast", minReadCount = 1,
                useMedian = FALSE, returnPosterior = FALSE)
  res <- gcn.mops::calcIntegerCopyNumbers(res)
}


experiment1 # to sync
experiment2 # to sync
message()
message("****** results for GPU 0 ******")
print(experiment1)
message()
message("****** results for GPU 1 ******")
print(experiment2)

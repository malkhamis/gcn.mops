# uncomment relevant code in gcn.mops::cn.mops.R and gcn.mops::gcn.mops.R

num_tolerance = 1E-100
message("*** checking correctness of results ***")
#check lambda
message("Lambda numeric accuracy:")
message(all.equal(cpu_L, gpu_L, tolerance = num_tolerance, 
                  check.attributes = TRUE, use.names = TRUE, 
                  all.names = TRUE, check.names = TRUE ))
message("Lambda Identical? ", identical(cpu_L, gpu_L))

#check alpha
message("Alphs numeric accuracy:")
message(all.equal(cpu_A, gpu_A, tolerance = num_tolerance, 
                  check.attributes = TRUE, use.names = TRUE, 
                  all.names = TRUE, check.names = TRUE ))
message("Alpha Identical? ", identical(cpu_A, gpu_A))

#check expectedCN
message("expectedCN accuracy:")
message(all.equal(cpu_CN, gpu_CN, 
                  check.attributes = TRUE, use.names = TRUE, 
                  all.names = TRUE, check.names = TRUE ))
message("expCN Identical? ", identical(cpu_CN, gpu_CN))

#check sini
message("sini numeric accuracy:")
message(all.equal(cpu_sINI, gpu_sINI, tolerance = num_tolerance,
                  check.attributes = TRUE, use.names = TRUE, 
                  all.names = TRUE, check.names = TRUE ))
message("sini Identical? ", identical(cpu_sINI, gpu_sINI))

#check ini
message("ini numeric accuracy:")
message(all.equal(cpu_INI, gpu_INI, tolerance = num_tolerance,
                  check.attributes = TRUE, use.names = TRUE, 
                  all.names = TRUE, check.names = TRUE ))
message("ini Identical? ", identical(cpu_INI, gpu_INI))

#check post
message("post numeric accuracy:")
message(all.equal(cpu_POST, gpu_POST, tolerance = num_tolerance,
                  check.attributes = TRUE, use.names = TRUE, 
                  all.names = TRUE, check.names = TRUE ))
message("post Identical? ", identical(cpu_POST, gpu_POST))

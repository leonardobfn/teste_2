

setwd("C:/Users/leona/OneDrive/Documentos/GitHub/RThesis")
ALPHA <- c ("alpha35",
            "alpha50",
            "alpha65",
            "alpha80",
            "alpha95")
MODEL = c("Model_5","Model_6")
Method = c("1","2")
N = c(84, 144, 168)
for (model in MODEL) {
  path1 = paste0("Data_simulation/", model)
  path2 = paste0(path1, "/simulations")
  dir.create(path1)
  dir.create(path2)
  for (n in N) {
    path3 = paste0("Data_simulation/", model, "/simulations/n", n)
    dir.create(path3)
    for (alpha in ALPHA) {
      path4 = paste0(path3, "/", alpha)
      dir.create(path4)
    }
  }
}

for (model in MODEL) {
  path2_0 = paste0("Data_simulation/", model)
  path2_1 = paste0(path2_0, "/estimates/")
  dir.create(path2_1)
  for (Met in Method) {
    path2_2 = paste0(path2_1, "Method_", Met)
    path2_3 = paste0(path2_1, "Method_", Met,"/estimates")
    dir.create(path2_2)
    dir.create(path2_3)
    for (n in N) {
      path3 = paste0(path2_3, "/n", n)
      dir.create(path3)
      for (alpha in ALPHA) {
        path4 = paste0(path3, "/", alpha)
        dir.create(path4)
      }
    }
  }
}


for (model in MODEL) {
  path2_0 = paste0("Data_simulation/", model)
  path2_1 = paste0(path2_0, "/estimates/")
  for (Met in Method) {
    path2_2 = paste0(path2_1, "Method_", Met)
    path2_3 = paste0(path2_1, "Method_", Met,"/hessian")
    dir.create(path2_3)
    for (n in N) {
      path3 = paste0(path2_3, "/n", n)
      dir.create(path3)
      for (alpha in ALPHA) {
        path4 = paste0(path3, "/", alpha)
        dir.create(path4)
      }
    }
  }
}

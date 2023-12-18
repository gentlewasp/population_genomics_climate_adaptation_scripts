library(OptM)
library(RColorBrewer)
library(R.utils)
library(R.oo)
library(R.methodsS3)


linear = optM("./sum_m") 
plot_optM(linear, method = "Evanno", plot = F, pdf = "OptM.pdf")



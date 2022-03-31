library(nnet)

model_file = commandArgs(trailingOnly=TRUE)[1]
test_file = commandArgs(trailingOnly=TRUE)[2]
out_file = commandArgs(trailingOnly=TRUE)[3]

model <- readRDS (model_file)

test <- read.table(test_file,header=T,sep="\t")
pred_sv <- predict(model, test, type="probs")

write.table(pred_sv, out_file, quote=F, col.names=F)

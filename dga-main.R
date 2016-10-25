source("dga-prep.R")
library("caret")
suppressPackageStartupMessages(library(caret))
# make this repeatable
set.seed(1492)
# if we pass in a factor, it will do the stratified sampling on it.
# this will return the row numbers to include in the training data
datasample=andrew
trainindex <- createDataPartition(datasample$subclass, p=0.80, list=F)
# only train with these fields:
fields <- c("class", "length", "entropy", "dict", "gram345", "onegram",
            "twogram", "threegram", "fourgram","fivegram")





# Now you can create a training and test data set 
traindga <- datasample[trainindex, fields]
# going to leave all the fields in the test data
testdga <- datasample[-trainindex,  ]




# set up the training control attributes:
ctrl_fast <- trainControl(method="cv", 
                     repeats=1,
                     number=10, 
                     summaryFunction=twoClassSummary,
                     verboseIter=T,
                     classProbs=TRUE,
                     allowParallel = TRUE)
                     
library(doMC)
registerDoMC(cores=4)
ctrl <- trainControl(method="repeatedcv", 
                     repeats=5, 
                     verboseIter =T,
                     summaryFunction=twoClassSummary,
                     classProbs=TRUE,
                     allowParallel = TRUE)

rfFit <- train(class ~ .,
               data = traindga,
               metric="ROC",
               method = "rf",
               trControl = ctrl_fast)




svmFit <- train(class ~ .,
                data = traindga,
                method = "svmRadial",
                preProc = c("center", "scale"),
                metric="ROC",
                tuneLength = 10,
                trControl = ctrl_fast)



registerDoMC(cores=1)
c45Fit <- train(class ~ .,
                data = traindga,
                method = "J48",
                metric="ROC",
                trControl = ctrl_fast)

glmFit <- train(class ~ .,
                data = traindga,
                method = "glm",
                family=binomial(link='logit'),
                metric="ROC",
                trControl = ctrl_fast)




save(glmFit,c45Fit,svmFit,rfFit,file="models/andrew_models.rda",compress = "xz")
predsglm=predict(glmFit,testdga[,-1])
confusionMatrix(predsglm,testdga$class) #first level failure, second level success

predsc45=predict(c45Fit,testdga)
confusionMatrix(predsc45,testdga$class) #first level failure, second level success

predssvm=predict(svmFit,testdga[,-1])
confusionMatrix(predssvm, testdga$class)

predsrf=predict(rfFit,testdga,type='raw')
confusionMatrix(predsrf,testdga$class) #first level failure, second level success

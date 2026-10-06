# ----------------------------------------------------------------------------
# Legacy script from the MSc thesis annex (June 2019).
# Thesis: "Predicción de respuesta a fármacos quimioterapéuticos a partir de
# datos genómicos", MSc in Bioinformatics and Biostatistics (UOC-UB),
# deposited at https://hdl.handle.net/10609/97486
# Original file: ridgelassoelasticdiscretizadosrandomsearch.R (dated 2019-06-04 in the deposited annex),
# original encoding ISO-8859-1 (Latin-1), CRLF line endings.
# Unchanged except for text encoding (converted to UTF-8) and this header.
# Original SHA-256: 0a2c29ffd61072f4ccfe2463ccafa0e51192a52cecc783fa5a5e15a2455c40e1
# ----------------------------------------------------------------------------
#Preprocesamiento de los datos:
##Cargamos las librerias que vamos a usar:
library(caret)
library(ROCR)
library(glmnet)
library(resample)
library(plyr)
library(caTools)

##Leemos el archivo.rds que creamos anteriormente y representa la matriz de expresion:
data <- readRDS("finalmatrix.rds")
quantile(data$AUC)

##Creacion de subsets dependiendo de la Varianza:
##Calculamos las varianzas de cada columna del dataframe y las ordenamos en orden decreciente
##para seleccionar mas tarde los X genes de expresion con varianzas mas altas.
n <- 200
data.var <- sort(colVars(data, na.rm = TRUE), decreasing = TRUE) 
data.var.n <- data.var[1:n]


##Seleccionamos los nombres de los X genes con las varianzas mas altas:
vectornombres <- names(data.var.n)

##Filtramos el dataframe original que hemos llamado "data" con las X columnas pertenecientes
## a los genes que hemos filtrado por sus varianzas.
gen.var.n <- subset(data, select = vectornombres)

##Añadimos la columna de AUC que no habiamos seleccionado previamente pues no corresponde
##a ningun gen.
gen.var.AUC <- data.frame(gen.var.n, data$AUC)
gen.var.AUC = rename(gen.var.AUC, c(data.AUC= "class"))




#Particion 80/20:
##Tomamos un 80% de los datos como conjunto de entrenamiento y el 20% restante
##como conjunto de test:
##Preparamos los valores AUC para indicar cuando son sensibles o resistentes
##a un farmaco las lineas celulares de estudio, es decir, discretizarlos:
tercercuantil <- Quantile(gen.var.AUC$class, 0.75)
gen.var.AUC$class[gen.var.AUC$class<tercercuantil] <- 0
gen.var.AUC$class[gen.var.AUC$class>=tercercuantil] <- 1
gen.var.AUC$class <- factor(gen.var.AUC$class, levels = c("0", "1"), labels= c("Resistente", "Sensible"))
predictors <- names(gen.var.AUC)[!names(gen.var.AUC) %in% "class"]


##Tomamos un 80% de los datos como conjunto de entrenamiento y el 20% restante
##como conjunto de test:
set.seed(1234)
split <- sample.split(gen.var.AUC$class, SplitRatio = 0.80)
data_train <- subset(gen.var.AUC, split==TRUE)
data_test_total <- subset(gen.var.AUC, split==FALSE)
data_test <- data_test_total[,1:n]
data_test_class <- data_test_total$class
x <- model.matrix(class ~., data_train)[,-1]
y <- data_train$class
x_test <- model.matrix(class ~., data = data_test_total) [,-1]
y_test <- data_test_total$class


#Busqueda de hiperparametros y ajuste directo de regresion logistica de Ridge:
fit_ridge <- cv.glmnet(x,
                       y,
                       family = "binomial",
                       alpha = 0, 
                       type.measur='auc',
                       nfolds= 10)
plot(fit_ridge)

##Preccion del modelo ajustado con los datos de prueba o testSet, matriz de confusion y curva ROC:
pred_ridge <- predict(fit_ridge, newx = x_test, 
               s = "lambda.min", type = "class")
pred_ridge <- factor(pred_ridge, levels = c("Resistente", "Sensible"), labels= c("Resistente", "Sensible"))
confusionMatrix(pred_ridge, y_test)
prob_std <- predict(fit_ridge, newx = x_test, type="response", "lambda.min")
pred_ridge <- prediction(prob_std , y_test)
perf_ridge <- performance(pred_ridge, measure = "tpr", x.measure = "fpr")
plot(perf_ridge)

##Guardado del modelo ajustado de Ridge:
save(fit_ridge, file = "fit_randomsearch_Ridge_discretizado_erlotinib_8020.RData")

#Busqueda de hiperparametros y ajuste directo de LASSO:
fit_lasso <- cv.glmnet(x,
                       y,
                       family = "binomial",
                       alpha = 1, 
                       type.measur='auc',
                       nfolds= 10)
plot(fit_lasso)

##Preccion del modelo ajustado con los datos de prueba o testSet, matriz de confusion y curva ROC:
pred_lasso <- predict(fit_lasso, newx = x_test, 
                      s = "lambda.min", type = "class")
pred_lasso <- factor(pred_lasso, levels = c("Resistente", "Sensible"), labels= c("Resistente", "Sensible"))
confusionMatrix(pred_lasso, y_test)
prob_stdlasso <- predict(fit_lasso, newx = x_test, type="response", "lambda.min")
pred_lasso <- prediction(prob_stdlasso , y_test)
perf_lasso <- performance(pred_lasso, measure = "tpr", x.measure = "fpr")
plot(perf_lasso)

##Guardado del modelo ajustado de LASSO:
save(fit_lasso, file = "fit_randomsearch_LASSO_discretizado_erlotinib_8020.RData")

#Busqueda de hiperparametros y ajuste directo de Elastic net:
fit_Elasticnet <- cv.glmnet(x,
                            y,
                            family = "binomial",
                            type.measur='auc',
                            nfolds= 10)
print(fit_Elasticnet)

##Preccion del modelo ajustado con los datos de prueba o testSet, matriz de confusion y curva ROC:
pred_Elasticnet <- predict(fit_Elasticnet, newx = x_test, 
                           s = "lambda.min", type = "class")
pred_Elasticnet <- factor(pred_Elasticnet, levels = c("Resistente", "Sensible"), labels= c("Resistente", "Sensible"))
confusionMatrix(pred_Elasticnet, y_test)
prob_stdElasticnet <- predict(fit_Elasticnet, newx = x_test, type="response", "lambda.min")
pred_Elasticnet <- prediction(prob_stdElasticnet , y_test)
perf_Elasticnet <- performance(pred_Elasticnet, measure = "tpr", x.measure = "fpr")
plot(perf_Elasticnet)

##Guardado del modelo ajustado de Elastic net:
save(fit_Elasticnet, file = "fit_randomsearch_Elasticnet_discretizado_erlotinib_8020.RData")

#60/20/20 (3:1:1 training, validation, testing)#####
##Division de los datos en 80% y 20%:
set.seed(1234)
inTrainingSet602020 <- createDataPartition(gen.var.AUC$class, p = 0.8, list = FALSE)
train602020 <- gen.var.AUC[inTrainingSet602020, ]
test602020 <- gen.var.AUC[-inTrainingSet602020, ]
validation602020 <- sample(1:nrow(train602020),.75*nrow(train602020))
validationdata602020 <- train602020[validation602020,]
x602020 = validationdata602020[,predictors]
y602020 = validationdata602020$class
xtrain602020 = train602020[,predictors]
ytrain602020 =train602020$class
x602020 = as.matrix(x602020)
xtrain602020 = as.matrix(xtrain602020)

##Comprobacion de las particiones:
Reduce("intersect",list(train602020 ,validationdata602020,test602020))



#Busqueda de hiperparametros y ajuste directo de Ridge:
control <- trainControl(method = "repeatedcv", 
                        number = 10, 
                        repeats = 5)
lambda <- 10^seq(-3, 3, length = 100)

ridge_caret602020_validation <- train(x602020, 
                                      y602020, 
                                      method = "glmnet", 
                                      metric = "Accuracy",
                                      tuneGrid =expand.grid(alpha=0, lambda = lambda),
                                      trControl = control)
print(ridge_caret602020_validation)


##ajuste del modelo de regresion Ridge:
fit_ridge602020 <- glmnet(xtrain602020,
                       ytrain602020,
                       family = "binomial",
                       alpha = 0,
                       lambda = 1000)
print(fit_ridge602020)

##Preccion del modelo ajustado con los datos de prueba o testSet, matriz de confusion y curva ROC:
pred_ridge602020 <- predict(fit_ridge602020, newx = as.matrix(test602020[,predictors]), 
                      s = "lambda.min", type = "class")
pred_ridge602020 <- factor(pred_ridge602020, levels = c("Resistente", "Sensible"), labels= c("Resistente", "Sensible"))
confusionMatrix(pred_ridge602020, test602020$class)
prob_std602020 <- predict(fit_ridge602020, newx = as.matrix(test602020[,predictors], type="response", "lambda.min"))
pred_ridge602020 <- prediction(prob_std602020 , test602020$class)
perf_ridge602020 <- performance(pred_ridge602020, measure = "tpr", x.measure = "fpr")
plot(perf_ridge602020)

##Guardado del modelo ajustado de Ridge:
save(fit_ridge602020, file = "fit_randomsearch_Ridge_discretizado_erlotinib_602020.RData")

#Busqueda de hiperparametros LASSO:
lasso_caret602020_validation <- train(x602020, 
                                      y602020, 
                                      method = "glmnet", 
                                      metric = "Accuracy",
                                      tuneGrid =expand.grid(alpha=1, lambda = lambda),
                                      trControl = control)
print(lasso_caret602020_validation)

##ajuste del modelo de regresion LASSO:
fit_lasso602020 <- glmnet(xtrain602020,
                             ytrain602020,
                             family = "binomial",
                             alpha = 1,
                          lambda = 0.04977024)
print(fit_lasso602020)

##Preccion del modelo ajustado con los datos de prueba o testSet, matriz de confusion y curva ROC:
pred_lasso602020 <- predict(fit_lasso602020, newx = as.matrix(test602020[,predictors]), 
                            s = "lambda.min", type = "class")
pred_lasso602020 <- factor(pred_lasso602020, levels = c("Resistente", "Sensible"), labels= c("Resistente", "Sensible"))
confusionMatrix(pred_lasso602020, test602020$class)
prob_stdlasso602020 <- predict(fit_lasso602020, newx = as.matrix(test602020[,predictors], type="response", "lambda.min"))
pred_lasso602020 <- prediction(prob_stdlasso602020 , test602020$class)
perf_lasso602020 <- performance(pred_lasso602020, measure = "tpr", x.measure = "fpr")
plot(perf_lasso602020)

save(fit_lasso602020, file = "fit_randomsearch_LASSO_discretizado_erlotinib_602020.RData")

#Busqueda de hiperparametros Elastic Net:
Elasticnet_caret602020_validation <- train(x602020, 
                                      y602020, 
                                      method = "glmnet", 
                                      metric = "Accuracy",
                                      tuneLength = 20,
                                      trControl = control)
print(Elasticnet_caret602020_validation)

##ajuste del modelo de regresion Elastic net:
fit_Elasticnet602020 <- glmnet(xtrain602020,
                             ytrain602020,
                             family = "binomial",
                             alpha =0.1473684,
                             lambda =0.2952464)
print(fit_Elasticnet602020)

##Preccion del modelo ajustado con los datos de prueba o testSet, matriz de confusion y curva ROC:
pred_Elasticnet602020 <- predict(fit_Elasticnet602020, newx = as.matrix(test602020[,predictors]), 
                            s = "lambda.min", type = "class")
pred_Elasticnet602020 <- factor(pred_Elasticnet602020, levels = c("Resistente", "Sensible"), labels= c("Resistente", "Sensible"))
confusionMatrix(pred_Elasticnet602020, test602020$class)
prob_stdElasticnet602020 <- predict(fit_Elasticnet602020, newx = as.matrix(test602020[,predictors], type="response", "lambda.min"))
pred_Elasticnet602020 <- prediction(prob_stdElasticnet602020 , test602020$class)
perf_Elasticnet602020 <- performance(pred_Elasticnet602020, measure = "tpr", x.measure = "fpr")
plot(perf_Elasticnet602020)

##Guardado del modelo ajustado de Elastic net:
save(fit_Elasticnet602020, file = "fit_randomsearch_Elasticnet_discretizado_erlotinib_602020.RData")
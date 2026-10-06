# ----------------------------------------------------------------------------
# Legacy script from the MSc thesis annex (June 2019).
# Thesis: "Predicción de respuesta a fármacos quimioterapéuticos a partir de
# datos genómicos", MSc in Bioinformatics and Biostatistics (UOC-UB),
# deposited at https://hdl.handle.net/10609/97486
# Original file: ridgelassoelasticcontinuosrandomsearch.R (dated 2019-06-06 in the deposited annex),
# original encoding ISO-8859-1 (Latin-1), CRLF line endings.
# Unchanged except for text encoding (converted to UTF-8) and this header.
# Original SHA-256: db20b30e843c1346b5fb3063b96992da75c72275550a9c4b41f87b882ebb30ae
# ----------------------------------------------------------------------------
#Preprocesamiento de los datos:
##Cargamos las librerias que vamos a usar:
library(caret)
library(miscTools)
library(glmnet)
library(resample)
library(plyr)
library(caTools)
library(ggplot2)

##Leemos el archivo.rds que creamos anteriormente y representa la matriz de expresion:
data <- readRDS("finalmatrix.rds")

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

#Busqueda de hiperparametros y ajuste directo de regresion de Ridge:
set.seed(123)
fit_ridge <- cv.glmnet(x,
                       y,
                       alpha = 0, 
                       type.measur='mse',
                       nfolds= 10)
print(fit_ridge)
##Preccion del modelo ajustado con los datos de prueba o testSet:
prediction_ridge <- predict(fit_ridge,x_test)
prediction_ridge <- as.numeric(prediction_ridge)

##Calculo de la R^2:
rsq <- function (x, y) cor(x, y) ^ 2
r2_ridge <- rsq(y_test, y_test-prediction_ridge)

##Creacion de la gráfica R^2:
p <- ggplot(aes(x=actual, y=pred),
            data=data.frame(actual=y_test, pred=prediction_ridge))
p + geom_point() +
  geom_abline(color="red") +
  ggtitle(paste("Ridge Regression in R r^2=", r2_ridge, sep=""))

##Guardado del modelo ajustado de Ridge:
save(fit_ridge, file = "fit_randomsearch_Ridge_continuo_erlotinib_8020.RData")


#Busqueda de hiperparametros y ajuste directo de LASSO:
fit_Lasso <- cv.glmnet(x,
                       y,
                       alpha = 1, 
                       type.measur='mse',
                       nfolds= 10)
print(fit_Lasso)
##Preccion del modelo ajustado con los datos de prueba o testSet:
prediction_Lasso <- predict(fit_Lasso,x_test)
prediction_Lasso <- as.numeric(prediction_Lasso)

##Calculo de la R^2:

r2_Lasso <- rsq(y_test, y_test-prediction_Lasso)
##Creacion de la gráfica R^2:
p <- ggplot(aes(x=actual, y=pred),
            data=data.frame(actual=y_test, pred=prediction_Lasso))
p + geom_point() +
  geom_abline(color="red") +
  ggtitle(paste("LASSO Regression in R r^2=", r2_Lasso, sep=""))

##Guardado del modelo ajustado de LASSO:
save(fit_Lasso, file = "fit_randomsearch_LASSO_continuo_erlotinib_8020.RData")


#Busqueda de hiperparametros y ajuste directo de Elastic net:
fit_Elasticnet <- cv.glmnet(y =y,
                    x = x,
                    type.measur='mse',
                    nfolds= 10)
print(fit_Elasticnet)
##Preccion del modelo ajustado con los datos de prueba o testSet:
prediction_Elasticnet <- predict(fit_Elasticnet,x_test)
prediction_Elasticnet <- as.numeric(prediction_Elasticnet)

##Calculo de la R^2:
r2_Elasticnet<- rsq(y_test, y_test-prediction_Elasticnet)

##Creacion de la gráfica R^2:
p <- ggplot(aes(x=actual, y=pred),
            data=data.frame(actual=y_test, pred=prediction_Elasticnet))
p + geom_point() +
  geom_abline(color="red") +
  ggtitle(paste("Elastic net Regression in R r^2=", r2_Elasticnet, sep=""))

##Guardado del modelo ajustado de Elastic net:
save(fit_Elasticnet, file = "fit_randomsearch_Elasticnet_continuo_erlotinib_8020.RData")

#60/20/20 (3:1:1 training, validation, testing)#####
##Division de los datos en 80% y 20%:
samp <- sample(1:nrow(gen.var.AUC),.8*nrow(gen.var.AUC))
remain80 <- gen.var.AUC[samp,] ##  80% 
testing20 <- gen.var.AUC[-samp,]  ## 20%
class_testing20 <- testing20[,"class"]


##Division del 80% en train y validation (3:1)
samp2 <- sample(1:nrow(remain80),.75*nrow(remain80))
train60 <- remain80[samp2,] ## 60%
validation20 <- remain80[-samp2,] ## 20%

##Creacion del test set para predicciones:
x <- model.matrix(class ~., remain80)[,-1]
y <- remain80$class
x_test <- model.matrix(class ~., data = testing20) [,-1]
y_test <- testing20$class


##Comprobacion de las particiones:
Reduce("intersect",list(train60,validation20,testing20))



#Busqueda de hiperparametros Ridge:
control <- trainControl(method = "repeatedcv", 
                        number = 10, 
                        repeats = 5)
lambda <- 10^seq(-3, 3, length = 100)


set.seed(123)
ridge_caret602020_validation <- train(class ~., 
                                    data = validation20, 
                                    method = "glmnet", 
                                    metric = "RMSE",
                                    tuneGrid =expand.grid(alpha=0, lambda = lambda),
                                    trControl = control)
print(ridge_caret602020_validation)

##ajuste del modelo de regresion Ridge:
fit_ridge602020 <- glmnet(y =y,
                    x = x, 
                    alpha = 0,
                    lambda = 4.328761)

##Preccion del modelo ajustado con los datos de prueba o testSet:
predictionridge602020 <-predict(fit_ridge602020, x_test)
predictionridge602020 <- as.numeric(predictionridge602020)

##Calculo de la R^2:
rsq <- function (x, y) cor(x, y) ^ 2
r2_ridge602020 <- rsq(y_test, y_test - predictionridge602020)


##Creacion de la gráfica R^2:
p602020 <- ggplot(aes(x=actual, y=pred),
                  data=data.frame(actual=y_test, pred=predictionridge602020))
p602020 + geom_point() +
  geom_abline(color="red") +
  ggtitle(paste("Ridge Regression in R r^2=", r2_ridge602020, sep=""))

##Guardado del modelo ajustado de Ridge:
save(fit_ridge602020, file = "fit_randomsearch_ridge_continuo_erlotinib_602020.RData")


#Busqueda de hiperparametros LASSO:
Lasso_caret602020_validation <- train(class ~., 
                                      data = validation20, 
                                      method = "glmnet", 
                                      metric = "RMSE",
                                      tuneGrid =expand.grid(alpha=1, lambda = lambda),
                                      trControl = control)
print(Lasso_caret602020_validation)

##ajuste del modelo de regresion LASSO:
fit_Lasso602020 <- glmnet(y =y,
                          x = x, 
                          alpha = 1,
                          lambda = 0.02477076)

##Preccion del modelo ajustado con los datos de prueba o testSet:
predictionLasso602020 <-predict(fit_Lasso602020, x_test)
predictionLasso602020 <- as.numeric(predictionLasso602020)

##Calculo de la R^2:
r2_Lasso602020 <-  rsq(y_test, y_test - predictionLasso602020)

##Creacion de la gráfica R^2:
pLasso602020 <- ggplot(aes(x=actual, y=pred),
                  data=data.frame(actual=y_test, pred=predictionLasso602020))
pLasso602020 + geom_point() +
  geom_abline(color="red") +
  ggtitle(paste("Lasso Regression in R r^2=", r2_Lasso602020, sep=""))

##Guardado del modelo ajustado de LASSO:
save(fit_Lasso602020, file = "fit_randomsearch_LASSO_continuo_erlotinib_602020.RData")

#Busqueda de hiperparametros Elastic Net:
set.seed(123)
Elasticnet_caret602020_validation <- train(class ~., 
                                      data = validation20, 
                                      method = "glmnet", 
                                      metric = "RMSE",
                                      tuneLength= 20,
                                      trControl = control)
print(Elasticnet_caret602020_validation)

##ajuste del modelo de regresion Elastic net:
fit_Elasticnet602020 <- glmnet(y =y,
                          x = x, 
                          alpha = 0.8105263,
                          lambda = 0.07328962)

##Preccion del modelo ajustado con los datos de prueba o testSet:
predictionElasticnet602020 <-predict(fit_Elasticnet602020, x_test)
predictionElasticnet602020 <- as.numeric(predictionElasticnet602020)

##Calculo de la R^2:
r2_Elasticnet602020 <-  rsq(y_test, y_test - predictionElasticnet602020)

##Creacion de la gráfica R^2:
pElasticnet602020 <- ggplot(aes(x=actual, y=pred),
                       data=data.frame(actual=y_test, pred=predictionElasticnet602020))
pElasticnet602020 + geom_point() +
  geom_abline(color="red") +
  ggtitle(paste("Elastic net Regression in R r^2=", r2_Elasticnet602020, sep=""))

##Guardado del modelo ajustado de Elastic net:
save(fit_Elasticnet602020, file = "fit_randomsearch_Elasticnet_continuo_erlotinib_602020.RData")
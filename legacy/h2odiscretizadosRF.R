# ----------------------------------------------------------------------------
# Legacy script from the MSc thesis annex (June 2019).
# Thesis: "Predicción de respuesta a fármacos quimioterapéuticos a partir de
# datos genómicos", MSc in Bioinformatics and Biostatistics (UOC-UB),
# deposited at https://hdl.handle.net/10609/97486
# Original file: h2odiscretizadosRF.R (dated 2019-06-04 in the deposited annex),
# original encoding ISO-8859-1 (Latin-1), CRLF line endings.
# Unchanged except for text encoding (converted to UTF-8) and this header.
# Original SHA-256: 84996636b473d83e1924e258725d646442db7eb1abc38f29543679c474bb40f6
# ----------------------------------------------------------------------------
#Preprocesamiento de los datos:
##Cargamos las librerias que vamos a usar:
library(randomForest)
library(caret)
library(h2o)
library(ROCR)
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

##Tomamos un 80% de los datos como conjunto de entrenamiento y el 20% restante
##como conjunto de test:
set.seed(1234)
split <- sample.split(gen.var.AUC$class, SplitRatio = 0.80)
data_train <- subset(gen.var.AUC, split==TRUE)
data_test_total <- subset(gen.var.AUC, split==FALSE)
data_test <- data_test_total[,1:n]
data_test_class <- data_test_total$class

##Para confirmar que la proporcion de resistentes y sensibles es aproximadamente
##la misma en ambos conjuntos de datos, vemos la distribucion de los tipos en cada set:
table(data_train$class)
table(data_test_class)


##Busqueda de hiperparámetros con H2o:
##Iniciamos el software:
h2o.init(max_mem_size = "5g")
y <- "class"
x <- setdiff(names(data_train), y)
train.h2o <- as.h2o(data_train)

##Creacion de lista con intervalos de hiperparametros para buscar:
hyper_grid.h2o <- list(
  ntrees      = seq(200, 500, by = 150),
  mtries      = seq(15, 35, by = 10),
  max_depth   = seq(20, 40, by = 5),
  min_rows    = seq(1, 5, by = 2),
  nbins       = seq(10, 30, by = 5),
  sample_rate = c(.55, .632, .75)
)

##Criterios de seleccion de hiperparametros:
search_criteria <- list(
  strategy = "RandomDiscrete",
  stopping_metric = "mse",
  stopping_tolerance = 0.005,
  stopping_rounds = 10,
  max_runtime_secs = 30*60
)

##Construccion del algoritmo de random Forest sobre el 
##que se hara la busqueda de hiperparámetros:
random_grid <- h2o.grid(
  algorithm = "randomForest",
  grid_id = "rf_grid2",
  x = x, 
  y = y, 
  training_frame = train.h2o,
  hyper_params = hyper_grid.h2o,
  search_criteria = search_criteria
)

##Guardado de los mejores hiperparametros encontrados en H2o:
grid_perf2 <- h2o.getGrid(
  grid_id = "rf_grid2", 
  sort_by = "Accuracy", 
  decreasing = FALSE
)
print(grid_perf2)

##Ajuste del modelo de Random Forest:
fit_rf10 <- randomForest(class~.,
                        data = data_train,
                        ntree = 500,
                        mtry= 15,
                        max_depth = 35,
                        min_row = 5,
                        nbin = 10,
                        sample_rate =0.632)

prediction <-predict(fit_rf10, data_test)
confusionMatrix(prediction, data_test_class)

##Guardado del modelo ajustado de random Forest:
save(fit_rf10, file = "fit_rf_erlotinib_8020_discretizado_h2o.RData")

##Prediccion con el modelo ajustado de los 
##datos reservados como testSet y cracion de la matriz de confusion y curva Roc.
probs <- predict(fit_rf10, data_test, type= "prob")
pred <- prediction(probs[,2] , data_test_class)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)


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

##Comprobacion de las particiones:
Reduce("intersect",list(train60,validation20,testing20))

##Busqueda de hiperparámetros con H2o:
##Iniciamos el software:
h2o.init(max_mem_size = "5g")
y <- "class"
x <- setdiff(names(validation20), y)
train.h2o <- as.h2o(validation20)

##Creacion de lista con intervalos de hiperparametros para buscar:
hyper_grid.h2o <- list(
  ntrees      = seq(200, 500, by = 150),
  mtries      = seq(15, 35, by = 10),
  max_depth   = seq(20, 40, by = 5),
  min_rows    = seq(1, 5, by = 2),
  nbins       = seq(10, 30, by = 5),
  sample_rate = c(.55, .632, .75)
)

##Criterios de seleccion de hiperparametros:
search_criteria <- list(
  strategy = "RandomDiscrete",
  stopping_metric = "mse",
  stopping_tolerance = 0.005,
  stopping_rounds = 10,
  max_runtime_secs = 30*60
)

##Construccion del algoritmo de random Forest sobre el 
##que se hara la busqueda de hiperparámetros:
random_grid <- h2o.grid(
  algorithm = "randomForest",
  grid_id = "rf_grid2",
  x = x, 
  y = y, 
  training_frame = train.h2o,
  hyper_params = hyper_grid.h2o,
  search_criteria = search_criteria
)

##Guardado de los mejores hiperparametros encontrados en H2o:
grid_perf2 <- h2o.getGrid(
  grid_id = "rf_grid2", 
  sort_by = "Accuracy", 
  decreasing = FALSE
)
print(grid_perf2)

##Ajuste del modelo de Random Forest con los datos de entrenamiento:
fit_rf602020 <- randomForest(class~.,
                        data = remain80,
                        ntree = 500,
                        mtry= 15,
                        max_depth = 35,
                        min_row = 5,
                        nbin = 10,
                        sample_rate =0.632)
##Prediccion con el modelo ajustado de los 
##datos reservados como testSet y cracion de la matriz de confusion y curva Roc.
prediction602020 <-predict(fit_rf602020, testing20)
confusionMatrix(prediction602020, class_testing20)
probs602020 <- predict(fit_rf602020, testing20, type= "prob")
pred602020 <- prediction(probs602020[,2] , class_testing20)
perf602020 <- performance(pred602020, measure = "tpr", x.measure = "fpr")
plot(perf602020)

##Guardado del modelo ajustado de random Forest:
save(fit_rf602020, file = "fit_rf_erlotinib_602020_discretizado_h2o.RData")



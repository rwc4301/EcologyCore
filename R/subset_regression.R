#ACHTUNG: From the RStudio menu, click on "Session" and then "Set Working Directory" to "To Source File Location"
#Script for subset regression of dependent variable against environmental data
#v1.2 Multiple files/model

library(leaps)
library(xtable)
library(jtools)
library(sjPlot)
library(dplyr)
library(tidyverse)
library(caret)


#PARAMETERS ###########################
meta_table<-read.csv("../../data/meta_table.csv",header=T,row.names=1,check.names=FALSE)
regression_method="forward" #exhaustive, backward, forward, seqrep
really_big=FALSE #TRUE/FALSE

#Use selected_explanatory_variables[!selected_explanatory_variables %in% colnames(meta_table)] to debug
#Make sure there are no hyphens "-" in column names, remove them or convert them to underscores "_"

label="Shannon_Hypothesis1"
dependent_table<-read.csv("../alpha_diversity/Diversity_Otus_Hypothesis1.csv",header=T,row.names=1)
dependent_variable<-c("Shannon")


selected_explanatory_variables<-c(
  "pH",
  "Temp",
  "conductivity",
  "turbidity",
  "D50_grain_size",
  "Mean_Gradient_of_Catchment"
)



#/PARAMETERS ###########################

#Ensure that we are only selecting samples for which the meta_table[,selected_explanatory_variables] is complete and
#the samples also exist in the dependent_table
meta_table[,selected_explanatory_variables] <- lapply(meta_table[,selected_explanatory_variables], function(x) as.numeric(as.character(x)))
meta_table<-meta_table[complete.cases(meta_table[,selected_explanatory_variables]),]
meta_table<-meta_table[rownames(meta_table) %in% rownames(dependent_table),]
dependent_table<-dependent_table[rownames(meta_table),,drop=F]
dependent_table<-dependent_table[complete.cases(dependent_table[,dependent_variable]),,drop=F]
meta_table<-meta_table[rownames(dependent_table),]


#HELPER FUNCTION#######################

#Reference: https://stackoverflow.com/questions/19340277/converting-r-formula-format-to-mathematical-equation
# define a function to take a linear regression
#  (anything that supports coef() and terms() should work)
expr.from.lm <- function (fit) {
  # the terms we're interested in
  con <- names(coef(fit))
  # current expression (built from the inside out)
  expr <- quote(epsilon)
  # prepend expressions, working from the last symbol backwards
  for (i in length(con):1) {
    if (con[[i]] == '(Intercept)')
      expr <- bquote(beta[.(i-1)] + .(expr))
    else
      expr <- bquote(beta[.(i-1)] * .(as.symbol(con[[i]])) + .(expr))
  }
  # add in response
  expr <- bquote(.(terms(fit)[[2]]) == .(expr))
  # convert to expression (for easy plotting)
  as.expression(expr)
}


lm.dat<-data.frame(dependent_table[,dependent_variable,drop=F],meta_table[,selected_explanatory_variables,drop=F])

library(dplyr)
library(tidyverse)
library(caret)
library(leaps)

models<-regsubsets(as.formula(paste(dependent_variable,"~",paste(selected_explanatory_variables,collapse=" + "))),
                   data=lm.dat,
                   nvmax=length(selected_explanatory_variables),really.big=really_big,method=regression_method)

nvmax=length(summary(models)[[2]])

res.sum <- summary(models)
data.frame(
  Adj.R2 = which.max(res.sum$adjr2),
  CP = which.min(res.sum$cp),
  BIC = which.min(res.sum$bic)
)


#Reference http://www.sthda.com/english/articles/37-model-selection-essentials-in-r/155-best-subsets-regression-essentials-in-r/
# id: model id
# object: regsubsets object
# data: data used to fit regsubsets
# outcome: outcome variable
get_model_formula <- function(id, object, outcome){
  # get models data
  models <- summary(object)$which[id,-1]
  # Get outcome variable
  form <- as.formula(object$call[[2]])
  outcome <- all.vars(form)[1]
  # Get model predictors
  predictors <- names(which(models == TRUE))
  predictors <- paste(predictors, collapse = "+")
  # Build model formula
  as.formula(paste0(outcome, "~", predictors))
}

get_cv_error <- function(model.formula, data){
  set.seed(1)
  train.control <- trainControl(method = "cv", number = 5)
  cv <- train(model.formula, data = data, method = "lm",
              trControl = train.control)
  cv$results$RMSE
}


# Compute cross-validation error
model.ids <- 1:nvmax
cv.errors <-  purrr::map(model.ids, get_model_formula, models, dependent_variable) %>%
  map(get_cv_error, data = lm.dat) %>%
  unlist()

best_variable_model<-which.min(cv.errors)

#Generate a CV table
CV_table<-data.frame(as.character(map(model.ids, get_model_formula, models, dependent_variable)))
names(CV_table)<-c("Model")
CV_table$`Cross-validation Errors`<-cv.errors
CV_table<-CV_table[order(CV_table$`Cross-validation Errors`,decreasing=FALSE),]
print(xtable(CV_table,display=c("s","s","f"),digits=5), type="html", file=paste(label,"_CV_errors.html",sep=""),html.table.attributes = "border = '1', align = 'center', cellspacing='0', cellpadding='0'")

best_model<-lm(get_model_formula(best_variable_model,models,dependent_variable),data=lm.dat)

#Get visualisations for all best models
for (i in 1:nvmax){
  tmp<-get_model_formula(i,models,dependent_variable)
  current_model<-lm(get_model_formula(i,models,dependent_variable),data=lm.dat)
  if(length(current_model$coefficients)>sum(complete.cases(current_model$coefficients))){
    current_model<-lm(as.formula(paste(dependent_variable,"~",paste(names(current_model$coefficients)[complete.cases(current_model$coefficients)][-1],collapse="+"))),data=lm.dat)
  }
  
  
  #Reference: https://cran.r-project.org/web/packages/sjPlot/vignettes/tab_model_estimates.html
  #Reference: https://www.r-bloggers.com/beautiful-tables-for-linear-model-summaries-rstats/
  #if it fails, then use the other one depending on the version of tab_model
  q<-tab_model(current_model, p.style="scientific_stars", digits=5,show.se = TRUE, show.std = TRUE, show.df=TRUE, show.stat = TRUE,file=NULL)
  #q<-tab_model(current_model, p.style="both", digits=5,show.se = TRUE, show.std = TRUE, show.df=TRUE, show.stat = TRUE,file=NULL)
  p<-as.data.frame(gsub("\n","",q$knitr))
  colnames(p)<-c(" ")
  write.csv(p,file=paste(label,"_M",i,".html",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE)
  
  #Reference: https://cran.r-project.org/web/packages/jtools/vignettes/summ.html
  pdf(paste(label,"_M",i,".pdf",sep=""),width=8,height=4)
  p<-plot_summs(current_model,
                model.names=c(paste("M",i,sep="")),
                plot.distributions=TRUE,rescale.distributions=TRUE,
                omit.coefs=NULL,
                color.class="Rainbow")
  print(p)
  dev.off()
}
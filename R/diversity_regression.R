#ACHTUNG: From the RStudio menu, click on "Session" and then "Set Working Directory" to "To Source File Location"
#Script for subset regression of dependent variable against environmental data
#v1.2 Multiple files/model

# library(leaps)
# library(xtable)
# library(jtools)
# library(sjPlot)
# library(dplyr)
# library(tidyverse)
# library(caret)
# library(dplyr)
# library(tidyverse)
# library(caret)
# library(leaps)

diversity_regression <- function(
  dependent_table,
  meta_table,
  dependent_variable,
  explanatory_variables,
  regression_method = "forward", #exhaustive, backward, forward, seqrep
  really_big = FALSE, #TRUE/FALSE
  num_top_models = 5
) {
  #Use explanatory_variables[!explanatory_variables %in% colnames(meta_table)] to debug
  #Make sure there are no hyphens "-" in column names, remove them or convert them to underscores "_"

  label="Shannon_Hypothesis1"

  #/PARAMETERS ###########################

  #Ensure that we are only selecting samples for which the meta_table[,explanatory_variables] is complete and
  #the samples also exist in the dependent_table
  meta_table[,explanatory_variables] <- lapply(meta_table[,explanatory_variables], function(x) as.numeric(as.character(x)))
  meta_table<-meta_table[complete.cases(meta_table[,explanatory_variables]),]
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


  lm.dat<-data.frame(dependent_table[,dependent_variable,drop=F],meta_table[,explanatory_variables,drop=F])

  models<-leaps::regsubsets(as.formula(paste(dependent_variable,"~",paste(explanatory_variables,collapse=" + "))),
                     data=lm.dat,
                     nvmax=length(explanatory_variables),really.big=really_big,method=regression_method)

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
    purrr::map(get_cv_error, data = lm.dat) %>%
    unlist()

  best_variable_model<-which.min(cv.errors)

  #Generate a CV table
  CV_table<-data.frame(as.character(purrr::map(model.ids, get_model_formula, models, dependent_variable)))
  names(CV_table)<-c("Model")
  CV_table$`Cross-validation Errors`<-cv.errors
  CV_table<-CV_table[order(CV_table$`Cross-validation Errors`,decreasing=FALSE),]
  #print(xtable(CV_table,display=c("s","s","f"),digits=5), type="html", file=paste(label,"_CV_errors.html",sep=""),html.table.attributes = "border = '1', align = 'center', cellspacing='0', cellpadding='0'")

  best_model<-lm(get_model_formula(best_variable_model,models,dependent_variable),data=lm.dat)

  model_formulas <- lapply(1:nvmax, function(x) get_model_formula(x, models, response_variable))
  model_details <- lapply(model_formulas, function(x) {
    m <- lm(x, data = lm.dat)
    if (length(m$coefficients) > sum(complete.cases(m$coefficients))) {
      m <- lm(as.formula(paste(response_variable,"~",paste(names(m$coefficients)[complete.cases(m$coefficients)][-1],collapse="+"))),data=lm.dat)
    }

    return(m)
  })

  model_summaries <- lapply(model_details, function(x) {
      # Get standardized coefficients
    standardized_model <- lm(scale(model.frame(x)[,1]) ~ scale(model.frame(x)[,-1]))
    std_coef <- coef(standardized_model)
    std_se <- summary(standardized_model)$coefficients[,2]

    # Create summary data frame
    summary_df <- data.frame(
      # Term = names(coef(x)),
      Estimate = coef(x),
      Std_Error = summary(x)$coefficients[,2],
      Std_Beta = c(NA, std_coef[-1]),  # NA for intercept
      Std_Beta_SE = c(NA, std_se[-1]), # NA for intercept
      t_value = summary(x)$coefficients[,3],
      p_value = summary(x)$coefficients[,4]
    )

    # Format p-values with stars
    summary_df$p_value <-
      ifelse(summary_df$p_value < 0.001, paste(summary_df$p_value, "***"),
      ifelse(summary_df$p_value < 0.01, paste(summary_df$p_value, "**"),
      ifelse(summary_df$p_value < 0.05, paste(summary_df$p_value, "*"), "")))

    return(summary_df)
  })

  return (structure(list(
    CV_table = CV_table,
    best_model = best_model,
    models = models,
    model_summaries = model_summaries
  ), class = "ECSubsetRegression"))
}

get_cv_table <- function(value) {

}

get_summary_table <- function(current_model) {
  #Reference: https://cran.r-project.org/web/packages/sjPlot/vignettes/tab_model_estimates.html
  #Reference: https://www.r-bloggers.com/beautiful-tables-for-linear-model-summaries-rstats/
  #if it fails, then use the other one depending on the version of tab_model
  q<-tab_model(current_model, p.style="scientific_stars", digits=5,show.se = TRUE, show.std = TRUE, show.df=TRUE, show.stat = TRUE,file=NULL)
  #q<-tab_model(current_model, p.style="both", digits=5,show.se = TRUE, show.std = TRUE, show.df=TRUE, show.stat = TRUE,file=NULL)
  p<-as.data.frame(gsub("\n","",q$knitr))
  colnames(p)<-c(" ")
  #write.csv(p,file=paste(label,"_M",i,".html",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE)

  return (p)
}

summary.ECSubsetRegression <- function(value) {
  #print(xtable(CV_table,display=c("s","s","f"),digits=5), type="html", file=paste(label,"_CV_errors.html",sep=""),html.table.attributes = "border = '1', align = 'center', cellspacing='0', cellpadding='0'")
}

plot_model_summaries <- function(value) {
  #Reference: https://cran.r-project.org/web/packages/jtools/vignettes/summ.html
  #pdf(paste(label,"_M",i,".pdf",sep=""),width=8,height=4)
  p <- jtools::plot_summs(
    value$model_summaries,
    model.names=c(paste("M",i,sep="")),
    plot.distributions=TRUE,
    rescale.distributions=TRUE,
    omit.coefs=NULL,
    color.class="Rainbow"
  )

  return (p)
  #print(p)
  #dev.off()
}

plot_all <- function() {
  all_models_coef <- all_models_coef[all_models_coef$Term != "(Intercept)", ]

  # Custom transformation function for square root including negative values
  custom_sqrt_transform <- function(x, c = 1) {
    sign(x) * sqrt(abs(x) + c)
  }

  custom_cubed_root_transform <- function(x) {
    sign(x) * abs(x)^(1/3)  # Cube root transformation
  }

  # Create a transformation object
  sqrt_transform <- scales::trans_new(
    name = "custom_sqrt",
    transform = custom_sqrt_transform,
    inverse = function(x) {
      sign(x) * (x^2 - 1)  # Inverse transformation
    }
  )

  cubed_root_transform <- scales::trans_new(
    name = "custom_cubed_root",
    transform = custom_cubed_root_transform,
    inverse = function(x) {
      sign(x) * (x^3)  # Inverse transformation
    }
  )

  # Apply the transformation to the heatmap
  p <- ggplot(all_models_coef, aes(x = factor(Model), y = Term, fill = Estimate)) +
    geom_tile() +
    geom_text(aes(label = Stars), color = "black", size = 3) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", 
                        midpoint = 0,
                        trans = sqrt_transform,  # Use the custom transformation
                        name = "Transformed\nCoefficient") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0),
          axis.text.y = element_text(size = 8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(x = "Model Number", y = "Variable",
        title = "Square Root Transformed Coefficients Across Models",
        subtitle = "* p<0.05, ** p<0.01, *** p<0.001")

  return(p)
}
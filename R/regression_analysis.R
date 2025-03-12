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
  df,
  response_variable,
  explanatory_variables,
  regression_method = "forward", #exhaustive, backward, forward, seqrep
  really.big = FALSE #TRUE/FALSE
) {
  # Ensure that we are only selecting samples for which the meta_table[,explanatory_variables] is complete
  # TODO: add interpolation options?
  df[, explanatory_variables] <- lapply(df[,explanatory_variables], function(x) as.numeric(as.character(x)))
  df <- df[complete.cases(df[, explanatory_variables]),]
  df <- df[complete.cases(df[, response_variable]), , drop = FALSE]

  lm.dat <- data.frame(df[, response_variable, drop = FALSE], df[, explanatory_variables, drop = FALSE])

  f <- as.formula(paste(response_variable, "~", paste(explanatory_variables, collapse = " + ")))
  nvmax <- length(explanatory_variables)

  return(leaps::regsubsets(f, data = lm.dat, nvmax = nvmax, really.big = really.big, method = regression_method))
}

best_model <- function(models, response_variable, lm.dat) {
  nvmax <- length(summary(models)[[2]])

  # Compute cross-validation error
  model.ids <- 1:nvmax
  cv.errors <- purrr::map(model.ids, get_model_formula, models, response_variable) %>%
    purrr::map(get_cv_error, data = lm.dat) %>%
    unlist()

  return(lm(get_model_formula(which.min(cv.errors), models, response_variable), data = lm.dat))
}

# Compute cross-validation error
cross_validate <- function(models, response_variable, lm.dat, count = 5) {
  nvmax <- min(count, length(summary(models)[[2]]))

  cv_table <- purrr::map(1:nvmax, get_model_formula, models, response_variable)

  cv_errors <- cv_table %>% purrr::map(get_cv_error, data = lm.dat) %>% unlist()

  # res.sum <- summary(models)
  # data.frame(
  #   Adj.R2 = which.max(res.sum$adjr2),
  #   CP = which.min(res.sum$cp),
  #   BIC = which.min(res.sum$bic)
  # )

  cv_table <- data.frame(as.character(cv_table))
  names(cv_table) <- c("Model")
  cv_table$`Cross-validation Errors` <- cv_errors
  cv_table <- cv_table[order(cv_table$`Cross-validation Errors`, decreasing=FALSE),]

  return(cv_table)
}

get_model_details <- function(models, response_variable, lm.dat, count = 5) {
  cv_table <- cross_validate(models, response_variable, lm.dat, count)

  order <- as.numeric(rownames(cv_table))
  model_names <- sapply(1:count, function(x) paste("Model", order[x]))  # Create the name for the formula

  model_formulas <- lapply(1:count, function(x) get_model_formula(order[x], models, response_variable))
  model_formulas <- setNames(model_formulas, model_names)

  model_details <- lapply(model_formulas, function(x) {
    m <- lm(x, data = lm.dat)
    if (length(m$coefficients) > sum(complete.cases(m$coefficients))) {
      m <- lm(as.formula(paste(response_variable,"~",paste(names(m$coefficients)[complete.cases(m$coefficients)][-1],collapse="+"))),data=lm.dat)
    }

    return(m)
  })
  model_details <- setNames(model_details, model_names)
  return(model_details)
}

get_model_summaries <- function(model_details) {
  model_names <- names(model_details)

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
    summary_df$stars <-
      ifelse(summary_df$p_value < 0.001, "***",
      ifelse(summary_df$p_value < 0.01, "**",
      ifelse(summary_df$p_value < 0.05, "\\*", "")))

    return(summary_df)
  })
  model_summaries <- setNames(model_summaries, model_names)

  return(model_summaries)
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

heatmap_all_models <- function(model_details) {
  model_names <- names(model_details)

  # Create a combined dataframe of all models' coefficients and significance
  all_models_coef <- do.call(rbind, lapply(seq_along(model_details), function(i) {
    x <- model_details[[i]]

    # Get model frame and clean variable names
    mf <- model.frame(x)
    var_names <- colnames(mf)[-1] # Exclude response variable

    # Get standardized coefficients
    standardized_model <- lm(scale(model.frame(x)[,1]) ~ scale(model.frame(x)[,-1]))
    std_coef <- coef(standardized_model)

    names(std_coef) <- c("(Intercept)", var_names)

    # Get p-values
    p_values <- summary(x)$coefficients[,4]

    # Create significance stars
    stars <- ifelse(p_values < 0.001, "***",
                    ifelse(p_values < 0.01, "**",
                           ifelse(p_values < 0.05, "*", "")))

    # Include all explanatory variables
    all_vars <- data.frame(
      Model = model_names[i],
      Term = c(var_names, setdiff(explanatory_variables, var_names)),  # Include all explanatory variables
      Estimate = c(std_coef[-1], rep(NA, length(setdiff(explanatory_variables, var_names)))),  # NA for non-included vars
      Stars = c(stars[-1], rep("", length(setdiff(explanatory_variables, var_names)))),  # No stars for non-included vars
      Response = as.character(x$terms[[2]]),
      stringsAsFactors = FALSE
    )

    return(all_vars)
  }))
  all_models_coef$Model <- factor(all_models_coef$Model, levels = model_names)

  all_models_coef <- all_models_coef[all_models_coef$Term != "(Intercept)", ]

  p <- ggplot(all_models_coef, aes(x = factor(Model), y = Term, fill = Estimate)) +
    geom_tile() +
    geom_text(aes(label = Stars), color = "black", size = 3) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "Coefficient", na.value = "white") +
    theme_minimal() +
    facet_wrap(~ Response) +
    theme(axis.text.x = element_text(angle = 0),
          axis.text.y = element_text(size = 8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(x = "Model Number", y = "Variable",
         title = "Coefficients Across Top Models",
         subtitle = "* p<0.05, ** p<0.01, *** p<0.001")

  return(p)
}

get_model_formula <- function(i, object, response_variable) {
  models <- summary(object)$which
  model <- models[i, ]

  explanatory_vars <- colnames(models)[-1]
  included_vars <- explanatory_vars[model[-1] == TRUE]

  # Construct the formula
  if (length(included_vars) == 0) {
    formula_str <- paste(response_variable, "~ 1")  # Intercept-only model
  } else {
    formula_str <- paste(response_variable, "~", paste(included_vars, collapse = " + "))
  }

  return(as.formula(formula_str))
}

get_cv_error <- function(model.formula, data) {
  set.seed(1)
  train.control <- caret::trainControl(method = "cv", number = 5)
  cv <- caret::train(model.formula, data = data, method = "lm", trControl = train.control)

  return(cv$results$RMSE)
}

#Reference: https://stackoverflow.com/questions/19340277/converting-r-formula-format-to-mathematical-equation
# define a function to take a linear regression
#  (anything that supports coef() and terms() should work)
expr.from.lm <- function(fit) {
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

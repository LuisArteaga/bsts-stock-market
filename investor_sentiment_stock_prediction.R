#====#====#====#====#====#====#====# Initialization #====#====#====#====#====#====#====#==== 
library(tidyverse) # contains: ggplot2, dplyr, tidyr, readr, purr, tibble, stringr, forcats
library(rstudioapi)
library(lubridate)
library(plotly)
library(bsts)
library(feather)

options(stringsAsFactors = FALSE,  # Strings are not represented as a label in an integer form
        scipen = 999)              # Scientific Notation is deactivated

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


#====#====#====#====#====#====#====# Funktion f�r Bayesische Structural Time Series #====#====#====#====#====#====#====#====

predict_stock_prices <- function(stock_name, 
                                 start_train_date, 
                                 end_train_date, 
                                 prediction_length_days, 
                                 multivariate = FALSE, 
                                 seasonality = FALSE,
                                 seasons = 52,
                                 seed = 120784,
                                 num_iterations = 500,
                                 random_sentiment = FALSE,
                                 mape_only = FALSE) {
  
  feather_path = paste0('./import/',stock_name,'.feather')
  df <- feather::read_feather(feather_path)

  if (start_train_date < min(df$date)) {
    start_train_date <- '2018-01-01'
  }
  
  if (end_train_date >= max(df$date)) {
    end_train_date <- '2020-05-20'
    prediction_length_days = 10
  }
  
  end_prediction_date = base::as.Date(prediction_length_days, origin = end_train_date)
  start_train_year = lubridate::year(start_train_date)
  start_train_year_day = lubridate::yday(start_train_date)
  end_train_year = lubridate::year(end_train_date)
  end_train_year_day = lubridate::yday(end_train_date)
  
  
  new_df <- feather::read_feather(feather_path) %>% 
    dplyr::filter((date >= lubridate::as_date(start_train_date)) & (date <= lubridate::as_date(end_prediction_date))) %>% 
    dplyr::mutate(date = lubridate::as_date(date) ,
                  total_emotional_score = (total_positive_news_sentiment - total_negative_news_sentiment),
                  random_values = stats::rnorm(n = dplyr::n(), mean = 1, sd = 1)) %>%
    dplyr::select(date,
                  stock_price_close_ffill,
                  total_positive_news_sentiment,
                  total_negative_news_sentiment,
                  total_neutral_news_sentiment,
                  total_news_sentiment,
                  total_emotional_score,
                  stock_trading_volume,
                  random_values)
  
  if (multivariate == FALSE) {
      new_ts <- stats::ts(data = new_df$stock_price_close_ffill, 
                          frequency = 365, 
                          start = c(start_train_year,
                                    start_train_year_day + 1))
  } else {
      stock <- stats::ts(data = new_df$stock_price_close_ffill,
                         frequency = 365,
                         start = c(start_train_year,
                                   start_train_year_day + 1))
      
      if (random_sentiment == FALSE) {
        sentiment <- stats::ts(data = new_df$total_emotional_score, 
                               frequency = 365,
                               start = c(start_train_year,
                                         start_train_year_day + 1))
      } else {
        sentiment <- stats::ts(data = new_df$random_values, 
                               frequency = 365,
                               start = c(start_train_year,
                                         start_train_year_day + 1))
      }
      new_ts <- stats::ts.intersect(stock, 
                                    sentiment)
    }
    
    train_ts <- stats::window(x = new_ts, 
                              end = c(end_train_year,
                                      end_train_year_day))
    
    print('Daten wurden als DataFrame geladen.')
    
    ss <- bsts::AddSemilocalLinearTrend(list(), 
                                        y = train_ts)
    
  if (seasonality == TRUE) {
    ss <- bsts::AddSeasonal(state.specification = ss, 
                            y = train_ts, 
                            nseasons = seasons)
  }
  
  if (multivariate == FALSE) {
    model <- bsts::bsts(formula = train_ts, 
                        state.specification = ss, 
                        niter = num_iterations,
                        seed = seed)
  } else {
    model <- bsts::bsts(formula = stock ~ .,
                          state.specification = ss, 
                          niter = num_iterations, 
                          seed = seed,
                          data = train_ts,
                          expected.model.size = 2)
  }
  
  print('Modell Training ist abgeschlossen.')
  
  burn <- bsts::SuggestBurn(proportion = 0.1, 
                            bsts.object = model)
  
  prediction <- bsts::predict.bsts(model, 
                                   horizon = prediction_length_days, 
                                   newdata = train_ts,
                                   burn = burn, 
                                   quantiles = c(.025, .975))
  
  print('Modell Inferenz ist abgeschlossen.')
  
  if (multivariate == FALSE) {
    
  final_df <- data.frame(
    c(as.numeric(-colMeans(model$one.step.prediction.errors[-(1:burn),]) + train_ts),  
      as.numeric(prediction$mean[0:prediction_length_days])),
    as.numeric(new_ts),
    new_df$date)
  
  names(final_df) <- c("Fitted", "Actual", "Date")
  
  mape <- dplyr::filter(final_df, Date > end_train_date) %>% 
    summarise(MAPE = mean(abs(Actual - Fitted)/Actual))
  
  print("Der mittlere absolute prozentuale Fehler (MAPE) wurde berechnet.")
  
  prediction_intervals <- cbind.data.frame(
    as.numeric(prediction$interval[1,]),
    as.numeric(prediction$interval[2,]), 
    subset(final_df, Date >= end_train_date)$Date)
  
  names(prediction_intervals) <- c("LowerBound", "UpperBound", "Date")
  
  print("Die Vorhersageintervalle wurden berechnet.")
  
  plot_df <- left_join(final_df, prediction_intervals, by = "Date")
  
  } else {
    
    if (seasonality == FALSE) {
      final_df <- data.frame(
        c(as.numeric(-colMeans(model$one.step.prediction.errors[-(1:burn),]) + train_ts[,1]),  
          as.numeric(prediction$mean[0:prediction_length_days])),
        as.numeric(new_ts[,1]),
        new_df$date)
    } else {
      final_df <- data.frame(
        c(as.numeric(-mean(model$one.step.prediction.errors[-(1:burn)]) + train_ts[,1]),  
          as.numeric(prediction$mean[0:prediction_length_days])),
        as.numeric(new_ts[,1]),
        new_df$date)      
    }
    names(final_df) <- c("Fitted", "Actual", "Date")
    
    mape <- dplyr::filter(final_df, Date > end_train_date) %>% 
      summarise(MAPE = mean(abs(Actual - Fitted)/Actual))
    
    print("Der mittlere absolute prozentuale Fehler (MAPE) wurde berechnet.")
    
    prediction_intervals <- cbind.data.frame(
      as.numeric(prediction$interval[1,0:prediction_length_days]),
      as.numeric(prediction$interval[2,0:prediction_length_days]), 
      subset(final_df, Date >= end_train_date)$Date)
    
    names(prediction_intervals) <- c("LowerBound", "UpperBound", "Date")
    
    print("Die Vorhersageintervalle wurden berechnet.")
    
    plot_df <- left_join(final_df, prediction_intervals, by = "Date")
    
  }
  

  if (mape_only == TRUE){
    return(round(100 * mape, 2))
  }

  ggplot(data = plot_df, aes(x = Date)) +
    geom_line(aes(y = Actual, colour = "Actual"), size = 1.2) +
    geom_line(aes(y = Fitted, colour = "Fitted"), size = 1.2, linetype = 2) +
    theme_bw() + theme(legend.title = element_blank()) + ylab("") + xlab("") +
    geom_vline(xintercept = as.numeric(as.Date(end_train_date)), linetype = 2) + 
    geom_ribbon(aes(ymin = LowerBound, ymax = UpperBound), fill = "grey", alpha = 0.5) +
    ggtitle(paste0(stock_name," | ",prediction_length_days," Tage | ","Multivariate = ", multivariate , " | Seasonality = ", seasonality," | MAPE = ", round(100 * mape, 2), "%")) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0))
  
  
}


#====#====#====#====#====#====#====# MAPE Calculation #====#====#====#====#====#====#====#====


model_test <- function(stock, prediction_days){

  results1 <- c()
  
  dates <- c('2019-07-01', '2019-08-01', '2019-09-01', '2019-10-01', '2019-11-01', '2019-12-01',
             '2020-01-01', '2020-02-01', '2020-03-01')
  
  for (value in dates) {
   temp <-  predict_stock_prices(stock_name = stock, 
                                 start_train_date = '2018-01-01', 
                                 end_train_date = value, 
                                 prediction_length_days = prediction_days, 
                                 multivariate = FALSE, 
                                 seasonality = FALSE,
                                 seasons = 12,
                                 seed = 120784,
                                 num_iterations = 50,
                                 random_sentiment = FALSE,
                                 mape_only = TRUE)
   results1 <- rbind(results1, paste0(value, ": ", temp))
   
  }
  
  
  results2 <- c()
  
  
  for (value in dates) {
    temp <-  predict_stock_prices(stock_name = stock, 
                                  start_train_date = '2018-01-01', 
                                  end_train_date = value, 
                                  prediction_length_days = prediction_days, 
                                  multivariate = FALSE, 
                                  seasonality = FALSE,
                                  seasons = 12,
                                  seed = 120784,
                                  num_iterations = 250,
                                  random_sentiment = FALSE,
                                  mape_only = TRUE)
    results2 <- rbind(results2, paste0(value, ": ", temp))
    
  }
  
  
  results3 <- c()
  
  
  for (value in dates) {
    temp <-  predict_stock_prices(stock_name = stock, 
                                  start_train_date = '2018-01-01', 
                                  end_train_date = value, 
                                  prediction_length_days = prediction_days, 
                                  multivariate = TRUE, 
                                  seasonality = FALSE,
                                  seasons = 12,
                                  seed = 120784,
                                  num_iterations = 50,
                                  random_sentiment = FALSE,
                                  mape_only = TRUE)
    results3 <- rbind(results3, paste0(value, ": ", temp))
    
  }
  
  
  
  results4 <- c()
  
  
  for (value in dates) {
    temp <-  predict_stock_prices(stock_name = stock, 
                                  start_train_date = '2018-01-01', 
                                  end_train_date = value, 
                                  prediction_length_days = prediction_days, 
                                  multivariate = TRUE, 
                                  seasonality = FALSE,
                                  seasons = 12,
                                  seed = 120784,
                                  num_iterations = 50,
                                  random_sentiment = TRUE,
                                  mape_only = TRUE)
    results4 <- rbind(results4, paste0(value, ": ", temp))
  
  }
  print("Single Variable | 50 Iterations")
  print(results1)
  print("")
  print("Single Variable | 250 Iterations")
  print(results2)
  print("")
  print("Multi Variate | 50 Iterations | Stock Price and News Sentiment")
  print(results3)
  print("")
  print("Multi Variate | 50 Iterations | Stock Price and Random Numbers")
  print(results4)
}

#====#====#====#====#====#====#====# MAPE Results #====#====#====#====#====#====#====#====

stock = 'adidas'
model_test(stock = stock, prediction_days = 7)
#====#====#====#====#====#====#====# Tests #====#====#====#====#====#====#====#====

#====#====#====#====#====#====#====# Test Short Start #====#====#====#====#====#====#====#====
resultsx <- c()

datesx <- c('2019-07-01', '2019-08-01', '2019-09-01', '2019-10-01', '2019-11-01', '2019-12-01',
           '2020-01-01', '2020-02-01', '2020-03-01')

for (value in datesx) {
  temp <-  predict_stock_prices(stock_name = 'adidas', 
                                start_train_date = '2018-01-01', 
                                end_train_date = value, 
                                prediction_length_days = 30, 
                                multivariate = FALSE, 
                                seasonality = FALSE,
                                seasons = 12,
                                seed = 120784,
                                num_iterations = 500,
                                random_sentiment = FALSE,
                                mape_only = TRUE)
  
  resultsx <- rbind(resultsx, paste0(value, ": ", temp))
  
}

print(resultsx)

#====#====#====#====#====#====#====# Test Short End #====#====#====#====#====#====#====#====

predict_stock_prices(stock_name = 'adidas', 
                                                       start_train_date = '2018-01-01', 
                                                       end_train_date = '2019-08-01', 
                                                       prediction_length_days = 30, 
                                                       multivariate = TRUE, 
                                                       seasonality = TRUE,
                                                       seasons = 12,
                                                       seed = 120784,
                                                       num_iterations = 50,
                                                       random_sentiment = FALSE)





#====#====#====#====#====#====#====# Modellvergleich #====#====#====#====#====#====#====#====

bsts::CompareBstsModels(list("Modell 1 - Nur Aktienpreise" = adidas_20190901_30d_mvf_sf_rsf,
                             "Modell 2 - Aktienpreise und Sentiment" = adidas_20190901_30d_mvt_sf_rsf),
                        colors = c("red", "blue", "green"))


#====#====#====#====#====#====#====# plot_time_series #====#====#====#====#====#====#====#====
plot_time_series <- function(model, ask = TRUE) {
  ## Make all the plots callable by plot.bsts.
  opar <- par(ask = ask)
  on.exit(par(opar))
  plot.types <- c("state", "components", "residuals",
                  "prediction.errors", "forecast.distribution")
  for (plot.type in plot.types) {
    plot(model, plot.type)
  }
  if (model$has.regression) {
    regression.plot.types <- c("coefficients", "predictors", "size")
    for (plot.type in regression.plot.types) {
      plot(model, plot.type)
    }
  }
}

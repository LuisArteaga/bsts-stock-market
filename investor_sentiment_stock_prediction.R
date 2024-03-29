#====#====#====#====#====#====#====# Initialization #====#====#====#====#====#====#====#==== 
library(tidyverse) # enthalten: ggplot2, dplyr, tidyr, readr, purr, tibble, stringr, forcats
library(rstudioapi)
library(lubridate)
library(plotly)
library(bsts)
library(feather)

options(stringsAsFactors = FALSE,  # Zeichen enthalten eine Labels in Ganzzahlen Form
        scipen = 999)              # Wissenschaftliche Notation ist deaktiviert

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


#====#====#====#====#====#====#====# Funktion: Bayessche strukturelle Zeitreihenanalyse #====#====#====#====#====#====#====#====

predict_stock_prices <- function(stock_name, # Aktienname
                                 start_train_date, # Start des Trainingszeitraums
                                 end_train_date,  # Ende des Trainingszeitraums
                                 prediction_length_days, # Länge des Prognosezeitraums in Tagen
                                 multivariate = FALSE,  # Neben dem Aktienkurs eine weitere Variable, wie die Investorenstimmung berücksichtigen; Standardmäßig deaktiviert
                                 seasonality = FALSE, # Seasonalität in der Zeitreihe berücksichtigen; Standardmäßig deaktiviert
                                 seasons = 52, # Seasonalität der Zyklen für ein Jahr definieren, z.B. 52 Wochen in einem Jahr
                                 seed = 120784, # Ermöglicht die gleichen Ergebnisse bei jeder Wiederholung 
                                 num_iterations = 500, # Anzahl der Wiederholungen des Trainings
                                 random_sentiment = FALSE, # Zufallszahlen als zweite Variable im Modell nutzen
                                 mape_only = FALSE, # Lediglich die Ausgabe des MAPE Wert vom Modell forcieren
                                 model_only = FALSE) { # Lediglich die Ausgabe des Modell (als Dateiformat) forcieren
  
  feather_path = paste0('./import/',stock_name,'.feather') # Ergebnisse aus der Sentimentanalyse inkl. Aktienkurse
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
                  # Definition des Emotional Index
                  emotional_index = dplyr::if_else(is.na(total_positive_news_sentiment - total_negative_news_sentiment),
                                                          0,
                                                          total_positive_news_sentiment - total_negative_news_sentiment),
                  # Benutzerdefinierter Sentiment Index (veraltet); nicht mehr genutzt
                  positive_sentiment_share = dplyr::if_else(is.na(total_positive_news_sentiment / (total_neutral_news_sentiment + 
                                                                                                     total_positive_news_sentiment + 
                                                                                                     total_negative_news_sentiment)),
                                                            0,
                                                            round(total_positive_news_sentiment / (total_neutral_news_sentiment + total_positive_news_sentiment + total_negative_news_sentiment), digits = 2)),
                  # Definition des BSI Index (Sentiment Index)
                  bsi_index = dplyr::if_else(is.na((total_positive_news_sentiment - total_negative_news_sentiment) / total_news_sentiment),
                                             0,
                                             (total_positive_news_sentiment - total_negative_news_sentiment) / total_news_sentiment),
                  # Definition des Sent Doc (Sentiment Index)
                  sent_doc = dplyr::if_else(is.na((total_positive_news_sentiment - total_negative_news_sentiment) / (total_positive_news_sentiment + total_negative_news_sentiment)),
                                                  0,
                                                  (total_positive_news_sentiment - total_negative_news_sentiment) / (total_positive_news_sentiment + total_negative_news_sentiment)),
                  random_values = stats::rnorm(n = dplyr::n(), mean = 1, sd = 1)) %>%
    dplyr::select(date,
                  stock_price_close_ffill,
                  total_positive_news_sentiment,
                  total_negative_news_sentiment,
                  total_neutral_news_sentiment,
                  total_news_sentiment,
                  emotional_index,
                  positive_sentiment_share,
                  bsi_index,
                  sent_doc,
                  stock_trading_volume,
                  random_values)
  # Filtereinstellung, ob mehr als eine Variable berücksichtigt werden soll
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
      
      # Filtereinstellung, ob Zufallsgenerierte Zahlen genutzt werden sollen
      if (random_sentiment == FALSE) {
        sentiment <- stats::ts(data = new_df$emotional_index, #new_df$emotional_index, # new_df$bsi_index , #new_df$sent_doc,
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
    
    # Filtereinstellung, ob Seasonalität in der Zeitreihe berücksichtigt bzw. angenommen wird.
  if (seasonality == TRUE) {
    ss <- bsts::AddSeasonal(state.specification = ss, 
                            y = train_ts, 
                            nseasons = seasons)
  }
  
  # Korrektur der Datentransformation in Abhängigkeit der Eingabevariablen
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
  
  # Zur Optimierung des Modells verwendet
  burn <- bsts::SuggestBurn(proportion = 0.1, 
                            bsts.object = model)
  
  # Prognosemodell wird erstellt
  prediction <- bsts::predict.bsts(model, 
                                   horizon = prediction_length_days, 
                                   newdata = train_ts,
                                   burn = burn, 
                                   quantiles = c(.025, .975))
  
  print('Modell Inferenz ist abgeschlossen.')
  
  if (multivariate == FALSE) {
    
    # Historische Aktienkursentwicklung wird mit Prognosewerten in einer Zeitreihe verbunden
  final_df <- data.frame(
    c(as.numeric(-colMeans(model$one.step.prediction.errors[-(1:burn),]) + train_ts),  
      as.numeric(prediction$mean[0:prediction_length_days])),
    as.numeric(new_ts),
    new_df$date)
  
  names(final_df) <- c("Fitted", "Actual", "Date")
  
  # Berechnung des MAPE Wertes
  mape <- dplyr::filter(final_df, Date > end_train_date) %>% 
    summarise(MAPE = mean(abs(Actual - Fitted)/Actual))
  
  print("Der mittlere absolute prozentuale Fehler (MAPE) wurde berechnet.")
  
  # Konfidenzintervalle werden aus dem Prognosemodell extrahiert für die spätere Visualisierung
  prediction_intervals <- cbind.data.frame(
    as.numeric(prediction$interval[1,]),
    as.numeric(prediction$interval[2,]), 
    subset(final_df, Date >= end_train_date)$Date)
  
  names(prediction_intervals) <- c("LowerBound", "UpperBound", "Date")
  
  print("Die Vorhersageintervalle wurden berechnet.")
  
  plot_df <- left_join(final_df, prediction_intervals, by = "Date")
  
  } else {
    
    if (seasonality == FALSE) { 

      # Historische Aktienkursentwicklung wird mit Prognosewerten in einer Zeitreihe verbunden
      final_df <- data.frame(
        c(as.numeric(-mean(model$one.step.prediction.errors[-(1:burn)]) + train_ts[,1]),
        #c(as.numeric(-colMeans(model$one.step.prediction.errors[-(1:burn),]) + train_ts[,1]),  
          as.numeric(prediction$mean[0:prediction_length_days])),
        as.numeric(new_ts[,1]),
        new_df$date)
    } else {
      # Historische Aktienkursentwicklung wird mit Prognosewerten in einer Zeitreihe verbunden
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
  if (model_only == TRUE) {
    return(model)
  }

  # Visualiserung des Prognosemodells als Zeitreihe
  ggplot(data = plot_df, aes(x = Date)) +
    geom_line(aes(y = Actual, colour = "Actual"), size = 1.2) +
    geom_line(aes(y = Fitted, colour = "Fitted"), size = 1.2, linetype = 2) +
    theme_bw() + theme(legend.title = element_blank()) + ylab("") + xlab("") +
    geom_vline(xintercept = as.numeric(as.Date(end_train_date)), linetype = 2) + 
    geom_ribbon(aes(ymin = LowerBound, ymax = UpperBound), fill = "grey", alpha = 0.5) +
    ggtitle(paste0(stock_name," | ",prediction_length_days," Tage | ","Multivariate = ", multivariate , " | Seasonalit�t = ", seasonality," | MAPE = ", round(100 * mape, 2), "%")) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0))
  
  
}

#====#====#====#====#====#====#====# Nutzung der Funktion #====#====#====#====#====#====#====#====
resultsx <- c()

# Mögliche Eingabe Parmeter für das Startdatum des Prognosezeitraums
# datesx <- c('2019-07-01', '2019-08-01', '2019-09-01', '2019-10-01', '2019-11-01', '2019-12-01',
#           '2020-01-01', '2020-02-01', '2020-03-01')

datesx <- c('2019-12-01')

for (value in datesx) {
  temp <-  predict_stock_prices(stock_name = 'vonovia', 
                                start_train_date = '2018-01-01', 
                                end_train_date = value, 
                                prediction_length_days = 7, 
                                multivariate = TRUE, 
                                seasonality = TRUE,
                                seasons = 4,
                                seed = 120784,
                                num_iterations = 250,
                                random_sentiment = FALSE,
                                mape_only = FALSE,
                                model_only = TRUE)
  
  resultsx <- rbind(resultsx, paste0(value, ": ", temp))
  
}

print(resultsx)

#====#====#====#====#====#====#====# Modellvergleich #====#====#====#====#====#====#====#====

bsts::CompareBstsModels(list("Modell 1 - Nur Aktienpreise" = temp,
                             "Modell 2 - Aktienpreise und Sentiment" = temp2),
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


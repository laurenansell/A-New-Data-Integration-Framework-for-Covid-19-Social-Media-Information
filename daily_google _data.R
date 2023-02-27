## Code used to collect the daily google trends data for Covid-19

## Load in the libraries
library(gtrendsR)
library(tidyverse)
library(lubridate)


get_daily_gtrend <- function(keyword = c('',''), geo = 'UK', from = '2020-01-01', to = '2021-06-20') {
  if (ymd('2021-06-20') >= floor_date(Sys.Date(), 'month')) {
    to <- floor_date(ymd(to), 'month') - days(1)
    
    if (to < from) {
      stop("Specifying \'to\' date in the current month is not allowed")
    }
  }
  
  UK_data<-gtrends(keyword = "variant", geo = "GB", time = paste('2020-02-01', '2021-06-15'))
  
  
  
  mult_m <- UK_data$interest_over_time %>%
     group_by(month = floor_date(date, 'month'), keyword) %>% 
    summarise(hits = sum(as.numeric(hits))) %>%
     ungroup() %>%
     mutate(ym = format(month, '%Y-%m'),
            mult = hits / max(hits)) %>%
     select(month, ym, keyword, mult) %>%
    as_tibble()
  
  pm <- tibble(s = seq(ymd('2020-02-01'), ymd('2021-06-20'), by = 'month'), 
               e = seq(ymd('2020-02-01'), ymd('2021-06-20'), by = 'month') + months(1) - days(1))
  
  raw_trends_m <- tibble()
  
  for (i in seq(1, nrow(pm), 1)) {
    curr <- gtrends("variant", geo = "GB", time = paste(pm$s[i], pm$e[i]))
    print(paste('for', pm$s[i], pm$e[i], 'retrieved', count(curr$interest_over_time), 'days of data (all keywords)'))
    raw_trends_m <- rbind(raw_trends_m,
                          curr$interest_over_time)
  }
  
  trend_m <- raw_trends_m %>%
    select(date, keyword, hits) %>%
    mutate(ym = format(date, '%Y-%m')) %>%
    as_tibble()
  
  trend_res <- trend_m %>%
    left_join(mult_m) %>%
    mutate(est_hits = hits * mult) %>%
    select(date, keyword, est_hits) %>%
    as_tibble() %>%
    mutate(date = as.Date(date))
  
  return(trend_res)
}

get_daily_gtrend(keyword = c('Coronavirus','covid-19'), geo = 'GB', from = '2020-01-01', to = '2021-06-20')


write.csv(trend_res,"Google_trends_UK_variant.csv",row.names = FALSE)



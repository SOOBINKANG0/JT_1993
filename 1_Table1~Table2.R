## setting
rm(list=ls())
pkg = c("data.table", "tidyverse", "RPostgres", "DBI", 'glue', "frenchdata", 
        "nanoparquet", "arrow", "zoo", "lubridate", "writexl")

for(i in 1:length(pkg)){
  
  if(!require(pkg[i], character.only = T)){
    install.packages(pkg[i], dependencies = TRUE)
    require(pkg[i], character.only = T)
  }
}; rm(i, pkg)

if (!dir.exists("results")) {
  dir.create("results")
}

dat = read_parquet("raw_data/dat_monthly.parquet")

## param
J <- 6
K <- 6

sample_start <- as.Date("1965-01-01")
sample_end   <- as.Date("1989-12-31")

## custom function
calculate_cum_ret <- function(x, J) {
  
  zoo::rollapplyr(
    x,
    width = J,
    FUN = function(z) {
      if (sum(!is.na(z)) < J) {
        NA_real_
      } else {
        prod(1 + z) - 1
      }
    },
    fill = NA_real_,
    partial = FALSE
  )
}

## make signal 1
J_temp = c(3,6,9,12)

dat <- dat %>%
  arrange(permno, date) %>%
  group_by(permno) %>%
  group_modify(~ {
    
    out <- .x
    
    for (J in J_temp) {
      out[[paste0("cum_ret_", J)]] <-
        calculate_cum_ret(out$mthret, J)
    }
    
    out
  }) %>%
  ungroup()

rm(calculate_cum_ret)

## make signal(J)
dat = dat %>% group_by(date) %>% mutate(
  st3 = ntile(cum_ret_3, 10),
  st6 = ntile(cum_ret_6, 10),
  st9 = ntile(cum_ret_9, 10),
  st12 = ntile(cum_ret_12, 10)
) %>% ungroup()

dat = dat %>% mutate(date = floor_date(date, "month"))

temp = frenchdata::download_french_data('Fama/French 3 Factors')
temp = as.data.frame(temp[3]$subsets$data[1])

temp = temp %>% mutate(Mkt.RF = Mkt.RF/100,
                       SMB = SMB/100,
                       HML = HML/100,
                       RF = RF/100)

temp <- temp %>% rename(ym = date) %>% mutate(date = as.Date(paste0(ym, "01"), format = "%Y%m%d"))
temp = temp %>% select(-ym) %>% select(date, Mkt.RF, RF)

dat = dat %>% left_join(temp, by = "date")
rm(temp)

##
cohort <- dat %>%
  select(date, permno, st3, st6, st9, st12) %>%
  pivot_longer(
    cols = c(st3, st6, st9, st12),
    names_to = "J",
    names_pattern = "st(\\d+)",
    names_transform = list(J = as.integer),
    values_to = "momr"
  ) %>%
  filter(!is.na(momr)) 

##---------------------------
## Table 1, Panel A
##---------------------------
res_table = data.frame()
J_temp = c(3,6,9,12)
K_temp = c(3,6,9,12)

for(j in K_temp){

  for(i in J_temp){

    cohort_6_6 <- cohort %>%
      filter(J == i) %>%
      mutate(K = j) %>%
      uncount(K, .remove = FALSE, .id = "holding_month")

    cohort_6_6 = cohort_6_6 %>% arrange(date) %>%
      mutate(holding_date = date %m+% months(holding_month))

    temp = dat %>% select(date, permno, mthret) %>% rename(holding_date = date)
    cohort_6_6 = cohort_6_6 %>% left_join(temp, by = c("holding_date", "permno"))
    rm(temp)

    ##
    cohort_6_6 = cohort_6_6 %>% group_by(holding_date, date, momr) %>%
      summarise(cohort_ret = if (
        all(is.na(mthret))
      ) {
        NA_real_
      } else {
        mean(mthret, na.rm = TRUE)
      }, .groups = "drop")

    portfolio_ret <- cohort_6_6 %>%
      group_by(holding_date, momr) %>%
      summarise(
        portfolio_ret = mean(cohort_ret, na.rm = TRUE),
        n_cohorts = sum(!is.na(cohort_ret)),
        .groups = "drop"
      ) %>% ungroup()

    rm(cohort_6_6)

    ##
    results = portfolio_ret %>% filter(holding_date >= as.Date(sample_start)) %>%
      filter(holding_date <= as.Date(sample_end)) %>% filter(momr == 1 | momr == 10) %>%
      group_by(momr) %>% select(-n_cohorts) %>%
      pivot_wider(names_from = momr, values_from = portfolio_ret) %>%
      rename(Winner=`10`, Loser = `1`) %>% mutate(LS = Winner - Loser) %>% select(-holding_date)

    res = apply(results, 2, FUN = function(x){
      test = t.test(x)
      c(
        mean = unname(round(test$estimate,4)),
        t_value = unname(round(test$statistic,4)))
    } ) %>% t() %>% as.data.frame() %>%
      tibble::rownames_to_column("portfolio") %>%
      mutate(
        J = i,
        K = j,
        .before = 1
      )

    res_table = bind_rows(res_table, res)
  }
  rm(res, results)
}

res_table <- res_table %>%
  pivot_wider(
    id_cols = c(J, portfolio),
    names_from = K,
    values_from = c(mean, t_value),
    names_glue = "{.value}_K{K}"
  )

write_xlsx(res_table, "results/Table1.xlsx")
rm(j,J,J_temp,K,K_temp)
##---------------------------
## Table 2
##---------------------------

cohort_dat <- cohort %>%
  filter(J == 6) %>%
  mutate(K = 6) %>%
  uncount(K, .remove = FALSE, .id = "holding_month") %>% arrange(date) %>% 
  mutate(holding_date = date %m+% months(holding_month))

temp = dat %>% select(date, permno, mthret, mthcap, mthprevcap, Mkt.RF, RF) %>% rename(holding_date = date)
cohort_dat = cohort_dat %>% left_join(temp, by = c("holding_date", "permno"))
rm(temp)

cohort_dat = cohort_dat %>% mutate(ri_rf = mthret - RF)

##
cohort_66 = cohort_dat %>% group_by(holding_date, date, momr) %>% 
  summarise(cohort_ret = if (
    all(is.na(mthret))
  ) {
    NA_real_
  } else {
    mean(mthret, na.rm = TRUE)
  }, 
  mthprevcap = mean(mthprevcap, na.rm = T)
  ,
  .groups = "drop")

portfolio_ret <- cohort_66 %>%
  group_by(holding_date, momr) %>%
  summarise(
    portfolio_ret = mean(cohort_ret, na.rm = TRUE),
    mthprevcap = mean(mthprevcap, na.rm =T),
    n_cohorts = sum(!is.na(cohort_ret)),
    .groups = "drop"
  ) %>% ungroup() %>% filter(holding_date >= sample_start) %>% filter(holding_date <= sample_end)

## Fama french data loading
fama = frenchdata::download_french_data('Fama/French 3 Factors')
fama = as.data.frame(fama[3]$subsets$data[1])

fama = fama %>% mutate(Mkt.RF = Mkt.RF/100,
                       SMB = SMB/100,
                       HML = HML/100,
                       RF = RF/100)

fama <- fama %>% rename(ym = date) %>% mutate(date = as.Date(paste0(ym, "01"), format = "%Y%m%d"))
fama = fama %>% select(-ym) %>% select(date, Mkt.RF, RF) %>% rename(holding_date = date)

##
res_table = data.frame(momr = seq(1:10), beta = NA)

for(i in 1:10){
  temp = portfolio_ret %>% filter(momr == i) %>% select(-momr,-n_cohorts)
  temp = temp %>% left_join(fama, by = "holding_date") %>% mutate(ri_rf = portfolio_ret - RF) %>% 
    mutate(mkt_lag = lag(Mkt.RF,1),
           mkt_lead = lead(Mkt.RF,1)) 
  
  res = lm(data= temp, formula = ri_rf ~ Mkt.RF)$coef[2]
  res_table[i,2] <- round(res, 2)
}

Avg_cap = portfolio_ret %>% group_by(momr) %>% summarise(avg_mthcap = mean(mthprevcap/1000000, na.rm = T))

table2 = left_join(res_table, Avg_cap, by = "momr")
write_xlsx(table2, "results/table2.xlsx")

portfolio_ret = portfolio_ret %>% left_join(res_table, by = "momr") %>% select(-n_cohorts)
rm(temp, res_table, table2, Avg_cap, table2, Avg_cap, fama)
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
daily_dat = read_parquet("raw_data/crsp_data.parquet")

dat = dat %>% mutate(date = floor_date(date, "month"))

## param
sample_start <- as.Date("1965-01-01")
sample_end   <- as.Date("1989-12-31")

##
## Fama french data loading
fama = frenchdata::download_french_data('Fama/French 3 Factors')
fama = as.data.frame(fama[3]$subsets$data[1])

fama = fama %>% mutate(Mkt.RF = Mkt.RF/100,
                       SMB = SMB/100,
                       HML = HML/100,
                       RF = RF/100) %>% select(date, Mkt.RF, RF) %>% mutate(MKT = Mkt.RF + RF) %>% 
  rename(ym = date) %>% mutate(date = as.Date(paste0(ym, "01"), format = "%Y%m%d")) %>% 
  select(-ym)

dat = dat %>% left_join(fama, by = "date")


fama = frenchdata::download_french_data('Fama/French 3 Factors [Daily]' )
fama = as.data.frame(fama[3]$subsets$data[1])

fama <- fama %>%
  mutate(
    date = ymd(
      as.character(date)
    ),
    Mkt.RF = Mkt.RF / 100,
    SMB = SMB / 100,
    HML = HML / 100,
    RF = RF / 100,
    MKT = Mkt.RF + RF
  ) %>%
  select(
    date,
    Mkt.RF,
    RF,
    MKT
  )

daily_dat = daily_dat %>% left_join(fama, by = "date")
daily_dat = daily_dat %>% filter(is.na(dlyret) == F)

rm(fama)

##-----------------------------
## Table 3
##-----------------------------
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

J_temp = c(6)

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

##
dat = dat %>% arrange(date, permno) %>% 
  group_by(date) %>% mutate(st6 = ntile(cum_ret_6, 10)) %>% ungroup() # 전체표본 mom
  
temp = dat %>% filter(is.na(cum_ret_6) == F, is.na(mthcap) == F, mthcap >0) %>% arrange(date, permno) %>% group_by(date) %>% mutate(size_tile = ntile(mthcap, 3)) %>% 
  ungroup() %>% group_by(size_tile, date) %>% mutate(size_momr = ntile(cum_ret_6, 10)) %>% 
  select(date, permno, size_tile, size_momr)

dat = dat %>% left_join(temp, by = c("date", "permno"), relationship = "one-to-one"); rm(temp)

## holding return
temp = dat %>% select(date, permno, mthret) %>% rename(holding_date = date, holding_ret = mthret)

cohort_dat = dat %>% filter(!is.na(st6)) %>% select(date, permno, st6, size_tile, size_momr) %>% mutate(K = 6L) %>%
  uncount(K, .remove = T, .id = "holding_month") %>% 
  mutate(holding_date = date %m+% months(holding_month)) %>% ungroup()

cohort_dat <- cohort_dat %>%
  left_join(
    temp,
    by = c(
      "holding_date",
      "permno"
    )
  )
rm(temp)

cohort_66_all <- cohort_dat %>%
  filter(
    !is.na(st6)
  ) %>%
  group_by(
    holding_date,
    date,
    st6
  ) %>%
  summarise(
    cohort_ret = if (
      all(!is.finite(holding_ret))
    ) {
      NA_real_
    } else {
      mean(
        holding_ret[
          is.finite(holding_ret)
        ]
      )
    },
    .groups = "drop"
  )


## portfolio 단위
portfolio_ret_all <- cohort_66_all %>% group_by(
    holding_date,
    st6
  ) %>% summarise(
    portfolio_ret = if (
      all(!is.finite(cohort_ret))
    ) {
      NA_real_
    } else {
      mean(
        cohort_ret[
          is.finite(cohort_ret)
        ]
      )
    },
    .groups = "drop"
  ) %>%
  filter(
    holding_date >= sample_start,
    holding_date <= sample_end
  ) %>%
  arrange(
    holding_date,
    st6
  )

table3_all <- portfolio_ret_all %>%
  group_by(st6) %>%
  summarise(
    mean = mean(
      portfolio_ret
    ),
    
    t_value = mean(
      portfolio_ret
    ) /
      (
        sd(portfolio_ret) /
          sqrt(n())
      ),
    .groups = "drop"
  ) %>%
  mutate(
    mean = round(
      mean,
      4
    ),
    
    t_value = round(
      t_value,
      4
    )
  )

table3_all_spread <- portfolio_ret_all %>%
  filter(
    st6 %in% c(1, 10)
  ) %>%
  select(
    holding_date,
    st6,
    portfolio_ret
  ) %>%
  pivot_wider(
    id_cols = holding_date,
    names_from = st6,
    values_from = portfolio_ret,
    names_prefix = "P"
  ) %>%
  mutate(
    `P10-P1` = P10 - P1
  ) %>%
  filter(
    is.finite(`P10-P1`)
  ) %>%
  arrange(holding_date)

test <- t.test(
  table3_all_spread$`P10-P1`,
  mu = 0
)

table3_all_spread_res <- data.frame(
  st6 = "P10-P1",
  
  mean = round(
    unname(test$estimate),
    4
  ),
  
  t_value = round(
    unname(test$statistic),
    4
  )
)

temp = rbind(table3_all, table3_all_spread_res)


table3_all_wide =portfolio_ret_all %>% 
  pivot_wider(id_cols = holding_date, names_from = st6, names_prefix = "P", 
              values_from = portfolio_ret) %>% 
  arrange(holding_date)

diff_mat <- table3_all_wide %>%
  transmute(across(P2:P10, ~ .x - P1)) %>%
  as.matrix()

## 
fit <- stats::manova(diff_mat ~ 1)
fit2 <- summary(fit, test = "Wilks", intercept = TRUE)  

temp = temp %>% rbind(data.frame(st6 = "F-stat", 
                          mean = fit2$stats["(Intercept)","approx F"],
                          t_value = round(fit2$stats["(Intercept)","Pr(>F)"],3))
)
write_xlsx(temp, "results/table3_all.xlsx")

rm(table3_all, table3_all_spread_res, test, table3_all_spread, 
   table3_all_wide,fit, fit2, temp, diff_mat, cohort_66_all)

##--------------------------------
## Table 3:Size
##--------------------------------
cohort_66_size = cohort_dat %>% group_by(size_tile, size_momr,holding_date, date) %>% 
  summarise(cohort_ret = mean(holding_ret, na.rm = T))

portfolio_ret_size <- cohort_66_size %>%
  group_by(
    holding_date,
    size_tile,
    size_momr
  ) %>%
  summarise(
    portfolio_ret = if (all(!is.finite(cohort_ret))) {
      NA_real_
    } else {
      mean(cohort_ret[is.finite(cohort_ret)])
    },
    n_cohorts = sum(is.finite(cohort_ret)),
    .groups = "drop"
  ) %>%
  filter(
    holding_date >= sample_start,
    holding_date <= sample_end
  )

table3_size <- portfolio_ret_size %>%
  filter(is.finite(portfolio_ret)) %>%
  group_by(size_tile, size_momr) %>%
  summarise(
    mean = mean(portfolio_ret),
    t_value = mean(portfolio_ret) /
      (sd(portfolio_ret) / sqrt(n())),
    n_months = n(),
    .groups = "drop"
  ) %>% select(-n_months)

##
temp = portfolio_ret_size %>% filter(size_tile %in% 1) %>% select(-n_cohorts, -size_tile) %>% 
  pivot_wider(
    id_cols = holding_date,
    names_from = size_momr,
    names_prefix = "P",
    values_from = 
      portfolio_ret
    )

diff_mat <- temp %>%
  transmute(across(P2:P10, ~ .x - P1)) %>%
  as.matrix()

fit <- stats::manova(diff_mat ~ 1)
fit2 <- summary(fit, test = "Wilks", intercept = TRUE)  


temp2 = temp %>% select(holding_date, P1, P10) %>% mutate(P10-P1)
test = t.test(temp2$`P10 - P1`)
temp3 = data.frame(size_tile = 1, size_momr = "P10-P1", 
                   mean = unname(round(test$estimate,4)), t_value = unname(round(test$statistic,2)))
table3_size = table3_size %>% rbind(temp3)

table3_size = table3_size %>% rbind(data.frame(size_tile = 1, size_momr = "F-stat", 
                                               mean = fit2$stats["(Intercept)","approx F"],
                                               t_value = round(fit2$stats["(Intercept)","Pr(>F)"],3))
)

rm(temp, temp2, temp3, tempy)
##
temp = portfolio_ret_size %>% filter(size_tile %in% 2) %>% select(-n_cohorts, -size_tile) %>% 
  pivot_wider(
    id_cols = holding_date,
    names_from = size_momr,
    names_prefix = "P",
    values_from = 
      portfolio_ret
  )

diff_mat <- temp %>%
  transmute(across(P2:P10, ~ .x - P1)) %>%
  as.matrix()

fit <- stats::manova(diff_mat ~ 1)
fit2 <- summary(fit, test = "Wilks", intercept = TRUE)  

temp2 = temp %>% select(holding_date, P1, P10) %>% mutate(P10-P1)
test = t.test(temp2$`P10 - P1`)
temp3 = data.frame(size_tile = 2, size_momr = "P10-P1", 
                   mean = unname(round(test$estimate,4)), t_value = unname(round(test$statistic,2)))
table3_size = table3_size %>% rbind(temp3)

table3_size = table3_size %>% rbind(data.frame(size_tile = 2, size_momr = "F-stat", 
                                               mean = fit2$stats["(Intercept)","approx F"],
                                               t_value = round(fit2$stats["(Intercept)","Pr(>F)"],3))
)

##
temp = portfolio_ret_size %>% filter(size_tile %in% 3) %>% select(-n_cohorts, -size_tile) %>% 
  pivot_wider(
    id_cols = holding_date,
    names_from = size_momr,
    names_prefix = "P",
    values_from = 
      portfolio_ret
  )

diff_mat <- temp %>%
  transmute(across(P2:P10, ~ .x - P1)) %>%
  as.matrix()

fit <- stats::manova(diff_mat ~ 1)
fit2 <- summary(fit, test = "Wilks", intercept = TRUE)  

temp2 = temp %>% select(holding_date, P1, P10) %>% mutate(P10-P1)
test = t.test(temp2$`P10 - P1`)
temp3 = data.frame(size_tile = 2, size_momr = "P10-P1", 
                   mean = unname(round(test$estimate,4)), t_value = unname(round(test$statistic,2)))
table3_size = table3_size %>% rbind(temp3)

table3_size = table3_size %>% rbind(data.frame(size_tile = 3, size_momr = "F-stat", 
                                               mean = fit2$stats["(Intercept)","approx F"],
                                               t_value = round(fit2$stats["(Intercept)","Pr(>F)"],3))
)
write_xlsx(table3_size, "results/table3_size.xlsx")
rm(temp, temp2, temp3, test, table3_size, fit, tit2, diff_mat)
##-----------------------------
## Table 3: beta
##-----------------------------
## review jt1993_table3_full_replication.R
## omit

##-----------------------------
## Table 4: ALL
##-----------------------------
Table4_all = portfolio_ret_all %>% filter(st6 %in% c(1,10)) %>% 
  pivot_wider(id_cols = holding_date, names_from = st6, names_prefix = "P", values_from = portfolio_ret) %>% 
  mutate(`P10-P1` = P10-P1) %>% select(holding_date, `P10-P1`) %>% group_by(months(holding_date)) %>% 
  summarise(ret = mean(`P10-P1`, na.mr = T),
            tstat = mean(`P10-P1`)/(sd(`P10-P1`)/sqrt(n()))) %>% arrange(`months(holding_date)`) %>% 
  rename(month = `months(holding_date)`)

Table4_all$month = str_remove_all(Table4_all$month, "월") %>% as.numeric() 
Table4_all = Table4_all %>% arrange(month)

Table4_spread = portfolio_ret_all %>% filter(st6 %in% c(1,10)) %>% 
  pivot_wider(id_cols = holding_date, names_from = st6, names_prefix = "P", values_from = portfolio_ret) %>% 
  mutate(`P10-P1` = P10-P1) %>% select(holding_date, `P10-P1`)

Table4_febdec = Table4_spread %>% filter(months(holding_date) != "1월") %>% 
  summarise(month = "Feb-Dec", ret = mean(`P10-P1`, na.mr = T),
            tstat = mean(`P10-P1`)/(sd(`P10-P1`)/sqrt(n())))

Table4_wide <- Table4_spread %>%
  mutate(
    year = year(holding_date),
    mon = factor(
      month(holding_date),
      levels = 1:12,
      labels = month.abb
    )
  ) %>%
  select(year, mon, `P10-P1`) %>%
  pivot_wider(
    id_cols = year,
    names_from = mon,
    values_from = `P10-P1`
  )  

diff_mat <- Table4_wide %>%
  transmute(across(Feb:Dec, ~ .x - Jan)) %>%
  as.matrix()

fit <- stats::manova(diff_mat ~ 1)
fit2 <- summary(fit, test = "Wilks", intercept = TRUE)  

diff_mat2 = Table4_wide %>% select(-Jan) %>%
  transmute(across(Mar:Dec, ~ .x - Feb)) %>%
  as.matrix()

fit3 <- stats::manova(diff_mat2 ~ 1)
fit4 <- summary(fit3, test = "Wilks", intercept = TRUE)  

temp = data.frame(month = "F-stat", ret = round(fit2$stats["(Intercept)","approx F"],4),
                                 tstat = round(fit2$stats["(Intercept)","Pr(>F)"],3))

temp2 = data.frame(month = "F-stat", ret = round(fit4$stats["(Intercept)","approx F"],4),
                  tstat = round(fit4$stats["(Intercept)","Pr(>F)"],3))

temp2 = rbind(Table4_all,Table4_febdec,temp,temp2)
write_xlsx(temp2, "results/table4_all.xlsx")

rm(fit, fit2, temp, temp2, diff_mat)
##--------------------------------------------------
## Table IV: S1, S2, S3
##--------------------------------------------------
Table4_size = portfolio_ret_size %>% filter(size_momr %in% c(1,10)) %>% 
  pivot_wider(id_cols = c(holding_date, size_tile), names_from = size_momr, names_prefix = "P", values_from = portfolio_ret) %>% 
  mutate(`P10-P1` = P10-P1) 

size1 = Table4_size %>% select(-P1, -P10) %>% group_by(months(holding_date), size_tile) %>% 
  summarise(ret = round(mean(`P10-P1`, na.mr = T),4),
            tstat = mean(`P10-P1`)/(sd(`P10-P1`)/sqrt(n()))) %>% arrange(`months(holding_date)`) %>% 
  rename(month = `months(holding_date)`)
size1$month = str_remove_all(size1$month, "월") %>% as.numeric() 
size1 = size1 %>% arrange(size_tile, month); size1$month = as.character(size1$month)

## Feb-dec
size_spread1 = Table4_size %>% filter(months(holding_date) != "1월") %>% select(-P1, -P10) %>% 
  group_by(size_tile) %>% summarise(month = "Feb-Dec", ret = round(mean(`P10-P1`, na.mr = T),4),
                                               tstat = mean(`P10-P1`)/(sd(`P10-P1`)/sqrt(n())))
size_spread1 = size_spread1 %>% relocate(size_tile, .after = month)

temp = bind_rows(size1, size_spread1); gc()

## f-test
Table4_size_f <- Table4_size %>%
  mutate(
    mon = month(holding_date)
  ) %>% group_by(size_tile) %>%
  group_modify(~ {
    
    fit_a <- aov(
      `P10-P1` ~ factor(mon),
      data = .x
    )
    result_a <- summary(fit_a)[[1]]

    fit_b <- aov(
      `P10-P1` ~ factor(mon),
      data = filter(.x, mon != 1)
    )
    result_b <- summary(fit_b)[[1]]
    
    tibble(
      test = c("Jan-Dec", "Feb-Dec"),
      ret = c(
        result_a[["F value"]][1],
        result_b[["F value"]][1]
      ),
      tstat = c(
        round(result_a[["Pr(>F)"]][1],3),
        round(result_b[["Pr(>F)"]][1],3)
      )
    )
  }) %>%
  ungroup() %>%
  arrange(size_tile, test)

Table4_size_f = Table4_size_f %>% rename(month = `test`) %>% relocate(month, .before = size_tile)
res = rbind(temp, Table4_size_f)
write_xlsx(res, "results/Table4_size.xlsx")

##------------------------------
## Table V
##------------------------------

t1 = Table4_spread %>% group_by(months(holding_date)) %>% summarise(win = mean(`P10-P1`>0))
t1$win = round(t1$win, 2); t1 = t1 %>% rename(month = `months(holding_date)`)
t2 = Table4_spread %>% filter(months(holding_date) != "1월") %>% summarise(win = mean(`P10-P1`>0))
t2$win = round(t2$win, 2); t2 = data.frame(month = "Feb-Dec", win = t2$win)

t3 = Table4_spread %>% summarise(win = mean(`P10-P1`>0))
t3$win = round(t3$win, 2); t3 = data.frame(month = "All", win = t3$win)

temp1 = rbind(t1, t2,t3)

t4 = Table4_size %>% select(-P1, -P10) %>% group_by(size_tile, months(holding_date)) %>% 
  summarise(win = round(mean(`P10-P1`>0),2))
t4 = t4 %>% rename(month = `months(holding_date)`)
t5 = Table4_size %>% select(-P1, -P10) %>% filter(months(holding_date) != "1월") %>% 
  group_by(size_tile)%>% summarise(win = round(mean(`P10-P1`>0),2))
t5 = data.frame(t5, month = "Feb-Dec")
t6 = Table4_size %>% select(-P1, -P10) %>%  group_by(size_tile) %>% summarise(win = mean(`P10-P1`>0))
t6 = data.frame(t6, month = "All")

temp2 = rbind(t4,t5,t6)

write_xlsx(temp1, "results/table5_all.xlsx")
write_xlsx(temp2, "results/table5_size.xlsx")
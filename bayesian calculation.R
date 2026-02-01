
setwd('C:/Users/user/Desktop/R files')

suppressPackageStartupMessages({
  library(readr); library(readxl); library(dplyr); library(stringr)
  library(lubridate); library(anytime); library(zoo)
  library(rstan); library(ggplot2); library(tidyr)
})
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

PATH_RT    <- "Rt.csv"                       # 날짜, Rt(관측)
PATH_R0T   <- "R(0,t).csv"                   # 날짜, R0_t
PATH_NPI   <- "RoOI_adjusted.xlsx"           # 날짜, RoOI_adjusted(=N)
PATH_VACC  <- "practical vaccination rate.xlsx" # 날짜, practical vaccination rate(=V)
PATH_HUMID <- "humidity_prepared_daily.xlsx" # 날짜, Hs 또는 humidity_pct



add_date <- function(df) {
  cn <- names(df)
  cand <- cn[grepl("date|날짜|일자", tolower(cn))]
  col  <- if (length(cand)) cand[1] else cn[1]
  d    <- suppressWarnings(anytime::anytime(df[[col]]))
  df$date <- as.Date(d)
  df
}
pick_numeric_col <- function(df, exclude = "date") {
  cs <- setdiff(names(df), exclude)
  if (!length(cs)) return(NULL)
  cs[ which.max(sapply(cs, function(c) sum(!is.na(suppressWarnings(as.numeric(df[[c]]))))) ) ]
}




load_rt <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE) |> add_date()
  val <- pick_numeric_col(df)
  tibble(date=df$date, Rt=as.numeric(df[[val]])) |> tidyr::drop_na()
}
load_r0t <- function(path) {
  df <- if (grepl("\\.csv$", path, TRUE)) readr::read_csv(path, show_col_types = FALSE) else readxl::read_excel(path)
  df <- add_date(df)
  hit <- names(df)[tolower(names(df)) %in% c("r0_t","r0t","r(0,t)","r0,t")]
  col <- if (length(hit)) hit[1] else pick_numeric_col(df)
  tibble(date=df$date, R0_t=as.numeric(df[[col]])) |> tidyr::drop_na()
}
load_n <- function(path) {
  df <- readxl::read_excel(path) |> add_date()
  ncol <- names(df)[ stringr::str_detect(tolower(names(df)), "rooi") ]
  ncol <- if (length(ncol)) ncol[1] else pick_numeric_col(df)
  n <- as.numeric(df[[ncol]])
  if (max(n, na.rm=TRUE) > 1.5) n <- n/100
  tibble(date=df$date, N=pmin(pmax(n,0),1)) |> filter(!is.na(date))
}



load_vacc_strict <- function(path, v_pattern="practical.*vaccination|\\bvaccination\\b") {
  hdr <- readxl::read_excel(path, n_max = 0)
  cols <- names(hdr)
  d_idx <- which(grepl("date|날짜|일자", tolower(cols))); d_idx <- if (length(d_idx)) d_idx[1] else 1L
  v_idx <- which(grepl(v_pattern, tolower(cols))); stopifnot("접종률 열을 찾지 못했습니다." = length(v_idx)>=1); v_idx <- v_idx[1]
  types <- rep("guess", length(cols)); types[d_idx] <- "date"; types[v_idx] <- "numeric"
  df <- readxl::read_excel(path, col_types = types) |> add_date()
  v_raw <- df[[v_idx]]
  if (inherits(v_raw, "Date") || inherits(v_raw, "POSIXct")) stop("접종률 열이 날짜형으로 읽혔습니다. 엑셀에서 숫자/백분율로 변경하세요.")
  v <- as.numeric(v_raw); if (max(v, na.rm=TRUE) > 1.5) v <- v/100
  tibble(date=df$date, V=v) |> filter(!is.na(date))
}
load_hs <- function(path) {
  df <- readxl::read_excel(path) |> add_date()
  cn <- tolower(names(df))
  if ("hs" %in% cn) {
    hs <- as.numeric(df[[ which(cn=="hs") ]])
  } else {
    idx <- which(stringr::str_detect(cn, "humidity_pct"))
    s <- as.numeric(df[[ if (length(idx)) idx[1] else pick_numeric_col(df) ]])
    mu <- mean(s, na.rm=TRUE); sd <- sd(s, na.rm=TRUE); if (is.na(sd)||sd==0) sd <- 1
    hs <- (s - mu)/sd
  }
  tibble(date=df$date, Hs=hs) |> filter(!is.na(date))
}



rt   <- load_rt(PATH_RT)
r0t  <- load_r0t(PATH_R0T)
npi  <- load_n(PATH_NPI)
vacc <- load_vacc_strict(PATH_VACC)   # ★ 숫자 강제
hum  <- load_hs(PATH_HUMID)

dat <- rt |>
  inner_join(r0t, by="date") |>
  left_join(npi,  by="date") |>
  left_join(vacc, by="date") |>
  left_join(hum,  by="date") |>
  filter(date >= as.Date("2020-05-01"), date <= as.Date("2022-07-31")) |>
  arrange(date)



dat <- dat |> mutate(V = if_else(date < as.Date("2021-03-01"), 0, V))



dat <- dat |>
  mutate(
    N  = zoo::na.locf(N, na.rm = FALSE),
    V  = zoo::na.locf(V, na.rm = FALSE),
    Hs = zoo::na.approx(Hs, maxgap = 7, na.rm = FALSE)
  )



stopifnot(sum(is.na(dat$Rt))==0, sum(is.na(dat$R0_t))==0)



dat <- dat |>
  mutate(month = floor_date(date, "month"),
         m_idx  = as.integer(as.factor(month)))
Tn <- nrow(dat); L <- length(unique(dat$m_idx))




stan_code <- "
data{
  int<lower=1> T;
  vector<lower=0>[T] Rt;
  vector<lower=0>[T] R0_t;
  vector<lower=0,upper=1>[T] N;
  vector<lower=0,upper=1>[T] V;
  vector[T] Hs;
  int<lower=1> L;
  int<lower=1,upper=L> m_idx[T];
}
parameters{
  real<lower=0,upper=1> u_alpha;
  real<lower=0,upper=1> u_beta;
  real<lower=0,upper=1> u_lambda;

  real<lower=0> alpha;
  real<lower=0> beta;
  real<lower=0> lambda;

  real<lower=0> k_phi;
  real phi;

  real<lower=0> sigma_D;
  vector[L] Delta;

  real<lower=0> kappa;
}

model{
  // hyperpriors
  u_alpha ~ uniform(0,1);
  u_beta  ~ uniform(0,1);
  u_lambda~ uniform(0,1);

  // priors
  alpha ~ gamma(u_alpha, 1);
  beta  ~ gamma(u_beta,  1);
  lambda~ gamma(u_lambda,1);

  k_phi  ~ normal(0,0.3);
  phi    ~ normal(0,k_phi);

  sigma_D ~ normal(0,0.3);
  Delta   ~ normal(0,sigma_D);

  kappa ~ normal(0,0.5);

  // likelihood
  for (t in 1:T){
    real log_mu = log(R0_t[t])
                 - alpha * N[t]
                 - beta  * V[t]
                 - lambda * N[t]*V[t]
                 - phi   * Hs[t]
                 - Delta[m_idx[t]];
    Rt[t] ~ gamma(kappa, kappa/exp(log_mu)); // mean-param
  }
}
generated quantities {
  vector[T] mu;
  for (t in 1:T){
    mu[t] = exp( log(R0_t[t])
                 - alpha * N[t]
                 - beta  * V[t]
                 - lambda * N[t]*V[t]
                 - phi   * Hs[t]
                 - Delta[m_idx[t]] );
  }
}
"

sm <- rstan::stan_model(model_code = stan_code)
stan_data <- list(
  T=Tn,
  Rt=as.vector(dat$Rt),
  R0_t=as.vector(dat$R0_t),
  N=as.vector(dat$N),
  V=as.vector(dat$V),
  Hs=as.vector(dat$Hs),
  L=L,
  m_idx=as.integer(dat$m_idx)
)

set.seed(2025)
fit <- rstan::sampling(
  sm, data=stan_data,
  iter=4000, warmup=2000, chains=4,
  control=list(adapt_delta=0.95, max_treedepth=12)
)





print(fit, pars=c("alpha","beta","lambda","phi","sigma_D","kappa"), probs=c(.025,.5,.975))
mu_draws <- rstan::extract(fit, pars="mu")$mu    # iter x T
mu_med <- apply(mu_draws, 2, median)
mu_lo  <- apply(mu_draws, 2, quantile, 0.025)
mu_hi  <- apply(mu_draws, 2, quantile, 0.975)

res <- dat |>
  transmute(date, Rt, R0_t, N, V, Hs,
            mu_med=mu_med, mu_lo=mu_lo, mu_hi=mu_hi)

readr::write_csv(res, "bayes_rt_results.csv")
coef_sum <- as.data.frame(rstan::summary(fit, pars=c("alpha","beta","lambda","phi","sigma_D","kappa"))$summary)
coef_sum$term <- rownames(coef_sum)
coef_sum <- coef_sum |> select(term, mean, `2.5%`, `50%`, `97.5%`, Rhat, n_eff=`n_eff`)
readr::write_csv(coef_sum, "bayes_coefficients.csv")

p <- ggplot(res, aes(date)) +
  geom_line(aes(y=Rt, colour="Observed Rt"), linewidth=0.7, alpha=0.85) +
  geom_ribbon(aes(ymin=mu_lo, ymax=mu_hi, fill="Posterior 95% CrI"), alpha=.25, colour=NA) +
  geom_line(aes(y=mu_med, colour="Posterior median μ_t"), linewidth=1.0) +
  geom_hline(yintercept=1, linetype="dashed") +
  scale_colour_manual(values=c("Observed Rt"="#333333","Posterior median μ_t"="#1B4F72"), name=NULL) +
  scale_fill_manual(values=c("Posterior 95% CrI"="#5DADE2"), name=NULL) +
  labs(title="South Korea: Bayesian fit of Rt (2020-05 to 2022-07)",
       subtitle="Median (line) with 95% CrI (band); model uses Hs (humidity), N, V, N×V, month effect; offset=log(R0_t)",
       x="Date", y="Rt") +
  theme_minimal(base_size = 12) +
  theme(legend.position="top")
ggsave("bayes_rt_timeseries.png", p, width=20, height=5, dpi=150)
print(p)

cat("Saved:\n  - bayes_rt_results.csv\n  - bayes_coefficients.csv\n  - bayes_rt_timeseries.png\n")

library(ggplot2)
library(stringr)

title_txt <- "South Korea: Bayesian fit of Rt (2020-05 to 2022-07)"
sub_txt   <- "Median (line) with 95% CrI (band); model uses Hs (humidity), N, V, N×V, month effect; offset=log(R0_t)"

p2 <- p +
  labs(
    title    = str_wrap(title_txt, 70),         
    subtitle = str_wrap(sub_txt,   90)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position      = "bottom",            
    legend.box           = "horizontal",
    plot.title.position  = "plot",              
    plot.title           = element_text(face = "bold", size = 15, margin = margin(b = 6)),
    plot.subtitle        = element_text(size = 11, margin = margin(b = 8)),
    plot.margin          = margin(t = 14, r = 14, b = 14, l = 14),  
    axis.text.x          = element_text(angle = 0, hjust = 0.5)
  )



print(p2)



ggsave("bayes_rt_timeseries_fixed.png", p2,
       width = 12, height = 6, dpi = 150, units = "in", limitsize = FALSE)



library(bayesplot)
mcmc_trace(fit, pars=c("alpha","beta","lambda","phi","sigma_D","kappa"))
ggsave("mcmc2.png")
mcmc_pairs(as.array(fit), pars=c("alpha","beta","lambda","phi"))
ggsave("mcmc3.png")









suppressPackageStartupMessages({
  library(rstan)
  library(loo)
  library(dplyr)
})



dr <- rstan::extract(fit, pars = c("mu","kappa"))   # default permuted=TRUE
mu    <- dr$mu                  # dim: S x T
kappa <- as.numeric(dr$kappa)   # length S
stopifnot(nrow(mu) == length(kappa))
S  <- nrow(mu)
Tn <- ncol(mu)



ex_mu_pf <- rstan::extract(fit, pars = "mu", permuted = FALSE)  # dim: iter x chains x T
iter   <- dim(ex_mu_pf)[1]
chains <- dim(ex_mu_pf)[2]
stopifnot(S == iter * chains)   # 총 드로우 수 일치 확인
chain_id <- rep(seq_len(chains), each = iter)



y <- dat$Rt
log_lik_mat <- matrix(NA_real_, nrow = S, ncol = Tn)
for (s in 1:S) {
  log_lik_mat[s, ] <- dgamma(y, shape = kappa[s], rate = kappa[s] / mu[s, ], log = TRUE)
}



r_eff  <- loo::relative_eff(exp(log_lik_mat), chain_id = chain_id)
loo_res <- loo::loo(log_lik_mat, r_eff = r_eff)
print(loo_res)



pk_tbl <- loo::pareto_k_table(loo_res)
print(pk_tbl)

png("loo_pareto_k.png", width = 1000, height = 500, res = 150)
plot(loo_res)  

dev.off()



sink("loo_summary.txt"); print(loo_res); print(pk_tbl); sink()
saveRDS(loo_res, "loo_result.rds")




suppressPackageStartupMessages({
  library(dplyr); library(lubridate); library(readr); library(tidyr)
  library(rstan)
})

#


if (exists("fit") && inherits(fit, "stanfit")) {
  saveRDS(fit, "bayes_fit.rds")
  message("Saved: bayes_fit.rds")
}




PATH_RES <- "bayes_rt_results.csv" 

stopifnot(file.exists(PATH_RES))
dat <- read_csv(PATH_RES, show_col_types = FALSE) %>%
  mutate(date = as.Date(date)) %>%
  arrange(date)




scale01 <- function(x){
  if (max(x, na.rm = TRUE) > 1.5) x <- x/100
  pmin(pmax(x, 0), 1)
}
if (!"N" %in% names(dat)) stop("dat에 N 열이 필요합니다(RoOI_adjusted).")
if (!"V" %in% names(dat)) stop("dat에 V 열이 필요합니다(practical vaccination rate).")
if (!"Rt" %in% names(dat)) stop("dat에 Rt 열이 필요합니다.")
if (!"R0_t" %in% names(dat)) stop("dat에 R0_t 열이 필요합니다.")




DATE_START <- as.Date("2020-05-01")
DATE_END   <- as.Date("2022-07-31")

dat_use <- dat %>%
  filter(date >= DATE_START, date <= DATE_END)



monthly_inputs <- dat_use %>%
  mutate(mon = floor_date(date, "month")) %>%
  group_by(mon) %>%
  summarise(
    N_bar   = mean(N, na.rm = TRUE),
    V_bar   = mean(V, na.rm = TRUE),
    Rt_bar  = mean(Rt, na.rm = TRUE),
    R0t_bar = mean(R0_t, na.rm = TRUE),
    R_rho   = 1 - (Rt_bar / R0t_bar),
    .groups = "drop"
  )

saveRDS(monthly_inputs, "monthly_inputs.rds")
message("Saved: monthly_inputs.rds")





if (file.exists("bayes_fit.rds")) {
  fit_rds <- readRDS("bayes_fit.rds")
  dr <- rstan::extract(fit_rds, pars = c("alpha","beta","lambda"))
  a <- as.numeric(dr$alpha); b <- as.numeric(dr$beta)
  lam <- if (!is.null(dr$lambda)) as.numeric(dr$lambda) else rep(0, length(a))
  S <- length(a)
  draws_tbl <- tibble(
    draw = seq_len(S),
    alpha = a, beta = b, lambda = lam
  )
  saveRDS(draws_tbl, "draws_alpha_beta_lambda.rds")
  message("Saved: draws_alpha_beta_lambda.rds")
} else {
  message("주의: bayes_fit.rds를 찾지 못해 드로우 RDS는 건너뜀.")
}




readr::write_csv(monthly_inputs, "monthly_inputs.csv")
message("Saved (csv): monthly_inputs.csv")













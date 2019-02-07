library(brms)
library(tidyverse)
# rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())
# Sys.setenv(LOCAL_CPPFLAGS = '-march=native')


wf <- read_delim("D:/datasets/_requests/whitefly.txt",
                 delim = " ",
                 skip = 10) %>%
  mutate(trt = as.factor(ifelse(trt == 5, 0, trt)))
# imm= number of immature whiteflies
# week= the week in which the measurement was taken
# rep= block (experiment was a randomized complete block design with 
# repeated measures over several weeks)
# trt = the treatment received (I believe 1-4 are increasing levels of 
# pesticide via subirrigation, 5= control and 6= hand watering)
# bindenom = the number of live adult insects placed on the plant
# nlive = the number of adult insects surviving (out of bindenom)
# plantid = an identifier for the plant on which the measurement was taken.

# logit(pbin_ijk) = mu + block_j + trt_i + b*week_k
# logit(pmix_ijk) = c

# wf %>%
#   ggplot(aes(x = nlive/bindenom)) +
#   geom_bar() +
#   facet_wrap(vars(week))

ei_binomial <- custom_family(
  "ei_binomial",
  dpars = c("mu", "po", "pm"),
  links = rep("identity", 3),
  lb = c(NA, 0, 0), ub = c(NA, 1, 1),
  type = "int", vars = c("yo[n]", "ym[n]", "trials[n]")
)

stan_funs <- "
  real ei_binomial_lpmf(int y, real mu, real po, real pm, int yo, int ym, int T) {
    return log(yo*po + ym*pm + (1-po-pm)*exp(binomial_logit_lpmf(y|T, mu)));
  }
  int ei_binomial_rng(real mu, real po, real pm, int T) {
    int which_component = categorical_rng([po, pm, (1-po-pm)]');

    if (which_component == 1) {
      return 0;
    }
    if (which_component == 2) {
      return T;
    }

    return binomial_rng(T, inv_logit(mu));
  }
"

stanvars <- stanvar(scode = stan_funs, block = "functions") +
  stanvar(as.integer(wf$bindenom), name = "trials") +
  stanvar(as.integer(wf$nlive == 0), name = "yo") +
  stanvar(as.integer(wf$nlive == wf$bindenom), name = "ym")

priors <- c(set_prior("normal(0,10)", class = "b"))
a <- list(b = as.array(rep(0, 6)), po = 0.33, pm = 0.33)

eib_fit <- brm(
  nlive ~ trt + week + (1|plantid) + (1|rep), data = wf,
  family = ei_binomial, stanvars = stanvars, prior = priors,
  cores = 3, control = list(adapt_delta = 0.95),
  inits = list(a, a, a, a)
)

#prior_summary(eib_fit)

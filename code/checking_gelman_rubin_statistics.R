rm(list = ls())

percent = function(post){
  summ = as.data.frame(post$BUGSoutput$summary)
  summ2 = summ |>
    filter(Rhat < 1.01)
  100*nrow(summ2) / nrow(summ)
}
load("out/ms_fuxian_1e5.rda")
percent(post.ms)
load("out/ms_fuxian_D47_1e5.rda")
percent(post.ms)
load("out/ms_fuxian_D47_MS_1e5.rda")
percent(post.ms)
load("out/ms_fuxian_D47_MS_ECS_1e5.rda")
percent(post.ms)

load("out/ts_fuxian_1e4.rda")
percent(post.ts)
load("out/ts_fuxian_D47_1e4.rda")
percent(post.ts)
load("out/ts_fuxian_D47_MS_1e4.rda")
percent(post.ts)
load("out/ts_fuxian_D47_MS_ECS_1e4.rda")
percent(post.ts)
load("out/ts_fuxian_D47_MS_ECS_1e5.rda")
percent(post.ts)

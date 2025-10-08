library(tidyverse)
library(asreml)
library(broom.asreml)
dat <- agridat::barrero.maize |> 
    transform(yearf = factor(year))


# m1 ----------------------------------------------------------------------
m1 <- asreml(
  yield ~ loc * yearf,
  data = dat,
  random = ~ gen + rep:env + gen:yearf + gen:loc + gen:env,
  residual = ~ dsum(~ units | env),
  workspace = "500mb"
)
m1$converge

saveRDS(m1, here::here("outputs/01-asreml-barrero-maize.rds"))

# install `broom.asreml` from https://github.com/emitanaka/broom.asreml
# install.packages9("pak")
# pak::pak("emitanaka/broom.asreml")
broom.asreml::tidy(m1, "fixed")
broom.asreml::tidy(m1, "random")



# m2 ----------------------------------------------------------------------
dat <- agridat::barrero.maize |> 
    filter(year %in% 2000:2001) |> 
    mutate(yearf = factor(year)) |> 
    select(-c("yor", "daystoflower", "plantheight", "earheight", "population", "lodged", "moisture", "testweight")) |> 
    na.omit()

m2 <- asreml(
    yield ~ env,
    data = dat,
    random = ~ gen + gen:env,
    residual = ~ dsum(~ units | env),
    workspace = "500mb"
)
m2$converge

saveRDS(m2, here::here("outputs/01-asreml-barrero-maize-m2.rds"))
saveRDS(tidy(m2, "random"), here::here("outputs/01-asreml-barrero-maize-m2-random.rds"))
saveRDS(tidy(m2, "fixed"), here::here("outputs/01-asreml-barrero-maize-m2-fixed.rds"))
saveRDS(tidy(m2, "vcomp"), here::here("outputs/01-asreml-barrero-maize-m2-vcomp.rds"))


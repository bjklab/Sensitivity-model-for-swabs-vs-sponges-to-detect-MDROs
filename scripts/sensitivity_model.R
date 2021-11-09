#' load libraries and set seed
library(tidyverse)
library(tidybayes)
library(brms)
library(gt)

set.seed(16)


#' load data
read_csv("./data/swab_comparison_iceman.csv") %>%
  pivot_wider(id_cols = c(id, day, site, media), names_from = swab_type, values_from = MicrobialSpecies, values_fn = function(x) paste(x, collapse = " ")) %>%
  mutate_at(.vars = c("eswab", "sponge"), .funs = ~ replace(.x, .x == "NA", NA)) %>%
  mutate_at(.vars = c("eswab", "sponge"), .funs = list("positive" = ~ !is.na(.x))) %>%
  rowwise() %>%
  mutate(either_positive = eswab_positive | sponge_positive) %>%
  ungroup() %>%
  identity() -> im
im

im %>%
  gt()



#' observed sensitivity per media/MDRO
im %>%
  group_by(media) %>%
  summarise(sponge_sensitivity = sum(sponge_positive) / sum(either_positive),
            eswab_sensitivity = sum(eswab_positive) / sum(either_positive),
            surfaces = n()) %>%
  ungroup() %>%
  mutate(dif_sensis = sponge_sensitivity - eswab_sensitivity) %>%
  identity() -> im_sensis
im_sensis

im_sensis %>%
  gt()



#' simple model: sponge swab
im %>%
  select(media, contains("positive")) %>%
  filter(complete.cases(.)) %>%
  # 2460 observations
  brm(formula = bf(sponge_positive ~ 1 + either_positive | media),
      data = .,
      family = bernoulli(),
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.99999, max_treedepth = 10),
      backend = "cmdstanr",
      seed = 16,
      file = glue::glue("./models/m_binomial_sponge_media_mixed"),
      file_refit = "on_change") -> m_binomial_sponge_media_mixed

m_binomial_sponge_media_mixed
rstan::check_hmc_diagnostics(m_binomial_sponge_media_mixed$fit)
m_binomial_sponge_media_mixed %>% pp_check()

m_binomial_sponge_media_mixed %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3)


#' fitted: sponge swab
m_binomial_sponge_media_mixed$data %>%
  as_tibble() %>%
  expand(either_positive = TRUE,
         media = unique(media)
  ) %>%
  add_epred_draws(m_binomial_sponge_media_mixed) %>%
  ungroup() %>%
  identity() -> m_binomial_sponge_media_mixed_fitted
m_binomial_sponge_media_mixed_fitted


m_binomial_sponge_media_mixed$data %>%
  as_tibble() %>%
  expand(either_positive = TRUE,
         media = "New Media"
  ) %>%
  add_epred_draws(m_binomial_sponge_media_mixed, allow_new_levels = TRUE) %>%
  ungroup() %>%
  identity() -> m_binomial_sponge_media_mixed_unknown
m_binomial_sponge_media_mixed_unknown


m_binomial_sponge_media_mixed_fitted %>%
  bind_rows(m_binomial_sponge_media_mixed_unknown) %>%
  mutate(media = factor(media, levels = c("CDIFF", "CRE", "CRPA", "ESBL", "MRSA", "VRE", "New Media"))) %>%
  identity() -> m_binomial_sponge_media_mixed_plot
  


m_binomial_sponge_media_mixed_plot %>%
  ggplot(data = ., aes(y = media, x = .epred, fill = media)) +
  tidybayes::stat_halfeye(.width = 0.95) +
  scale_fill_viridis_d(option = "turbo", begin = 0.2, end = 0.8) +
  scale_x_continuous(limits = c(0,1)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title.x = ggtext::element_markdown(color = "black"),
        axis.title.y = ggtext::element_markdown(color = "black"),
        legend.position = "bottom",
        plot.title.position = "plot") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Expected Sponge Swab Sensitivity<br>(posterior median and 95% credible interval)",
       y = "",
       fill = "Media",
       title = "Sponge Swab Sensitivity for MDRO Detection in the Healthcare Environment",
       subtitle = "Estimated from 2460 Hospital Surface Cultures") -> p_binomial_sponge_media_mixed_fitted
p_binomial_sponge_media_mixed_fitted


p_binomial_sponge_media_mixed_fitted %>%
  ggsave(plot = ., filename = "./figs/p_binomial_sponge_media_mixed_fitted.pdf", height = 6, width = 8, units = "in")
p_binomial_sponge_media_mixed_fitted %>%
  ggsave(plot = ., filename = "./figs/p_binomial_sponge_media_mixed_fitted.svg", height = 6, width = 8, units = "in")
p_binomial_sponge_media_mixed_fitted %>%
  ggsave(plot = ., filename = "./figs/p_binomial_sponge_media_mixed_fitted.png", height = 6, width = 8, units = "in", dpi = 600)





#' simple model: eswab
im %>%
  select(media, contains("positive")) %>%
  filter(complete.cases(.)) %>%
  # 2460 observations
  brm(formula = bf(eswab_positive ~ 1 + either_positive | media),
      data = .,
      family = bernoulli(),
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.99999, max_treedepth = 10),
      backend = "cmdstanr",
      seed = 16,
      file = glue::glue("./models/m_binomial_eswab_media_mixed"),
      file_refit = "on_change") -> m_binomial_eswab_media_mixed

m_binomial_eswab_media_mixed
rstan::check_hmc_diagnostics(m_binomial_eswab_media_mixed$fit)
m_binomial_eswab_media_mixed %>% pp_check()

m_binomial_eswab_media_mixed %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3)


#' fitted: eswab
m_binomial_eswab_media_mixed$data %>%
  as_tibble() %>%
  expand(either_positive = TRUE,
         media = unique(media)
  ) %>%
  add_epred_draws(m_binomial_eswab_media_mixed) %>%
  ungroup() %>%
  identity() -> m_binomial_eswab_media_mixed_fitted
m_binomial_eswab_media_mixed_fitted


m_binomial_eswab_media_mixed$data %>%
  as_tibble() %>%
  expand(either_positive = TRUE,
         media = "New Media"
  ) %>%
  add_epred_draws(m_binomial_eswab_media_mixed, allow_new_levels = TRUE) %>%
  ungroup() %>%
  identity() -> m_binomial_eswab_media_mixed_unknown
m_binomial_eswab_media_mixed_unknown


m_binomial_eswab_media_mixed_fitted %>%
  bind_rows(m_binomial_eswab_media_mixed_unknown) %>%
  mutate(media = factor(media, levels = c("CDIFF", "CRE", "CRPA", "ESBL", "MRSA", "VRE", "New Media"))) %>%
  identity() -> m_binomial_eswab_media_mixed_plot



m_binomial_eswab_media_mixed_plot %>%
  ggplot(data = ., aes(y = media, x = .epred, fill = media)) +
  tidybayes::stat_halfeye(.width = 0.95) +
  scale_fill_viridis_d(option = "turbo", begin = 0.2, end = 0.8) +
  scale_x_continuous(limits = c(0,1)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title.x = ggtext::element_markdown(color = "black"),
        axis.title.y = ggtext::element_markdown(color = "black"),
        legend.position = "bottom",
        plot.title.position = "plot") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Expected E-Swab Sensitivity<br>(posterior median and 95% credible interval)",
       y = "",
       fill = "Media",
       title = "E-Swab Sensitivity for MDRO Detection in the Healthcare Environment",
       subtitle = "Estimated from 2460 Hospital Surface Cultures") -> p_binomial_eswab_media_mixed_fitted
p_binomial_eswab_media_mixed_fitted


p_binomial_eswab_media_mixed_fitted %>%
  ggsave(plot = ., filename = "./figs/p_binomial_eswab_media_mixed_fitted.pdf", height = 6, width = 8, units = "in")
p_binomial_eswab_media_mixed_fitted %>%
  ggsave(plot = ., filename = "./figs/p_binomial_eswab_media_mixed_fitted.svg", height = 6, width = 8, units = "in")
p_binomial_eswab_media_mixed_fitted %>%
  ggsave(plot = ., filename = "./figs/p_binomial_eswab_media_mixed_fitted.png", height = 6, width = 8, units = "in", dpi = 600)



#' combined plots
m_binomial_eswab_media_mixed_plot %>%
  mutate(swab = "E-Swab") %>%
  bind_rows(mutate(m_binomial_sponge_media_mixed_plot, swab = "Sponge Swab")) %>%
  filter(media != "New Media") %>%
  ggplot(data = ., aes(y = swab, x = .epred, fill = media)) +
  facet_wrap(facets = ~ media) +
  tidybayes::stat_halfeye(.width = 0.95) +
  scale_fill_viridis_d(option = "turbo", begin = 0.2, end = 0.8) +
  scale_x_continuous(limits = c(0,1)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title.x = ggtext::element_markdown(color = "black"),
        axis.title.y = ggtext::element_markdown(color = "black"),
        legend.position = "bottom",
        plot.title.position = "plot") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Expected Sensitivity<br>(posterior median and 95% credible interval)",
       y = "",
       fill = "Media",
       title = "Sensitivity for MDRO Detection in the Healthcare Environment",
       subtitle = "Estimated from 2460 Hospital Surface Cultures") -> p_binomial_combined_swabs_media_mixed_fitted
p_binomial_combined_swabs_media_mixed_fitted



m_binomial_eswab_media_mixed_plot %>%
  mutate(swab = "E-Swab") %>%
  bind_rows(mutate(m_binomial_sponge_media_mixed_plot, swab = "Sponge Swab")) %>%
  filter(media == "New Media") %>%
  ggplot(data = ., aes(y = swab, x = .epred)) +
  #facet_wrap(facets = ~ media) +
  tidybayes::stat_halfeye(.width = 0.95, fill = "grey") +
  #scale_fill_viridis_d(option = "turbo", begin = 0.2, end = 0.8) +
  scale_x_continuous(limits = c(0,1)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title.x = ggtext::element_markdown(color = "black"),
        axis.title.y = ggtext::element_markdown(color = "black"),
        legend.position = "bottom",
        plot.title.position = "plot") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Expected Sensitivity Over All Media<br>(posterior median and 95% credible interval)",
       y = ""#,
       #fill = "Media",
       #title = "Sensitivity for MDRO Detection in the Healthcare Environment",
       #subtitle = "Estimated from 2460 Hospital Surface Cultures"
       ) -> p_binomial_combined_swabs_newmedia_mixed_fitted
p_binomial_combined_swabs_newmedia_mixed_fitted



library(patchwork)


((p_binomial_combined_swabs_media_mixed_fitted + theme(legend.position = "none")) /
  p_binomial_combined_swabs_newmedia_mixed_fitted) %>%
  identity() -> p_swab_sensitivity_models_combined
p_swab_sensitivity_models_combined




p_swab_sensitivity_models_combined %>%
  ggsave(plot = ., filename = "./figs/p_swab_sensitivity_models_combined.pdf", height = 6, width = 8, units = "in")
p_swab_sensitivity_models_combined %>%
  ggsave(plot = ., filename = "./figs/p_swab_sensitivity_models_combined.svg", height = 6, width = 8, units = "in")
p_swab_sensitivity_models_combined %>%
  ggsave(plot = ., filename = "./figs/p_swab_sensitivity_models_combined.png", height = 6, width = 8, units = "in", dpi = 600)



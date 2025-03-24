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



#' ###############################
#' separate sponge & e-swab models
#' ###############################

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
  mutate(media = factor(media, levels = c("CDIFF", "CRE", "KPC", "ESBL", "MRSA", "VRE", "New Media"))) %>%
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
  mutate(media = factor(media, levels = c("CDIFF", "CRE", "KPC", "ESBL", "MRSA", "VRE", "New Media"))) %>%
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
  labs(x = "Expected Flocked Swab Sensitivity<br>(posterior median and 95% credible interval)",
       y = "",
       fill = "Media",
       title = "Flocked Swab Sensitivity for MDRO Detection in the Healthcare Environment",
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
  mutate(swab = "Flocked Swab") %>%
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
  mutate(swab = "Flocked Swab") %>%
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
  p_binomial_combined_swabs_newmedia_mixed_fitted) +
  plot_layout(heights = c(2,1)) %>%
  identity() -> p_swab_sensitivity_models_combined
p_swab_sensitivity_models_combined




p_swab_sensitivity_models_combined %>%
  ggsave(plot = ., filename = "./figs/p_swab_sensitivity_models_combined.pdf", height = 6, width = 8, units = "in")
p_swab_sensitivity_models_combined %>%
  ggsave(plot = ., filename = "./figs/p_swab_sensitivity_models_combined.svg", height = 6, width = 8, units = "in")
p_swab_sensitivity_models_combined %>%
  ggsave(plot = ., filename = "./figs/p_swab_sensitivity_models_combined.png", height = 6, width = 8, units = "in", dpi = 600)






#' ###############################
#' combined sponge & e-swab model
#' ###############################

#' fully crossed random effects
im %>%
  select(media, contains("positive")) %>%
  filter(complete.cases(.)) %>%
  pivot_longer(cols = c(-media, -either_positive), names_to = "swab", values_to = "detected") %>%
  mutate(swab = gsub("_positive","",swab)) %>%
  # 2460 observations
  brm(formula = bf(detected ~ 1 + either_positive + (1 + either_positive | media) + (1 + either_positive | swab)),
      data = .,
      family = bernoulli(),
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.99999, max_treedepth = 10),
      backend = "cmdstanr",
      seed = 16,
      file = glue::glue("./models/m_binomial_swab_media_crossmixed"),
      file_refit = "on_change") -> m_binomial_swab_media_crossmixed

m_binomial_swab_media_crossmixed
rstan::check_hmc_diagnostics(m_binomial_swab_media_crossmixed$fit)
m_binomial_swab_media_crossmixed %>% pp_check()

m_binomial_swab_media_crossmixed %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3)


#' fitted: cross mixed model
m_binomial_swab_media_crossmixed$data %>%
  as_tibble() %>%
  expand(either_positive = TRUE,
         media = unique(media),
         swab = unique(swab)
  ) %>%
  add_epred_draws(m_binomial_swab_media_crossmixed) %>%
  ungroup() %>%
  identity() -> m_binomial_swab_media_crossmixed_fitted
m_binomial_swab_media_crossmixed_fitted



m_binomial_swab_media_crossmixed$data %>%
  as_tibble() %>%
  expand(either_positive = TRUE,
         media = "New Media",
         swab = unique(swab)
  ) %>%
  add_epred_draws(m_binomial_swab_media_crossmixed, allow_new_levels = TRUE, sample_new_levels = "uncertainty", re_formula = NULL) %>%
  ungroup() %>%
  identity() -> m_binomial_swab_media_crossmixed_unknown
m_binomial_swab_media_crossmixed_unknown


m_binomial_swab_media_crossmixed_fitted %>%
  bind_rows(m_binomial_swab_media_crossmixed_unknown) %>%
  mutate(media = factor(media, levels = c("CDIFF", "CRE", "KPC", "ESBL", "MRSA", "VRE", "New Media"))) %>%
  mutate(swab = case_when(swab == "sponge" ~ "Sponge<br>Stick",
                          swab == "eswab" ~ "Flocked<br>Swab")) %>%
  identity() -> m_binomial_swab_media_crossmixed_plot



m_binomial_swab_media_crossmixed_plot %>%
  filter(media != "New Media") %>%
  ggplot(data = ., aes(y = swab, x = .epred, fill = media)) +
  tidybayes::stat_halfeye(.width = 0.95) +
  facet_wrap(facets = ~ media) +
  scale_fill_viridis_d(option = "turbo", begin = 0.2, end = 0.8) +
  scale_x_continuous(limits = c(0,1)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        axis.title.x = ggtext::element_markdown(color = "black"),
        axis.title.y = ggtext::element_markdown(color = "black"),
        legend.position = "bottom",
        plot.title.position = "plot") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Expected Sensitivity<br>(posterior median and 95% credible interval)",
       y = "",
       fill = "Media"#,
       #title = "Sensitivity for MDRO Detection in the Healthcare Environment",
       #subtitle = "Estimated from 2460 Hospital Surface Cultures"
       ) -> p_binomial_swab_media_crossmixed_fitted
p_binomial_swab_media_crossmixed_fitted



m_binomial_swab_media_crossmixed_plot %>%
  filter(media == "New Media") %>%
  ggplot(data = ., aes(y = swab, x = .epred), fill = "grey") +
  tidybayes::stat_halfeye(.width = 0.95) +
  scale_fill_viridis_d(option = "turbo", begin = 0.2, end = 0.8) +
  scale_x_continuous(limits = c(0,1)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        axis.title.x = ggtext::element_markdown(color = "black"),
        axis.title.y = ggtext::element_markdown(color = "black"),
        legend.position = "bottom",
        plot.title.position = "plot") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Expected Sensitivity Over Any Media<br>(posterior median and 95% credible interval)",
       y = ""#,
       #fill = "Media",
       #title = "Sensitivity for MDRO Detection in the Healthcare Environment",
       #subtitle = "Estimated from 2460 Hospital Surface Cultures"
       ) -> p_binomial_swab_newmedia_crossmixed_fitted
p_binomial_swab_newmedia_crossmixed_fitted


library(patchwork)


((p_binomial_swab_media_crossmixed_fitted + theme(legend.position = "none")) /
    p_binomial_swab_newmedia_crossmixed_fitted) +
  plot_layout(heights = c(2,1)) |> 
  #plot_annotation(tag_levels = "A", title = NULL, subtitle = NULL) %>%
  patchwork::plot_annotation(tag_levels = "A", title = NULL, subtitle = NULL, caption = NULL) |> 
  identity() -> p_swab_sensitivity_models_together
p_swab_sensitivity_models_together




p_swab_sensitivity_models_together %>%
  ggsave(plot = ., filename = "./figs/p_swab_sensitivity_models_together.pdf", height = 6, width = 8, units = "in")
p_swab_sensitivity_models_together %>%
  ggsave(plot = ., filename = "./figs/p_swab_sensitivity_models_together.svg", height = 6, width = 8, units = "in")
p_swab_sensitivity_models_together %>%
  ggsave(plot = ., filename = "./figs/p_swab_sensitivity_models_together.png", height = 6, width = 8, units = "in", dpi = 600)




#' model contrasts for Results section
#' 
#' 

m_binomial_swab_media_crossmixed_fitted |> 
  #count(swab)
  group_by(media) |> 
  pivot_wider(id_cols = c(media, .draw), names_from = swab, values_from = .epred) |> 
  mutate(contrast_s_less_e = sponge - eswab) |>
  ungroup() |> 
  identity() -> m_binomial_swab_media_crossmixed_fitted_contrasts
m_binomial_swab_media_crossmixed_fitted_contrasts

m_binomial_swab_media_crossmixed_plot |> 
  #count(swab)
  group_by(media) |> 
  pivot_wider(id_cols = c(media, .draw), names_from = swab, values_from = .epred) |> 
  mutate(contrast_s_less_e = `Sponge<br>Stick` - `Flocked<br>Swab`) |>
  ungroup() |> 
  identity() -> m_binomial_swab_media_crossmixed_fitted_contrasts
m_binomial_swab_media_crossmixed_fitted_contrasts

m_binomial_swab_media_crossmixed_fitted_contrasts |> 
  group_by(media) |> 
  tidybayes::median_qi(contrast_s_less_e, .width = 0.95) |> 
  mutate(text = glue::glue("for {media}, sponge stick was {round(contrast_s_less_e,3)*100}% (95%CrI {round(.lower,3)*100}% to {round(.upper,3)*100}%) more sensitive than flocked swab;")) |> 
  identity() -> m_binomial_swab_media_crossmixed_fitted_contrasts_summary
m_binomial_swab_media_crossmixed_fitted_contrasts_summary |> 
  pull(text)

m_binomial_swab_media_crossmixed_fitted_contrasts_summary |> 
  write_csv("tabs/m_binomial_swab_media_crossmixed_fitted_contrasts_summary.csv")



# 2a. Analysis (0.5 km buffer; main text) --------------
library(tidyverse)
library(brms) # 2.17.0

# Data ------------------
gdata <- read.csv("SchmidtGarroway_PNAS_2022_analysisdata.csv", header = TRUE) %>% 
  mutate(species = as.factor(species), 
         class = as.factor(class), 
         author = as.factor(author))

# Data subsets ----------------
## 0.5 km --------------------

a10gd <- gdata %>% 
  drop_na(total.5) %>% # drop sites with no census data within 0.5km buffer
  # Compute & scale variables
  mutate(scalehfi = scale(hfi),
         pwhite.5 = scale(white.5/total.5),
         scaleAR = scale(allelic_richness), 
         scaleGD = scale(gene_diversity))

a10fst <- a10gd %>%
  drop_na(pop_fst) %>% # drop sites with no FST estimate
  mutate(scalehfi = scale(hfi), 
         pwhite.5 = scale(white.5/total.5),
         scalefst = scale(pop_fst))

a10ne <- a10gd %>%
  drop_na(Ne) %>% # drop sites with no Ne estimate
  mutate(scalehfi = scale(hfi), 
         pwhite.5 = scale(white.5/total.5),
         scalene = scale(log(Ne)))

# Model prep ------------------
brm_fun <- function(model, iter = 3000, prior = prior, data){
  brm(model, cores = 4, chains = 4, iter = iter, 
      prior = prior,
      control = list(adapt_delta = 0.999, max_treedepth = 15),
      data = data)
}

## Priors
### Normal(0.5, 0.25)
gan_pw <- prior(normal(0.5, 0.25), class = b, coef = pwhite.5)
gan_hfi <- prior(normal(-0.5, 0.25), class = b, coef = scalehfi)

gan_b <- c(prior(normal(0.5, 0.25), class = b, coef = pwhite.5),
           prior(normal(-0.5, 0.25), class = b, coef = scalehfi))

f_pw <- prior(normal(-0.5, 0.25), class = b, coef = pwhite.5)
f_hfi <- prior(normal(0.5, 0.25), class = b, coef = scalehfi)

f_b <- c(prior(normal(-0.5, 0.25), class = b, coef = pwhite.5),
         prior(normal(0.5, 0.25), class = b, coef = scalehfi))

### Student(0.5, 0.25), 3 DF
gan_pwS <- prior(student_t(3, 0.5, 0.25), class = b, coef = pwhite.5)
f_pwS <- prior(student_t(3, -0.5, 0.25), class = b, coef = pwhite.5)

### Normal(0, 1)
N01 <- prior(normal(0, 1), class = b, coef = pwhite.5)

# Models ---------------
## *Allelic richness --------------
# pwhite
# Model formula
for_arpw <- bf(scaleAR ~ pwhite.5 + (pwhite.5|species) + (1|class/species))

mod_arpw <- brm_fun(for_arpw, prior = gan_pw, data = a10gd)

# Plot model
plot(mod_arpw)
# Posterior predictive check
brms::pp_check(mod_arpw, type = "dens_overlay", nsamples = 100)
# Model summary
summary(mod_arpw)

# HFI
for_arhfi <- bf(scaleAR ~ scale_hfi + (scale_hfi|species) + (1|class/species))
mod_arhfi <- brm_fun(for_arhfi, prior = gan_hfi, data = a10gd)

# pwhite + HFI
for_arpwhfi <- bf(scaleAR ~ pwhite.5 + scalehfi + 
                    (pwhite.5 + scalehfi|species) + (1|class/species))
mod_arpwhfi <- brm_fun(for_arpwhfi, prior = gan_b, data = a10gd)

# Student prior:
mod_arpwS <- brm_fun(for_arpw, iter = 5000, prior = gan_pwS, data = a10gd)

# Normal prior:
mod_arpwN <- brm_fun(for_arpw, prior = N01, data = a10gd)

##### *Gene diversity --------------
# pwhite
for_gdpw <- bf(scaleGD ~ pwhite.5 + (pwhite.5|species) + (1|class/species))

mod_gdpw <- brm_fun(for_gdpw, prior = gan_pw, data = a10gd)

# HFI
for_gdhfi <- bf(scaleGD ~ scale_hfi + (scale_hfi|species) + (1|class/species))
mod_gdhfi <- brm_fun(for_gdhfi, prior = gan_hfi, data = a10gd)

# pwhite + HFI
for_gdpwhfi <- bf(scaleGD ~ pwhite.5 + scalehfi + 
                    (pwhite.5 + scalehfi|species) + (1|class/species))
mod_gdpwhfi <- brm_fun(for_gdpwhfi, prior = gan_b, data = a10gd)

# Student prior:
mod_gdpwS <- brm_fun(for_gdpw, prior = gan_pwS, data = a10gd)

# Normal prior:
mod_gdpwN <- brm_fun(for_gdpw, prior = N01, data = a10gd)

##### *Ne --------------
# pwhite
for_nepw <- bf(scalene ~ pwhite.5 + (pwhite.5|species) + (1|class/species))

mod_nepw <- brm_fun(for_nepw, prior = gan_pw, data = a10ne)

# HFI
for_nehfi <- bf(scalene ~ scale_hfi + (scale_hfi|species) + (1|class/species))
mod_nehfi <- brm_fun(for_nehfi, prior = gan_hfi, data = a10ne)

# pwhite + HFI
for_nepwhfi <- bf(scalene ~ pwhite.5 + scalehfi + 
                    (pwhite.5 + scalehfi|species) + (1|class/species))
mod_nepwhfi <- brm_fun(for_nepwhfi, prior = gan_b, data = a10ne)

# Student prior:
mod_nepwS <- brm_fun(for_nepw, prior = gan_pwS, data = a10ne)

# Normal prior:
mod_nepwN <- brm_fun(for_nepw, prior = N01, data = a10ne)

##### *FST --------------
# pwhite
for_fstpw <- bf(scalefst ~ pwhite.5 + (pwhite.5|species) + (1|class/species))

mod_fstpw <- brm_fun(for_fstpw, prior = f_pw, data = a10fst)

# HFI
for_fsthfi <- bf(scalefst ~ scale_hfi + (scale_hfi|species) + (1|class/species))
mod_fsthfi <- brm_fun(for_fsthfi, prior = f_hfi, data = a10fst)

# pwhite + HFI
for_fstpwhfi <- bf(scalefst ~ pwhite.5 + scalehfi + 
                     (pwhite.5 + scalehfi|species) + (1|class/species))
mod_fstpwhfi <- brm_fun(for_fstpwhfi, prior = f_b, data = a10fst)

# Student prior:
mod_fstpwS <- brm_fun(for_fstpw, iter = 5000, prior = f_pwS, data = a10fst)

# Normal prior:
mod_fstpwN <- brm_fun(for_fstpw, prior = N01, data = a10fst)

# Testing residual spatial autocorrelation --------
library(sf)
library(adespatial)
library(spdep)

# Gene diversity & allelic richness:
xyg <- dplyr::select(a10gd, lon, lat)
xysfg <- st_as_sf(xyg, coords = c("lon", "lat"), crs=4326)

## Generate neighborhood matrix
nbg <- dnearneigh(xysfg, 0, 450) 

# Moran test:
GDres <- residuals(mod_gdpw)[,1]
moran.test(GDres, nb2listw(nbg))

ARres <- residuals(mod_arpw)[,1]
moran.test(ARres, nb2listw(nbg))

# FST
xyf <- dplyr::select(a10fst, lon, lat)
xysff <- st_as_sf(xyf, coords = c("lon", "lat"), crs=4326)

## Generate neighborhood matrix
nbf <- dnearneigh(xysff, 0, 450)

# Moran test:
FSTres <- residuals(mod_fstpw)[,1]
moran.test(FSTres, nb2listw(nbf))

# Ne
xyn <- dplyr::select(a10ne, lon, lat)
xysfn <- st_as_sf(xyn, coords = c("lon", "lat"), crs=4326)

## Generate neighborhood matrix
nbn <- dnearneigh(xysfn, 0, 500)

# Moran test:
Neres <- residuals(mod_nepw)[,1]
moran.test(Neres, nb2listw(nbn))

# R2 --------
# Marginal & conditional R2
library(performance)
lapply(list(mod_arpw, mod_gdpw, mod_nepw, mod_fstpw), r2_bayes)

lapply(list(mod_arhfi, mod_gdhfi, mod_nehfi, mod_fsthfi), r2_bayes)

lapply(list(mod_arpwhfi, mod_gdpwhfi, 
            mod_nepwhfi, mod_fstpwhfi), r2_bayes)

# Adjusted R2
lapply(list(mod_arpwhfi, mod_gdpwhfi, 
            mod_nepwhfi, mod_fstpwhfi), loo_R2)

# Plot (Orchard plot) ----------
library(tidybayes)
library(extrafont)

summaryframe <- function(model, response){
  data.frame(Variable = rownames(summary(model)$fixed)[2],
             Coefficient = summary(model)$fixed[2, 1],
             CIlo95 = summary(model)$fixed[2,3],
             CIup95 = summary(model)$fixed[2,4],
             CIlo90 = posterior_interval(model, variable = "b_pwhite.5", prob = 0.90)[,1],
             CIup90 = posterior_interval(model, variable = "b_pwhite.5", prob = 0.90)[,2],
             Response_var = response)
}

# Model summaries
gdmod.df <- summaryframe(mod_gdpw, response = "gene diversity")
armod.df <- summaryframe(mod_arpw, response = "allelic richness")
nemod.df <- summaryframe(mod_nepw, response = "effective population size")
fstmod.df <- summaryframe(mod_fstpw, response = "FST")

# Combine:
allModelFrame <- data.frame(rbind(gdmod.df, armod.df, fstmod.df, nemod.df))

# Reorder factor levels:
allModelFrame$Response_var <- factor(allModelFrame$Response_var, 
                                     levels = c("FST",
                                                "allelic richness",
                                                "gene diversity",
                                                "effective population size"))
## Species effects
# Taxonomy:
taxoinfo <- gdata %>% 
  select(species, class) %>% 
  distinct()

# Summarise species effects from models & put in alphabetical order
spp_plotdat <- function(model){
  model %>% 
    spread_draws(b_pwhite.5, r_species[species, term]) %>% # get posterior draws
    mutate(species_mean = b_pwhite.5 + r_species) %>% # add offset to mean effect
    filter(term == "pwhite.5") %>% # keep only slope terms
    full_join(., taxoinfo, by = "species") %>% # join with taxonomic data
    drop_na() %>% 
    mutate(species = gsub("_", " ", species)) %>%
    group_by(class) %>% 
    arrange(species, .by_group = TRUE) %>%
    mutate(species = factor(species, levels = unique(.$species))) %>% 
    ungroup()
}

# Gene diversity:
gdplotdat <- spp_plotdat(mod_gdpw) %>% 
  group_by(class, species) %>% 
  summarise(sp_effect = mean(species_mean)) %>% 
  ungroup() %>% 
  mutate(Response_var = "gene diversity")

# Get species' sample sizes
sitesgdar <- a10gd %>% 
  group_by(species) %>% 
  summarise(sites = n()) %>% 
  ungroup() %>% 
  mutate(species = gsub('_', ' ', species))

# merge:
gdplotdat <- merge(gdplotdat, sitesgdar, by = "species")

# Allelic richness:
arplotdat <- spp_plotdat(mod_arpw) %>% 
  group_by(class, species) %>% 
  summarise(sp_effect = mean(species_mean)) %>% 
  ungroup() %>% 
  mutate(Response_var = "allelic richness")

# Merge with species' sample sizes (same as gene diversity)
arplotdat <- merge(arplotdat, sitesgdar, by = "species")

# Ne:
neplotdat <- spp_plotdat(mod_nepw) %>% 
  group_by(class, species) %>% 
  summarise(sp_effect = mean(species_mean)) %>% 
  ungroup() %>% 
  mutate(Response_var = "effective population size")

# Get species' sample sizes
sitesne <- a10ne %>% 
  group_by(species) %>% 
  summarise(sites = n()) %>% 
  ungroup() %>% 
  mutate(species = gsub('_', ' ', species))

# Merge
neplotdat <- merge(neplotdat, sitesne, by = "species")

# FST
fstplotdat <- spp_plotdat(mod_fstpw) %>% 
  group_by(class, species) %>% 
  summarise(sp_effect = mean(species_mean)) %>% 
  ungroup() %>% 
  mutate(Response_var = "FST")

# Get species' sample sizes
sitesfst <- a10fst %>% 
  group_by(species) %>% 
  summarise(sites = n()) %>% 
  ungroup() %>% 
  mutate(species = gsub('_', ' ', species))

# merge
fstplotdat <- merge(fstplotdat, sitesfst, by = "species")

# Combine all:
speciesfx <- rbind(gdplotdat, arplotdat, fstplotdat, neplotdat)

# Reorder factor levels:
speciesfx$Response_var <- factor(speciesfx$Response_var, 
                                 levels = c("FST",
                                            "allelic richness",
                                            "gene diversity",
                                            "effective population size"))
### Plot ------------
pal <- c("#274659", "#CAAE10", "#F2790F", "#B93102")

ggplot(allModelFrame) + 
  geom_hline(yintercept=seq(-1.00, 0.5, 0.25),  # x axis lines (will be flipped)
             lwd=0.25, colour="#DEDEDE") +
  geom_hline(yintercept = 0, colour = "black", lty = 2) +
  geom_jitter(data = speciesfx, aes(x = Response_var, y = sp_effect, size = sites, color = class), 
              width = 0.3, alpha = 0.6) +
  scale_color_manual(values = pal) +
  geom_linerange(aes(x = Response_var, ymin = CIlo90,
                     ymax = CIup90),
                 lwd = 2.5, position = position_dodge(width = 1), color = "#1F2749") +
  geom_pointrange(aes(x = Response_var, y = Coefficient, ymin = CIlo95,
                      ymax = CIup95),
                  lwd = 1, position = position_dodge(width = 1),
                  shape = 21, fill = "white", stroke = 2, color = "#1F2749") +
  geom_vline(xintercept=seq(1.5, length(unique(allModelFrame$Response_var))-0.5, 1), # y axis lines between groups
             lwd=0.25, colour="#DEDEDE") +
  scale_x_discrete(limits = levels(allModelFrame$Response_var),
                   labels = c(expression(F[ST]), expression("allelic richness"), 
                              expression("gene diversity"), 
                              expression("effective population size"))) +
  coord_flip() + 
  theme_classic(base_size = 14) +
  theme(axis.ticks.y = element_blank()) +
  labs(x= "", y = "model coefficients", title = "") +
  theme(text=element_text(family="Roboto Medium"))

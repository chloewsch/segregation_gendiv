# 2b. Analysis (1 & 5 km buffer; SI) --------------
library(brms)
library(tidyverse)
library(tidybayes)
library(performance)

# Data ------------------
gdata <- read.csv("SchmidtGarroway_PNAS_2022_analysisdata.csv", header = TRUE) %>% 
  mutate(species = as.factor(species), 
         class = as.factor(class), 
         author = as.factor(author))

# Data subsets ----------------
## 1 km --------------------
a10gd1 <- gdata %>% 
  drop_na(total1) %>% # drop sites with no census data within 1km buffer
  # Compute & scale variables
  mutate(scalehfi = scale(hfi),
         pwhite1 = scale(white1/total1),
         scaleAR = scale(allelic_richness),
         scaleGD = scale(gene_diversity))

a10fst1 <- a10gd1 %>%
  drop_na(pop_fst) %>% 
  mutate(scalehfi = scale(hfi), 
         pwhite1 = scale(white1/total1),
         scalefst = scale(pop_fst))

a10ne1 <- a10gd1 %>% 
  drop_na(Ne) %>% 
  mutate(scalehfi = scale(hfi),
         pwhite1 = scale(white1/total1),
         scalene = scale(log(Ne)))

## 5 km --------------------
a10gd5 <- gdata %>% 
  drop_na(total5) %>% 
  mutate(scalehfi = scale(hfi),
         pwhite5 = scale(white5/total5),
         scaleAR = scale(allelic_richness),
         scaleGD = scale(gene_diversity))

a10fst5 <- a10gd5 %>%
  drop_na(pop_fst) %>% 
  mutate(scalehfi = scale(hfi), 
         pwhite5 = scale(white5/total5),
         scalefst = scale(pop_fst))

a10ne5 <- a10gd5 %>% 
  drop_na(Ne) %>% 
  mutate(scalehfi = scale(hfi),
         pwhite5 = scale(white5/total5),
         scalene = scale(log(Ne)))

# Model prep ----------
brm_fun <- function(model, iter = 3000, prior = prior, adelta = 0.999, data){
  brm(model, cores = 4, chains = 4, iter = iter, 
      prior = prior,
      control = list(adapt_delta = adelta, max_treedepth = 15),
      data = data)
}

## Priors
gan_pw1 <- prior(normal(0.5, 0.25), class = b, coef = pwhite1)
gan_pw5 <- prior(normal(0.5, 0.25), class = b, coef = pwhite5)

f_pw1 <- prior(normal(-0.5, 0.25), class = b, coef = pwhite1)
f_pw5 <- prior(normal(-0.5, 0.25), class = b, coef = pwhite5)



# Models ------------------------
## *Allelic richness ----------
# 1 km
for_arpw1 <- bf(scaleAR ~ pwhite1 + (pwhite1|species) + (1|class/species))

mod_arpw1 <- brm_fun(for_arpw1, prior = gan_pw1, data = a10gd1)

# Plot model
plot(mod_arpw1)
# Posterior predictive check
brms::pp_check(mod_arpw1, type = "dens_overlay", nsamples = 100)
# Model summary
summary(mod_arpw1)

# 5 km
for_arpw5 <- bf(scaleAR ~ pwhite5 + (pwhite5|species) + (1|class/species))

mod_arpw5 <- brm_fun(for_arpw5, iter = 5000, prior = gan_pw5, data = a10gd5)

## *Gene diversity ------------
# 1 km
for_gdpw1 <- bf(scaleGD ~ pwhite1 + (pwhite1|species) + (1|class/species))

mod_gdpw1 <- brm_fun(for_gdpw1, iter = 4000, prior = gan_pw1, data = a10gd1)

# 5 km
for_gdpw5 <- bf(scaleGD ~ pwhite5 + (pwhite5|species) + (1|class/species))

mod_gdpw5 <- brm_fun(for_gdpw5, iter = 4000, prior = gan_pw5, data = a10gd5)

## *Ne ----------
# 1 km
for_nepw1 <- bf(scalene ~ pwhite1 + (pwhite1|species) + (1|class/species))

mod_nepw1 <- brm_fun(for_nepw1, prior = gan_pw1, data = a10ne1)


# 5 km
for_nepw5 <- bf(scalene ~ pwhite5 + (pwhite5|species) + (1|class/species))

mod_nepw5 <- brm_fun(for_nepw5, prior = gan_pw5, data = a10ne5)

## *FST -----------
# pwhite
# 1 km
for_fstpw1 <- bf(scalefst ~ pwhite1 + (pwhite1|species) + (1|class/species))

mod_fstpw1 <- brm_fun(for_fstpw1, iter = 4500, adelta = 0.9999, prior = f_pw1, data = a10fst1)

# 5 km
for_fstpw5 <- bf(scalefst ~ pwhite5 + (pwhite5|species) + (1|class/species))

mod_fstpw5 <- brm_fun(for_fstpw5, iter = 4500, adelta = 0.9999, prior = f_pw5, data = a10fst5)

## Testing residual spatial autocorrelation ---------------
## Gene diversity & allelic richness
xyg1 <- dplyr::select(a10gd1, lon, lat)
xysfg1 <- st_as_sf(xyg1, coords = c("lon", "lat"), crs=4326)

xyg5 <- dplyr::select(a10gd5, lon, lat)
xysfg5 <- st_as_sf(xyg5, coords = c("lon", "lat"), crs=4326)

### Generate neighborhood matrix
nb1g <- dnearneigh(xysfg1, 0, 450) 
nb5g <- dnearneigh(xysfg5, 0, 450) 

### Moran test
GDres1 <- residuals(mod_gdpw1)[,1]
GDres5 <- residuals(mod_gdpw5)[,1]
moran.test(GDres1, nb2listw(nb1g)) 
moran.test(GDres5, nb2listw(nb5g))

ARres1 <- residuals(mod_arpw1)[,1]
ARres5 <- residuals(mod_arpw5)[,1]
moran.test(ARres1, nb2listw(nb1g))
moran.test(ARres5, nb2listw(nb5g))

## FST
xy1f <- dplyr::select(a10fst1, lon, lat)
xysf1f <- st_as_sf(xy1f, coords = c("lon", "lat"), crs=4326)

xy5f <- dplyr::select(a10fst5, lon, lat)
xysf5f <- st_as_sf(xy5f, coords = c("lon", "lat"), crs=4326)

### Generate neighborhood matrix
nb1f <- dnearneigh(xysf1f, 0, 450)
nb5f <- dnearneigh(xysf5f, 0, 450)

### Moran test
FSTres1 <- residuals(mod_fstpw1)[,1]
FSTres5 <- residuals(mod_fstpw5)[,1]
moran.test(FSTres1, nb2listw(nb1f))
moran.test(FSTres5, nb2listw(nb5f))

## Ne
xy1n <- dplyr::select(a10ne1, lon, lat)
xysf1n <- st_as_sf(xy1n, coords = c("lon", "lat"), crs=4326)

xy5n <- dplyr::select(a10ne5, lon, lat)
xysf5n <- st_as_sf(xy5n, coords = c("lon", "lat"), crs=4326)

### Generate neighborhood matrix
nb1n <- dnearneigh(xysf1n, 0, 500)
nb5n <- dnearneigh(xysf5n, 0, 500)

### Moran test
Neres1 <- residuals(mod_nepw1)[,1]
Neres5 <- residuals(mod_nepw5)[,1]
moran.test(Neres1, nb2listw(nb1n))
moran.test(Neres5, nb2listw(nb5n))
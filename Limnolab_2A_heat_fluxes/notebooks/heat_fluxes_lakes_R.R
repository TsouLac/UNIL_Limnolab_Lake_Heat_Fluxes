#!/usr/bin/env Rscript
# ============================================================================
# Heat Flux Functions for Lake Analysis in R
# 
# Translated from Python to R by Damien Bouffard
# Original Python author: Damien Bouffard
# Email: damien.bouffard@eawag.ch
# 
# This file contains functions to calculate different components of heat
# fluxes for lake thermal analysis.
# ============================================================================

# Function for short-wave absorption
SW <- function(G, C, Adir) {
  # G: solar radiation (W/m2)
  # C: cloud coverage (fraction)
  # Adir: direct albedo
  
  Adiff <- 0.066  # diffuse albedo
  
  Fdir <- (1 - C) / ((1 - C) + 0.5 * C)
  Fdiff <- (0.5 * C) / ((1 - C) + 0.5 * C)
  
  sw <- (G * Fdir * (1 - Adir) + G * Fdiff * (1 - Adiff))
  return(sw)
}

# Function for absorption of atmospheric longwave radiation
LW_in <- function(Ta, Ea) {
  # Ta: air temperature (°C)
  # Ea: emissivity
  
  # Parameters
  A_L <- 0.03
  sigma <- 5.67e-8  # Stefan-Boltzmann constant
  
  # Absorption of atmospheric longwave radiation
  lwi <- (1 - A_L) * Ea * sigma * (Ta + 273.25)^4
  
  return(lwi)
}

# Function to calculate emissivity
emissivity_ <- function(Ta, ea, C) {
  # Ta: air temperature (°C)
  # ea: water vapor pressure (mbar)
  # C: cloud coverage
  
  # Parameter
  a <- 1.09
  
  # Emissivity
  Ea <- a * (1 + 0.17 * C^2) * 1.24 * (ea / (Ta + 273.25))^(1/7)
  
  return(Ea)
}

# Function for longwave emission from water surfaces
LW_out <- function(Tw) {
  # Tw: water temperature (°C)
  
  # Parameters
  E_w <- 0.972    # emissivity of water
  sigma <- 5.67e-8  # Stefan-Boltzmann constant
  
  # Longwave emission from water surfaces (negative, heat loss)
  lwo <- -E_w * sigma * (Tw + 273.25)^4
  
  return(lwo)
}

# Function for evaporation and condensation (latent heat flux)
He <- function(Ta, Tw, WS, ea) {
  # Ta: air temperature (°C)
  # Tw: water temperature (°C)
  # WS: wind speed (m/s)
  # ea: water vapor pressure (mbar)
  
  # Saturation vapor pressure over water
  es <- 6.112 * exp(17.62 * Tw / (243.12 + Tw))
  
  # Heat transfer coefficient
  f <- 4.8 + 1.98 * WS + 0.28 * (Tw - Ta)
  
  # Evaporation and condensation (negative, typically heat loss)
  he <- -f * (es - ea)
  
  return(he)
}

# Function for sensible heat flux
Hc <- function(Ta, Tw, WS) {
  # Ta: air temperature (°C)
  # Tw: water temperature (°C)
  # WS: wind speed (m/s)
  
  # Parameters
  cp <- 1005        # specific heat of air (J/(kg·K))
  p <- 1015         # atmospheric pressure (hPa)
  Lv <- 2.47e6      # latent heat of vaporization (J/kg)
  gamma <- cp * p / (0.622 * Lv)
  
  # Heat transfer coefficient
  f <- 4.8 + 1.98 * WS + 0.28 * (Tw - Ta)
  
  # Sensible heat flux (negative when Tw < Ta)
  hc <- -gamma * f * (Tw - Ta)
  
  return(hc)
}

# Function for throughflow heat flux
Hf <- function(Q, Tr, Tl, A0) {
  # Q: discharge (m3/s)
  # Tr: inflow temperature (°C)
  # Tl: lake temperature (°C)
  # A0: lake surface area (m2)
  
  # Parameters
  cp <- 4180      # specific heat of water (J/(kg·K))
  rho <- 1000     # density of water (kg/m3)
  
  # Throughflow heat flux
  hf <- -cp * rho * Q * (Tr - Tl) / A0
  
  return(hf)
}

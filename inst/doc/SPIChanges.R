## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(SPIChanges)

## ----eval=FALSE---------------------------------------------------------------
#  TSaggreg(daily.rain,
#    start.date,
#    TS = 4
#  )

## ----example 1----------------------------------------------------------------
library(SPIChanges)
daily.rain <- CampinasRain[, 2]
RainTS4 <- TSaggreg(daily.rain = daily.rain, start.date = "1980-01-01", TS = 4)
head(RainTS4)

## ----eval=FALSE---------------------------------------------------------------
#  SPIChanges(
#    rain.at.TS,
#    only.linear
#  )

## ----example 2----------------------------------------------------------------
library(SPIChanges)
daily.rain <- CampinasRain[, 2]
rainTS4 <- TSaggreg(daily.rain = daily.rain, start.date = "1980-01-01", TS = 4)
Changes_SPI <- SPIChanges(rain.at.TS = rainTS4, only.linear = "Yes")
head(Changes_SPI$data.week)
head(Changes_SPI$Statistics)
head(Changes_SPI$model.selection)
head(Changes_SPI$Changes.Freq.Drought)

## ----example 3----------------------------------------------------------------
library(SPIChanges)
daily.rain <- CampinasRain[, 2]
rainTS4 <- TSaggreg(daily.rain = daily.rain, start.date = "1980-01-01", TS = 4)
Changes_SPI <- SPIChanges(rain.at.TS = rainTS4, only.linear = "No")
head(Changes_SPI$data.week)
head(Changes_SPI$Statistics)
head(Changes_SPI$model.selection)
head(Changes_SPI$Changes.Freq.Drought)


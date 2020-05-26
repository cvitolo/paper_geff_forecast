# SET UP #######################################################################
# Install caliver from GitHub
# install.packages("devtools", repos = "http://cran.rstudio.com")
# devtools::install_github("ecmwf/caliver")

# Load relevant packages
library("caliver")
library("raster")
library("ggplot2")
library("rowr")
library("dplyr")
library("rgeos")
library("maptools")
library("sp")
library("sf")
library("RColorBrewer")
library("stringr")
library("colorspace")
library("easyVerification")
library("spatialEco")
library("doSNOW")
library("parallel")
library("gplots")

forecasts_folder <- "/hugetmp/forecasts/hres/2017"

dates2017 <- seq.Date(from = as.Date("2017-01-01"),
                      to = as.Date("2017-12-31"),
                      by = "day")

datesERA5 <- seq.Date(from = as.Date("1980-01-01"),
                      to = as.Date("2019-06-30"),
                      by = "day")

# Get all GFED4 regions
BasisRegions <- readRDS("data/BasisRegions_simplified.rds")

# Add the id column to the GFED4 regions for join
regionXFort <- fortify(BasisRegions, data = BasisRegions@data)
regionXPoly <- merge(regionXFort, BasisRegions@data,
                     by.x = "id", by.y = "ID")  # join data
# Re-order factors
regionXPoly$Region <- factor(as.character(regionXPoly$Region),
                             levels = BasisRegions$Region)

# GFED4 regions of interest
EURO <- BasisRegions[BasisRegions$Region == "EURO",]
TENA <- BasisRegions[BasisRegions$Region == "TENA",]
SHSA <- BasisRegions[BasisRegions$Region == "SHSA",]

# Countries of interest
pt <- raster::getData(name = "GADM", country = "Portugal", level = 1)[-c(2,13),]
chile <- raster::getData(name = "GADM", country = "Chile", level = 0)
cali <- raster::getData(name = "GADM", country = "USA", level = 1)[5, ]

# DATA PREPARATION #############################################################

# The following code uses forecast data which is restricted.
# A summary table is computed and saved for users at the end of this code chunk.
#
# Get data from previous publication
# https://github.com/cvitolo/GEFF-ERA5/blob/master/data/df_geff_erai_era5.rds
# save it in the "../data" folder
df <- readRDS(file = "df_geff_erai_era5.rds")
# insert 1-year fire season for tropics (lat = [-30, +30])
df$season[df$lat < 30 & df$lat > -30] <- "Dry"
# Remove wet season (we are only interested in the dry season)
df <- df[df$season == "Dry", ]
# How many stations remain?
stations <- unique(df[, c("id", "lat", "long")])

df$region <- sapply(strsplit(df$tzid, "/"), `[`, 1)
df$date <- as.Date(paste0(df$yr, "-", df$mon, "-", df$day))
df <- df[, c("id", "lat", "long", "region", "date", "OBS", "ERA5")]

days_2017 <- which(datesERA5 %in% c(dates2017, tail(dates2017)[6] + 1:9))
# system(paste0("cdo seltimestep,", paste0(days_2017, collapse = ","),
#                   " /scratch/rd/nen/perClaudia/era5/fwi_1980_2019.nc ",
#                   "/scratch/rd/nen/perClaudia/era5/fwi_2017_dayx.nc"),
#            ignore.stderr = TRUE)
ERA5_2017 <- raster::brick("/scratch/rd/nen/perClaudia/era5/fwi_2017_dayx.nc")

for (i in seq_along(dates2017)) {
  
  issue_date <- dates2017[i]
  
  # What stations have data on this date?
  dfx <- dplyr::filter(.data = df, date == issue_date)
  dfx <- dfx[complete.cases(dfx), ]
  
  # Define spatial point
  spdf <- sf::st_as_sf(x = dfx, coords = c("long", "lat"), crs = 4326)
  
  # Get HRES
  message(paste("Handling HRES for", issue_date))
  HRES <- raster::brick(file.path(forecasts_folder,
                                  paste0("ECMWF_FWI_",
                                         gsub("-", "",
                                              as.character(issue_date)),
                                         "_1200_hr_fwi.nc")))
  # HRES needs to be rotated
  df_hres <- raster::extract(x = raster::rotate(HRES), y = spdf)
  dfx <- cbind(dfx, df_hres)
  names(dfx)[8:17] <- paste0("HRES_d", 1:10)
  
  # Get ENS
  arr_ens <- array(NA, dim = c(dim(spdf)[1], 10, 51))
  message(paste("Handling ENS for", issue_date))
  for (j in 1:51){
    # Extract the modelled FWI from ENS
    ens_member <- sprintf("%02d", j - 1)
    ENS <- raster::brick(file.path(forecasts_folder,
                                   paste0("ECMWF_FWI_",
                                          gsub("-", "",
                                               as.character(issue_date)),
                                          "_1200_", ens_member,
                                          "_fwi.nc")))[[1:10]]
    arr_ens[, , j] <- raster::extract(x = raster::rotate(ENS), y = spdf)
  }
  
  # Get ERA5 data for 2017 - to be used as observation
  message(paste("Handling ERA5 for", issue_date, "plus 9 days"))
  ERA5_2017_day <- ERA5_2017[[which(dates2017 %in% issue_date) + 0:9]]
  day_in_2017 <- raster::extract(x = ERA5_2017_day, y = spdf)
  
  for (leadtime in 1:10){
    
    lead_day <- issue_date + leadtime - 1
    
    # Extract ERA5 from all the years (same month and day), excluding 2017
    idx_in_datesERA5 <- which(format(datesERA5, "%m-%d") %in%
                                format(lead_day, "%m-%d") &
                                lubridate::year(datesERA5) != "2017")
    system(paste0("cdo seltimestep,", paste0(idx_in_datesERA5, collapse = ","),
                  " /scratch/rd/nen/perClaudia/era5/fwi_1980_2019.nc ",
                  "/scratch/rd/nen/perClaudia/era5/fwi_1980_2019_lead_day.nc"),
           ignore.stderr = TRUE)
    ERA5 <- raster::brick("/scratch/rd/nen/perClaudia/era5/fwi_1980_2019_lead_day.nc")
    
    # ERA5 was stored with longitudes already rotated to [-180, +180],
    # HRES and ENS will need to be rotated
    df_era5_39 <- raster::extract(x = ERA5, y = spdf)
    # Sample ERA5 to get 51 random ensemble members
    df_era5_51 <- df_era5_39[, sample(x = 1:dim(df_era5_39)[2],
                                      size = 51, replace = TRUE)]
    
    # Populate CRPS
    crps_fc <- easyVerification::veriApply(verifun = 'EnsCrps',
                                           fcst = arr_ens[, leadtime,],
                                           obs = day_in_2017[, leadtime])
    crps_clim <- easyVerification::veriApply(verifun = 'EnsCrps',
                                             fcst = df_era5_51,
                                             obs = day_in_2017[, leadtime])
    crpss <- 1 - crps_fc/crps_clim
    
    dfx <- cbind(dfx, crps_fc, crps_clim, crpss)
  }
  names(dfx)[18:47] <- paste0(c("CRPS_fc_d", "CRPS_clim_d", "CRPSS_d"),
                              rep(1:10, each = 3))
  
  if (i == 1) {
    result <- dfx
  }else{
    result <- rbind(result, dfx)
  }
  
  saveRDS(result, "df_geff_era5_hres_crps.rds")
  
}

result <- readRDS("df_geff_era5_hres_crps.rds")

# Remove rows with NAs and update spatial points
rows2remove <- which(rowSums(is.na(result)) > 0)
if (length(rows2remove) > 0) result <- result[-rows2remove, ]

# Define spatial point
spdf <- sf::st_as_sf(x = result, coords = c("long", "lat"), crs = 4326)

synops <- spatialEco::point.in.poly(spdf, BasisRegions)

# FIGURE 1 #####################################################################

ggplot(data = map_data("world"), aes(x = long, y = lat)) +
  geom_polygon(aes(group = group), fill=NA, colour = "grey65") +
  coord_equal() +  theme_bw() + xlab("Longitude") + ylab("Latitude") +
  geom_polygon(data = regionXPoly,
               aes(x = long, y = lat, group = group, fill = factor(Region)),
               alpha = 0.5) +
  scale_fill_manual(name="GFED4 regions",
                    values = c("darkgrey",
                               RColorBrewer::brewer.pal(12,"Paired"),
                               "lightgrey"),
                    labels = BasisRegions$Region) +
  geom_polygon(data = map_data("world", region = "Portugal"),
               aes(x = long, y = lat, group = group, colour = "Portugal"),
               fill = NA) +
  geom_polygon(data = map_data("world", region = "Chile"),
               aes(x = long, y = lat, group = group, colour = "Chile"),
               fill = NA) +
  geom_polygon(data = map_data(map = "state", region = "California"),
               aes(x = long, y = lat, group = group, colour = "California"),
               fill = NA) +
  geom_point(data = data.frame(synops), aes(x = coords.x1, y = coords.x2),
             shape = 1, color = "gray20", size = 0.1) +
  scale_colour_manual(name = "Study areas",
                      values = c("Portugal" = "brown",
                                 "California" = "darkblue",
                                 "Chile" = "darkgreen"))

# FIGURE 2 #####################################################################

f2 <- data.frame(synops) %>%
  group_by(Region) %>%
  select(OBS, Region, ERA5, HRES_d1, HRES_d2, HRES_d3, HRES_d4, HRES_d5,
         HRES_d6, HRES_d7, HRES_d8, HRES_d9, HRES_d10) %>%
  reshape2::melt(id.vars = c("OBS", "Region")) %>%
  filter(complete.cases(.)) %>% # Remove NAs
  filter(OBS < 250, value < 250) %>%
  mutate(BIAS = OBS - value,
         MAE = abs(OBS - value)) %>%
  rename(modelled_type = variable, modelled_value = value) %>%
  reshape2::melt(id.vars = c("OBS", "Region", "modelled_type",
                             "modelled_value"))

ggplot(f2, aes(x = modelled_type, y = value)) +
  facet_grid(Region ~ variable) +
  geom_boxplot(color="black", outlier.shape = NA) +
  theme_bw() + ylim(-40, 50) +
  xlab("") + ylab("") +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("ERA5", "D1", "D2", "D3",
                              "D4", "D5", "D6", "D7",
                              "D8", "D9", "D10"))

# FIGURE 3 #####################################################################

f3 <- data.frame(synops) %>%
  select(Region, CRPS_fc_d1, CRPS_fc_d2, CRPS_fc_d3, CRPS_fc_d4, CRPS_fc_d5, CRPS_fc_d6,
         CRPS_fc_d7, CRPS_fc_d8, CRPS_fc_d9, CRPS_fc_d10) %>%
  reshape2::melt(id.vars = "Region") %>%
  filter(complete.cases(.)) %>%
  rename(CRPS = value) %>%
  group_by(Region, variable) %>%
  summarise(CRPS = mean(CRPS))

f3$variable <- as.numeric(f3$variable)
f3$Region <- factor(f3$Region)

ggplot(f3, aes(x = variable, y = CRPS, colour = Region)) +
  geom_line() +
  geom_point() +
  xlab("Lead time") +
  theme_bw() +
  scale_x_discrete(limits = 1:10) +
  scale_colour_manual(name="Regions",
                      values = c("darkgrey",
                                 RColorBrewer::brewer.pal(12,"Paired"),
                                 "lightgrey"),
                      labels = BasisRegions$Region)

# FIGURE 4 #####################################################################

f4 <- data.frame(synops) %>%
  filter(complete.cases(.)) %>%
  select(paste0(c("CRPS_fc_d", "CRPS_clim_d"), rep(1:10, each = 2))) %>%
  reshape2::melt(id.vars = c("CRPS_clim_d1", "CRPS_clim_d2", "CRPS_clim_d3",
                             "CRPS_clim_d4", "CRPS_clim_d5", "CRPS_clim_d6",
                             "CRPS_clim_d7", "CRPS_clim_d8", "CRPS_clim_d9",
                             "CRPS_clim_d10")) %>%
  rename(CRPS_fc = value) %>%
  rename(fc_name = variable) %>%
  reshape2::melt(id.vars = c("fc_name", "CRPS_fc")) %>%
  rename(CRPS_clim = value) %>%
  rename(clim_name = variable) %>%
  group_by(fc_name, clim_name) %>%
  summarise(CRPS_fc = mean(CRPS_fc),
            CRPS_clim = mean(CRPS_clim))

f4$fc_name <- as.numeric(f4$fc_name)
f4$clim_name <- as.numeric(f4$clim_name)

f4plot <- ggplot(f4) +
  geom_line(aes(x = fc_name, y = CRPS_fc, col = "Forecast"), size = 2) +
  geom_line(aes(x = clim_name, y = CRPS_clim, col = "Climatology"), size = 2) +
  xlab("Lead time") + ylab("CRPS") +
  scale_color_discrete(name = "") +
  theme_bw() +
  scale_x_discrete(limits = 1:10) +
  theme(text = element_text(size = 15))

ggsave(filename = "../images/Figure4.png", plot = f4plot)

# FIGURE 5 #####################################################################

f5 <- data.frame(synops) %>%
  filter(complete.cases(.)) %>%
  select(Region, paste0(c("CRPS_fc_d", "CRPS_clim_d"), rep(1:10, each = 2))) %>%
  reshape2::melt(id.vars = c("Region",
                             "CRPS_clim_d1", "CRPS_clim_d2", "CRPS_clim_d3",
                             "CRPS_clim_d4", "CRPS_clim_d5", "CRPS_clim_d6",
                             "CRPS_clim_d7", "CRPS_clim_d8", "CRPS_clim_d9",
                             "CRPS_clim_d10")) %>%
  rename(CRPS_fc = value) %>%
  rename(fc_name = variable) %>%
  reshape2::melt(id.vars = c("Region", "fc_name", "CRPS_fc")) %>%
  rename(CRPS_clim = value) %>%
  rename(clim_name = variable) %>%
  group_by(Region, fc_name, clim_name) %>%
  summarise(CRPS_fc = mean(CRPS_fc),
            CRPS_clim = mean(CRPS_clim))

f5$fc_name <- as.numeric(f5$fc_name)
f5$clim_name <- as.numeric(f5$clim_name)
f5$Region <- factor(f5$Region)

ggplot(f5) +
  geom_line(aes(x = fc_name, y = CRPS_fc, col = "Forecast")) +
  geom_line(aes(x = clim_name, y = CRPS_clim, col = "Climatology")) +
  facet_wrap( ~ Region, ncol = 2, scales = "free_y") +
  xlab("Lead time") + ylab("CRPS") +
  scale_color_discrete(name = "") +
  theme_bw() +
  scale_x_discrete(limits = 1:10)

# PREPARATION FIGURE 6-7 #######################################################

# Get ERA5
ERA5 <- raster::brick("/scratch/rd/nen/perClaudia/era5/fwi_1980_2019.nc")
# Get FRP
FRP <- raster::brick("/scratch/rd/nen/perClaudia/CAMS/CAMS_2017-01-01_2017-12-31_frpfire.nc")
myfilelist <- list.files(path = forecasts_folder,
                         pattern = "*hr",
                         full.names = TRUE)

ncores <- parallel::detectCores() - 1
cl <- makeCluster(ncores)
registerDoSNOW(cl)
iterations <- 365
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
result <- foreach(i = 10:iterations, .combine = rbind, 
                  .options.snow = opts) %dopar%
  {
    # Get FRP for 2017, keep only major fires (values above 0.5)
    frpx <- raster::rotate(FRP[[i]])
    p <- raster::rasterToPoints(frpx, fun = function(x){x > 0.5})
    spdf <- data.frame(p)
    spdf$date <- dates2017[i]
    names(spdf) <- c("long", "lat", "frp", "date")
    sp::coordinates(spdf) <- ~long + lat
    n_points <- dim(spdf@data)[1]
    
    # Get HRES for 2017, where FRP is above 0.5
    lead_day <- 0
    hres_df <- data.frame(matrix(NA, ncol = 10, nrow = n_points))
    for (j in i:(i - 9)){
      lead_day <- lead_day + 1
      hres <- raster::rotate(raster::brick(myfilelist[[j]])[[lead_day]])
      x <- raster::extract(x = hres, y = spdf)
      hres_df[, lead_day] <- x
    }
    names(hres_df) <- paste0("hres_day", 1:10)
    
    # Get ERA5 to calculate various warning levels
    # Low = 0.75, Moderate = 0.85, High = 0.90, Very high = 0.95, Extreme = 0.98
    median <- caliver::daily_clima(r = ERA5, dates = dates2017[i], probs = 0.50)
    moderate <- caliver::daily_clima(r = ERA5, dates = dates2017[i], probs = 0.85)
    high <- caliver::daily_clima(r = ERA5, dates = dates2017[i], probs = 0.90)
    very_high <- caliver::daily_clima(r = ERA5, dates = dates2017[i], probs = 0.95)
    CLIM <- raster::extract(x = median, y = spdf)
    thr85 <- raster::extract(x = high, y = spdf)
    thr90 <- raster::extract(x = high, y = spdf)
    thr95 <- raster::extract(x = very_high, y = spdf)
    
    # Append results to table
    df <- cbind(spdf@coords, spdf@data, CLIM, thr85, thr90, thr95, hres_df)
    df <- df[complete.cases(df), ]
    return(df)
  }
close(pb)
stopCluster(cl)
# saveRDS(result, "PODdataframe.rds")

# Load, if pre-calculate
# result <- readRDS("PODdataframe.rds")
# BasisRegions <- readRDS("BasisRegions_simplified.rds")

# FRP locations of interest
sp::coordinates(result) <- ~long + lat
crs(result) <- crs(BasisRegions)

# Add ID, Region and Zone to the dataset vias spatial intersection
fires <- spatialEco::point.in.poly(result, BasisRegions)
fires <- as.data.frame(fires)
fires <- fires[complete.cases(fires), ]

# Check coverage
# plot(BasisRegions)
# plot(result, col = "red", add = TRUE)

# If all looks good, we can remove previous datasets
# rm(result, BasisRegions)

# saveRDS(fires, "fires.rds")
# fires <- readRDS("fires.rds")

# Check if CLIM and HRES are above threshold (90th percentile)
my_thr <- fires$thr90
fires$CLIM <- fires$CLIM >= my_thr
for (colx in which(names(fires) %in% paste0("hres_day", 1:10))){
  x <- fires[, colx] > my_thr
  fires <- cbind(fires, x)
}
names(fires)[which(names(fires) == "x")] <- paste0("HDay", 1:10)
rm(x, colx)

# Group by Region and calculate POD
fires_binary <- fires %>%
  mutate(FRP_binary = TRUE) %>%
  select(Region, coords.x1, coords.x2, FRP_binary, CLIM,
         HDay1, HDay2, HDay3, HDay4, HDay5,
         HDay6, HDay7, HDay8, HDay9, HDay10) %>%
  reshape2::melt(id.vars = c("Region", "coords.x1", "coords.x2", "FRP_binary")) %>%
  group_by(Region, coords.x1, coords.x2, variable) %>%
  summarise(POD = verify(obs = FRP_binary, pred = value,
                         frcst.type = "binary", obs.type = "binary")$POD) %>%
  group_by(Region, variable) %>%
  summarise(PODmean = mean(POD))
# saveRDS(fires_binary, "fires_binary_90.rds")

# FIGURE 6 #####################################################################

fires_binary <- readRDS("fires_binary_90.rds")
fires_binary <- fires_binary[-which(fires_binary$variable == "CLIM"), ]

ggplot(fires_binary, aes(variable, Region, fill = PODmean)) + 
  geom_tile() +
  geom_text(aes(label = round(PODmean, 2))) +
  scale_fill_gradient(name = "POD", low = "white", high = "brown") +
  theme_bw() + xlab("") + ylab("")

# FIGURE 7 #####################################################################

ggplot(data = map_data("world"), aes(x = long, y = lat)) +
  geom_polygon(aes(group = group), fill=NA, colour = "grey65") +
  coord_equal() +  theme_bw() + xlab("Longitude") + ylab("Latitude") +
  geom_polygon(data = regionXPoly,
               aes(x = long, y = lat, group = group, fill = factor(Region)),
               alpha = 0.5) +
  scale_fill_manual(name="GFED4 regions",
                    values = c("darkgrey",
                               RColorBrewer::brewer.pal(12,"Paired"),
                               "lightgrey"),
                    labels = BasisRegions$Region) +
  geom_point(data = data.frame(result), aes(x = long, y = lat),
             shape = 1, color = "gray20", size = 0.1)

# FIGURE 8 #####################################################################

# Chile
dates_chile <- seq.Date(from = as.Date("2017-01-01"),
                        to = as.Date("2017-01-31"),
                        by = "day")
bbox_chile <- as(raster::extent(-74, -70, -37, -33), "SpatialPolygons")
raster::plot(chile, col = "lightgray")
raster::plot(bbox_chile, lwd = 2, border = "red", add = TRUE)
# Crop manually to mainland!

# Portugal
dates_pt <- seq.Date(from = as.Date("2017-06-01"),
                     to = as.Date("2017-06-30"),
                     by = "day")
bbox_pt <- as(raster::extent(-8.5, -7.6, 39.6, 40.3), "SpatialPolygons")
pt <- gUnaryUnion(pt, id = pt@data$NAME_0)
raster::plot(pt, col = "lightgray")
raster::plot(bbox_pt, lwd = 2, border = "red", add = TRUE)
# Crop manually to mainland!

# California
dates_cali <- seq.Date(from = as.Date("2017-09-21"),
                       to = as.Date("2017-10-20"),
                       by = "day")
bbox_cali <- as(raster::extent(-123.4, -117.6, 33.7, 42), "SpatialPolygons")
raster::plot(cali, col = "lightgray")
raster::plot(bbox_cali, lwd = 2, border = "red", add = TRUE)
# Crop manually to mainland!

obs_file <- "/scratch/rd/nen/perClaudia/CAMS/CAMS_2017-01-01_2017-12-31_frpfire.nc"
clima_file <- "/scratch/rd/nen/perClaudia/era5/fwi_era5_1980_2016_90th_daily_clima.nc"

# Chile
df_chile <- caliver:::make_forecast_summary(input_dir = "/hugetmp/forecasts/fwi",
                                            p = bbox_chile,
                                            event_dates = dates_chile,
                                            obs = obs_file,
                                            clima = clima_file)
saveRDS(df_chile, "df_chile.rds")
# Portugal
df_pt <- caliver:::make_forecast_summary(input_dir = "/hugetmp/forecasts/fwi",
                                         p = bbox_pt,
                                         event_dates = dates_pt,
                                         obs = obs_file,
                                         clima = clima_file)
saveRDS(df_pt, "df_pt.rds")
# California
df_cali <- caliver:::make_forecast_summary(input_dir = "/hugetmp/forecasts/fwi",
                                           p = bbox_cali,
                                           event_dates = dates_cali,
                                           obs = obs_file,
                                           clima = clima_file)
saveRDS(df_cali, "df_cali.rds")

# Chile
plot_chile <- caliver:::plot_forecast_summary(df_chile)
# Portugal
# df_pt <- readRDS("df_pt.rds")
plot_pt <- caliver:::plot_forecast_summary(df_pt)
# California
# df_cali <- readRDS("df_cali.rds")
plot_cali <- caliver:::plot_forecast_summary(df_cali)

# APPENDIX A ###################################################################

library("ncdf4")
library("raster")
library("sf")
library("spatialEco")

setwd("/perm/mo/moc0/repos/paper_fire_forecast")

dummy <- rotate(raster(paste0("/hugetmp/reanalysis/GEFF-ERA5/hres/fwi/ECMWF_FWI_20170101_1200_hr_fwi.nc")))
# dummy <- rotate(raster(paste0("~/ECMWF_FWI_20170101_1200_hr_fwi.nc")))

BasisRegions <- readRDS("data/BasisRegions_simplified.rds")
BasisRegions <- rasterize(x = BasisRegions, y = dummy)
BasisRegions <- raster::shift(raster::rotate(raster::shift(BasisRegions, 180)), 180)

dummy <- raster(paste0("/hugetmp/reanalysis/GEFF-ERA5/hres/fwi/ECMWF_FWI_20170101_1200_hr_fwi.nc"))
# dummy <- raster(paste0("~/ECMWF_FWI_20170101_1200_hr_fwi.nc"))
pts <- xyFromCell(object = dummy, cell = which(!is.na(dummy[])), )
pts <- SpatialPoints(coords = pts, proj4string = crs(dummy))

# Remove NAs
pts$region <- extract(BasisRegions, pts)
ptsx <- sp.na.omit(pts, margin = 1)

# check overlap
plot(dummy)
plot(BasisRegions, add = TRUE)
plot(ptsx, add = TRUE)

rm(dummy, pts)

issue_dates <- seq.Date(from = as.Date("2017-01-01"),
                        to = as.Date("2017-12-31"),
                        by = "day")

arr <- array(data = NA, dim = c(length(issue_dates), 10, 14))

for (i in seq_along(issue_dates)){
  issue_date <- issue_dates[i]
  print(issue_date)
  issuedate <- gsub(pattern = "-", replacement = "", x = issue_date)
  fc_file <- paste0("/hugetmp/forecasts/hres/2017/ECMWF_FWI_",
                    issuedate, "_1200_hr_fwi.nc")
  fc <- raster::brick(fc_file)
  for (leadtime in 1:10){
    lead_date <- gsub(pattern = "-", replacement = "",
                      x = issue_date + leadtime - 1)
    re_file <- paste0("/hugetmp/reanalysis/GEFF-ERA5/hres/fwi/ECMWF_FWI_",
                      lead_date, "_1200_hr_fwi.nc")
    re <- raster::raster(re_file)
    df_re <- extract(x = re, y = ptsx, df = TRUE)
    df_fc <- extract(x = fc[[leadtime]], y = ptsx, df = TRUE)
    bias <- df_fc[[2]] - df_re[[2]]
    for (region in 1:14){
      arr[i, leadtime, region] <- mean(bias[which(ptsx$region == region)],
                                       na.rm = TRUE)
    }
  }
}
saveRDS(arr, "data/bias_array.rds")

arr <- readRDS("data/bias_array.rds")
df <- data.frame(matrix(NA, ncol = 14, nrow = 10))
for (leadtime in 1:10){
  for (region in 1:14){
    df[leadtime, region] <- mean(arr[, leadtime, region])
  }
}
saveRDS(df, "data/bias_df.rds")
# df <- readRDS("data/bias_df.rds")
names(df) <- BasisRegions$Region
df$leadtime <- 1:10

df_melt <- tidyr::pivot_longer(df, cols = 1:14)
names(df_melt)[2] <- "Region"

# library
library(ggplot2)
library(viridis)
library(hrbrthemes)

palette <- rev(c("#E9E9E9", "#D8AB93", "#FFFECC", "#B49DCC", "#E4D8EA", "#FCBE7F",
             "#E4D8EA", "#F18C8D", "#FDCCCC", "#98CF95", "#D8EFC4", "#8EBBD9",
             "#D2E6F1", "#D4D4D4"))

# Grouped
ggplot(df_melt, aes(fill = Region, y = value, x = as.factor(leadtime))) + 
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = palette, labels = BasisRegions$Region) +
  # scale_fill_viridis(discrete = T) +
  #ggtitle("Mean bias by region and leadtime") +
  theme_bw() +
  xlab("Leadtime") + ylab("Bias")

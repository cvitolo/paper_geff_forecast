library("dplyr")
library("ggplot2")
library("colorspace")
library("easyVerification")
library("manipulate")

# Read FRP with longitudes in the range [-180, +180]
FRP <- rotate(brick("/scratch/rd/nen/perClaudia/CAMS/CAMS_2017-01-01_2017-12-31_frpfire.nc"))
for (i in seq_along(dates2017)){

  print(dates2017[i])


  # ERA5 was stored with longitudes already rotated to [-180, +180],
  # HRES and ENS will need to be rotated
  df_era5 <- raster::extract(x = ERA5, y = spdf)

  # Get HRES
  message(paste("Handling HRES for", dates2017[i], "day", leadtime))
  HRES <- raster::brick(file.path(forecasts_folder,
                                  paste0("ECMWF_FWI_",
                                         gsub("-", "",
                                              as.character(dates2017[i])),
                                         "_1200_hr_fwi.nc")))[[leadtime]]
  df_hres <- raster::extract(x = rotate(HRES), y = spdf)

}

df_used <- readRDS("/perm/mo/moc0/repos/paper_fire_forecast/data/df_geff_era5_hres_crps.rds")

plot(as.numeric(df_used@data[1, seq(8, 35, 3)]),
     main = "Growth of CRPS with leadtime", type = "l", xlab = "", ylab = "CRPS")

plot(ecdf(df_era5[1,]), main = "ENS vs ERA5clim", xlab = "")
plot(ecdf(df_ens[1,]), col = "red", add = TRUE)
legend("bottom",
       c("ERA5clim", "ENS"),
       pch = c(1, 1),
       lty = c(1, 1),
       lwd = c(1, 1),
       col = c("black", "red"),
       inset = c(0,1), xpd = TRUE, horiz = TRUE, bty = "n")

d <- df_used@data %>%
  group_by(region) %>%
  select(OBS, region, ERA5, HRES_d2, HRES_d6, HRES_d10) %>%
  reshape2::melt(id.vars = c("OBS", "region")) %>%
  filter(complete.cases(.)) %>% # Remove NAs
  filter(OBS < 250, value < 250) %>%
  mutate(BIAS = OBS - value,
         MAE = abs(OBS - value)) %>%
  rename(modelled_type = variable, modelled_value = value) %>%
  reshape2::melt(id.vars = c("OBS", "region", "modelled_type", "modelled_value"))

ggplot(d, aes(x = modelled_type, y = value, fill = modelled_type)) +
  facet_grid(region ~ variable) +
  geom_boxplot(color="black", size=0.2,  outlier.shape = NA) +
  theme_bw() + ylim(-40, 50) +
  xlab("") + ylab("") +
  scale_fill_discrete_qualitative() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("ERA5", "HRES_D2", "HRES_D6", "HRES_D10"))

d <- df_used@data %>%
  select(region, CRPS_d2, CRPS_d6, CRPS_d10) %>%
  rename(Day2 = CRPS_d2, Day6 = CRPS_d6, Day10 = CRPS_d10) %>%
  reshape2::melt(id.vars = "region") %>%
  filter(complete.cases(.)) %>%
  rename(CRPS = value, Region = region) %>%
  group_by(Region, variable) %>%
  summarise(CRPS = mean(CRPS))

ggplot(d, aes(x = variable, y = CRPS, colour = Region)) +
  geom_point() +
  xlab("") + #ylim(3, 7.5) +
  scale_fill_manual(name="Region")

easyVerification::veriApply(verifun = "EnsCrps",
                            fcst = mat2[idx_to_keep[1], ],
                            obs = dfx$OBS[idx_to_keep[1]])

#for (i in 1:51) {plot(ecdf(mat2[, i]), col = "red", add = TRUE)}

# which point?

manipulate({
minx <- min(mat2[mypoint, ], mat6[mypoint, ], mat10[mypoint, ],
            dfx$ERA5[mypoint], dfx$OBS[mypoint])
maxx <- max(mat2[mypoint, ], mat6[mypoint, ], mat10[mypoint, ],
            dfx$ERA5[mypoint], dfx$OBS[mypoint])
plot(ecdf(mat2[mypoint, ]), col = "yellow",
     xlim = c(floor(minx/5)*5, ceiling(maxx/5)*5), xlab = "",
     main = paste0("Fire Weather Index in ", dfx$region[mypoint],
                   " (OBS = ", round(dfx$OBS[mypoint]),
                   ", ERA5 = ", round(dfx$ERA5[mypoint]), ")"))
plot(ecdf(mat6[mypoint, ]), col = "orange", add = TRUE)
plot(ecdf(mat10[mypoint, ]), col = "brown", add = TRUE)
plot(ecdf(dfx$HRES_d2[mypoint]), add = TRUE, verticals=TRUE, do.points=FALSE, col = "yellow")
plot(ecdf(dfx$HRES_d6[mypoint]), add = TRUE, verticals=TRUE, do.points=FALSE, col = "orange")
plot(ecdf(dfx$HRES_d10[mypoint]), add = TRUE, verticals=TRUE, do.points=FALSE, col = "brown")
plot(ecdf(dfx$OBS[mypoint]), add = TRUE, verticals=TRUE, do.points=FALSE, lty = 2, lwd = 3)
plot(ecdf(dfx$ERA5[mypoint]), add = TRUE, verticals=TRUE, do.points=FALSE, col = "blue", lty = 3, lwd = 3)
legend("bottom",
       c("OBS", "ERA5", "HRES2", "HRES6", "HRES10",
         paste("ENS2", ), "ENS6", "ENS10"),
       pch = c(NA, NA, NA, NA, NA, 1, 1, 1),
       lty = c(2, 3, 1, 1, 1, NA, NA, NA),
       lwd = c(3, 3, 1, 1, 1, NA, NA, NA),
       col = c("black", "blue",
               "yellow", "orange", "brown",
               "yellow", "orange", "brown"),
       inset = c(0,1), xpd = TRUE, horiz = TRUE, bty = "n")
}, mypoint = slider(1, 160))


# Tests
x <- verification::crpsDecomposition(obs = dfx$OBS, eps = mat2)
nObs<-length(dfx$OBS)

x <- sort(mat2[1, ])
xa <- dfx$OBS[1]
crps_decomposed <- verification::crpsDecomposition(obs = dfx$OBS, eps = mat2)

for (k in 1:51) {probi <- k/51; ci <- crps_decomposed$alpha[,k]*probi^2+crps_decomposed$beta[,k]*(1-probi)^2}

#```{r Figure2, echo = FALSE, fig.width = 13, fig.height = 6, out.width = "100%", fig.cap = "\\label{fig:Figure2}Spatial distribution of weather stations from the synop network  which have at least 30 observations recorded at local noon in 2017."}
# Get data from previous publication
# https://github.com/cvitolo/GEFF-ERA5/blob/master/data/df_all.rds
# save them in the "../data" folder
# Use QGIS to beautify the map
knitr::include_graphics("../images/synops_col.png")
#```


### OLD CODE

# The following code uses forecast data which is restricted.
# A summary table is computed and saved for users at the end of this code chunk.
# Get data from previous publication
# https://github.com/cvitolo/GEFF-ERA5/blob/master/data/df_geff_erai_era5.rds
# save it in the "../data" folder
# df <- readRDS(file = "../data/df_geff_erai_era5.rds")
# # insert 1-year fire season for tropics (lat = [-30, +30])
# df$season[df$lat < 30 & df$lat > -30] <- "Dry"
# # Remove wet season (we are only interested in the dry season)
# df <- df[df$season == "Dry", ]
# # How many stations remain?
# stations <- unique(df[, c("id", "lat", "long")])
#
# df$region <- sapply(strsplit(df$tzid, "/"), `[`, 1)
# df$date <- as.Date(paste0(df$yr, "-", df$mon, "-", df$day))
# df <- df[, c("id", "lat", "long", "region", "date", "OBS", "ERA5")]
#
# dates2017 <- seq.Date(from = as.Date("2017-01-01"),
#                       to = as.Date("2017-12-31"),
#                       by = "day")
#
# df_used <- data.frame(matrix(NA, nrow = 0, ncol = 27))
# names(df_used) <- c("id", "lat", "long", "region", "date", "OBS", "ERA5",
#                     paste0(c("HRES_d", "CRPS_d"), rep(1:10, each =2)))
# for (i in seq_along(dates2017)){
#
#   print(dates2017[i])
#
#   # What stations have data on this date?
#   dfx <- df %>% filter(date == dates2017[i]) %>% filter(complete.cases(.))
#   # Define spatial point
#   spdf <- dfx[, c("id", "lat", "long")]
#   # Promote data.frame to spatial
#   coordinates(spdf) = ~long + lat
#   # Get data from ENS forecasts
#   nstations <- dim(spdf@data)[1]
#
#   # Get data from HRES forecasts, e.g. ECMWF_FWI_20170101_1200_hr_fwi.nc
#   if ((dates2017[i] - 9) %in% dates2017){
#
#     for (leadtime in 1:10){
#
#       lead_day <- dates2017[i] - leadtime + 1
#       # HRES
#       HRES <- raster::brick(file.path(forecasts_folder,
#                                       paste0("ECMWF_FWI_",
#                                              gsub("-", "",
#                                                   as.character(lead_day)),
#                                              "_1200_hr_fwi.nc")))[[leadtime]]
#
#       dfx <- cbind(dfx, raster::extract(x = rotate(HRES), y = spdf))
#       names(dfx)[ncol(dfx)] <- paste0("HRES_d", leadtime)
#
#       # ENS
#       ens_mat <- matrix(NA, nrow = nstations, ncol = 51)
#       for (j in 1:51){
#         # Extract the modelled FWI from ENS
#         ens_member <- sprintf("%02d", j - 1)
#
#         ENS <- raster::brick(file.path(forecasts_folder,
#                                        paste0("ECMWF_FWI_",
#                                               gsub("-", "",
#                                                    as.character(lead_day)),
#                                               "_1200_", ens_member,
#                                               "_fwi.nc")))[[leadtime]]
#         ens_mat[, j] <- raster::extract(x = rotate(ENS), y = spdf)
#       }
#       # Populate CRPS
#       dfx <- cbind(dfx, veriApply(verifun = "EnsCrps",
#                                   fcst = ens_mat, obs = dfx$ERA5))
#       names(dfx)[ncol(dfx)] <- paste0("CRPS_d", leadtime)
#
#     }
#     # Append to table
#     df_used <- dplyr::bind_rows(df_used, dfx)
#
#     saveRDS(df_used, "../data/df_geff_era5_hres_crps.rds")
#   }
# }

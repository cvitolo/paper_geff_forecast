library("dplyr")
library("ggplot2")
library("colorspace")
library("easyVerification")
library("manipulate")

df_used <- readRDS("/home/mo/moc0/paper_fire_forecast/data/df_geff_era5_hres_crps.rds")

d <- df_used %>%
  select(OBS, ERA5, HRES_d2, HRES_d6, HRES_d10) %>%
  filter(complete.cases(.)) %>% # Remove NAs
  reshape2::melt(id.vars = "OBS") %>%
  filter(OBS < 250, value < 250) %>%
  mutate(BIAS = OBS - value,
         MAE = abs(OBS - value)) %>%
  rename(modelled_type = variable, modelled_value = value) %>%
  reshape2::melt(id.vars = c("OBS", "modelled_type", "modelled_value"))

ggplot(d, aes(x = modelled_type, y = value, fill = modelled_type)) +
  facet_wrap( ~ variable, ncol = 2) +
  geom_boxplot(color="black", size=0.2,  outlier.shape = NA) +
  theme_bw() + ylim(-40, 50) +
  xlab("") + ylab("") +
  scale_fill_discrete_qualitative() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("ERA5", "HRES_D2", "HRES_D6", "HRES_D10"))

d <- df_used %>%
  select(region, CRPS_d2, CRPS_d6, CRPS_d10) %>%
  rename(Region = region, Day2 = CRPS_d2, Day6 = CRPS_d6, Day10 = CRPS_d10) %>%
  reshape2::melt(id.vars = "Region") %>%
  filter(complete.cases(.)) %>%
  rename(CRPS = value) %>%
  group_by(Region, variable) %>% summarise(CRPS = mean(CRPS))

ggplot(d, aes(x = variable, y = CRPS, colour = Region)) +
  geom_point() +
  xlab("") + ylim(3, 7.5) +
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
plot(ecdf(dfx$OBS[mypoint]), add = TRUE, verticals=TRUE, do.points=FALSE, lty = 2, lwd = 2)
plot(ecdf(dfx$ERA5[mypoint]), add = TRUE, verticals=TRUE, do.points=FALSE, col = "blue", lty = 3, lwd = 2)
legend("bottom",
       c("OBS", "ERA5", "HRES2", "HRES6", "HRES10", "ENS2", "ENS6", "ENS10"),
       pch = c(NA, NA, NA, NA, NA, 1, 1, 1), lty = c(1, 1, 1, 1, 1, NA, NA, NA),
       col = c("black", "blue",
               "yellow", "orange", "brown",
               "yellow", "orange", "brown"),
       inset = c(0,1), xpd = TRUE, horiz = TRUE, bty = "n")
}, mypoint = slider(1, 160))


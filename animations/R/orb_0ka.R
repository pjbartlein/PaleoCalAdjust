# orb_0ka.R
# orbital positions at 5-degree intervals at 0 ka

# functions

# degrees and radians
# input: angle in degrees or radians
# returns angle in radians or degrees
degree <- function(rad) {(rad * 180) / (pi)}
radians <- function(deg) {(deg * pi) / (180)}

# rotation
# rotates x and y through angle
rotate <- function(x, y, angle) {
  xp <- x * cos(angle) - y * sin(angle)
  yp <- y * cos(angle) + x * sin(angle)
  return(list(xp, yp))
}

# positive angles
# returns angles, replacing negative angles with positive ones
pos_angle <- function(angle_deg) {
  angle_deg[angle_deg < 0.0] <- angle_deg[angle_deg < 0.0] + 360.0
  return(angle_deg)
}

# xyr_pos
# finds x- and y-positions and radius from center of angular locations along orbit
# by comparing a specific angle with "e_angle", a set of small angualr increments
# input:  angle along orbit of the point, x-y location, small angles, and offset
# returns radius, x, y, and which e_angle was closest
xyr_pos <- function(eqnx_angle, xe, ye, e_angle, offset) {
  angle_diff <- pos_angle(e_angle) - eqnx_angle
  angle_diff <- pos_angle(angle_diff)
  min_index <- which.min(angle_diff)
  r_eqnx <- sqrt(xe[min_index]^2 + ye[min_index]^2) + offset
  x_eqnx <- r_eqnx * cos(radians(eqnx_angle))
  y_eqnx <- r_eqnx * sin(radians(eqnx_angle))
  return(list(r_eqnx, x_eqnx, y_eqnx, min_index))
}

# Kepler time-of-travel equation (used for debugging and prototyping)
# input:  eccentricity, orbital period, angle along orbit
kepler_time <- function (eccen, T, theta_deg) {
  print(eccen)
  print(T)
  print(theta_deg)
  theta_rad = radians(theta_deg)
  
  E <- 2.0 * atan( sqrt((1.0 - eccen) / (1.0 + eccen)) * tan(theta_rad / 2.0) )
  
  M <- E - eccen * sin(E)
  print(M)
  kepler_time <- ((M / (2 * pi)) * T)
  print(kepler_time)
  if (kepler_time < 0.0) kepler_time <- kepler_time + T
  print(kepler_time)
  return(kepler_time)
}

# set up

# generate angles
delta_1deg <- 1 # 1-degree steps
angle_1deg <-  seq(0, 365, by=delta_1deg) * (360/365)
angle1 <- radians(angle_1deg)
npts_1deg <- length(angle1)

delta_5deg <- 5 # 5-degree steps
angle_5deg<-  seq(0, 365, by=delta_5deg) * (360/365) 
angle5 <- radians(angle_5deg)
npts_5deg <- length(angle5)

# small steps
delta_sm <- 0.2
angle_smdeg <- seq(0, 360, by=delta_sm) * (360/365) 
angle_sm <- radians(angle_smdeg)
npts <- length(angle_smdeg)

# points along the circle
radius <- 1.0
x <- radius * cos(angle5)
y <- radius * sin(angle5)

# eccentricity exaggeration
ecc_factor <- 10

# get orbital parameters
path_prefix <- "e:/" # /Users/bartlein/"
orb_path <- "Projects/Calendar/data/Ellipse/data/"
orb_file <- "orb_elt_150ka_1kyr.csv"
orb_params <- read.csv(paste(path_prefix, orb_path, orb_file, sep=""))
names(orb_params)
nt <- length(orb_params$YearBP)

# get month-length data
monlen_path <- "Projects/Calendar/data/Ellipse/data/"
monlen_file <- "PMIP4_test_cal_noleap_rmonlen.csv"
monlen <- read.csv(paste(path_prefix, monlen_path, monlen_file, sep=""))
names(monlen)

# get day names
day_names_path <- "Projects/Calendar/data/Ellipse/data/"
day_names_file <- "day_names.txt"
day_names <- read.csv(paste(path_prefix, day_names_path, day_names_file, sep=""), as.is=1)
names(day_names)

# output .pngs path
png_path <- paste(path_prefix, "Projects/Calendar/data/Ellipse/orb_0ka/", sep="")

# plot limits
lim <- 1.6
xlim <- c(-1 * lim, lim); ylim <- xlim

# set age
i <- 151 - 0

# define elliptical orbit
ecc <- orb_params$Eccen[i] 
ecc <- ecc * ecc_factor
a <- 1.0
b <- a * sqrt(1.0 - ecc^2)
ae <- a * ecc
perih_angle_deg <- orb_params$Perih_deg[i] + ((80.5/365.0)*360) 
perih_angle <- radians(perih_angle_deg)
# print(c(i, orb_params$YearBP[i], orb_params$Eccen[i], ecc_factor, ecc, a, b, a/b, orb_params$Perih_deg[i], peri_angle_deg))

# generate points along the ellipse -- small increments 
E_anm_small <- radians(angle_smdeg)
theta_small <- 2.0 * atan(sqrt((1.0 + ecc)/(1.0 - ecc)) * tan(E_anm_small/2.0))
r_ellipse <- a * (1.0 - ecc^2) / (1.0 + ecc * cos(theta_small))
xe <- r_ellipse * cos(theta_small) + 2 * ae
ye <- r_ellipse * sin(theta_small) 

# rotate to reflect current longitude of perihelion
rot_xy <- rotate(xe, ye, perih_angle)
xe <- rot_xy[[1]]; ye <- rot_xy[[2]]

# e_angles
e_angle <- degree(atan2(ye,xe))
e_angle <- pos_angle(e_angle)

# equal-step true anomalies -- for position
theta <- 2.0 * atan(sqrt((1.0 + ecc)/(1.0 - ecc)) * tan(angle5/2.0))
rp <- a * (1.0 - ecc^2) / (1.0 + ecc * cos(angle5))
xp <- rp * cos(angle5) + 2 * ae
yp <- rp * sin(angle5) 

# rotate to reflect current longitude of perihelion
rot_xy <- rotate(xp, yp, perih_angle)
xp <- rot_xy[[1]]; yp <- rot_xy[[2]]
p_angle <- pos_angle(degree(atan2(yp, xp) + pi))

# equal-step true anomalies for orbital speed -- no eccentricity exaggeration
theta2 <- 2.0 * atan(sqrt((1.0 + orb_params$Eccen[i])/(1.0 - orb_params$Eccen[i])) * tan(angle5/2.0))
rv <- a * (1.0 - orb_params$Eccen[i]^2) / (1.0 + orb_params$Eccen[i] * cos(angle5))

# perihelion x and y
theta_perih <- pi + radians(1.0/(365.0/360.0))
r_perih =  (a * (1.0 - ecc^2) / (1.0 + ecc * cos(theta_perih))) + 0.1
x_perih <- r_perih * cos(theta_perih) + 2 * ae
y_perih <- r_perih * sin(theta_perih)
rxy <- rotate(x_perih, y_perih, perih_angle)
x_perih <- rxy[[1]]; y_perih <- rxy[[2]]

# equinox and solstice positions
ve_angle <- (80.5 * (360/365) + 180) 
eqnx <- xyr_pos(ve_angle, xe, ye, e_angle, 0.05) 
r_ve <- eqnx[[1]]; x_ve <- eqnx[[2]]; y_ve <- eqnx[[3]]; idx_ve <- eqnx[[4]]

ss_angle <- (173.0 * (360/365) + 180)
eqnx <- xyr_pos(ss_angle, xe, ye, e_angle, 0.05) 
r_ss <- eqnx[[1]]; x_ss <- eqnx[[2]]; y_ss <- eqnx[[3]]; idx_ss <- eqnx[[4]]

ae_angle <- (263.0 * (360/365) + 180) - 360
eqnx <- xyr_pos(ae_angle, xe, ye, e_angle, 0.05) 
r_ae <- eqnx[[1]]; x_ae <- eqnx[[2]]; y_ae <- eqnx[[3]]; idx_ae <- eqnx[[4]]

ws_angle <- (355.5 * (360/365) + 180) - 360
eqnx <- xyr_pos(ws_angle, xe, ye, e_angle, 0.05) 
r_ws <- eqnx[[1]]; x_ws <- eqnx[[2]]; y_ws <- eqnx[[3]]; idx_ws <- eqnx[[4]]

# quadrant-length label positions
VEtoSS_angle <- (80.5 * (360/365) + 180 + 45) 
x_VEtoSS <- 0.3; y_VEtoSS <- -0.2

SStoAE_angle <- (173.0 * (360/365) - 180 + 45)
x_SStoAE <- 0.3; y_SStoAE <- 0.2

AEtoWS_angle <- (263.0 * (360/365) + 180) - 360 + 45
x_AEtoWS <- -0.3; y_AEtoWS <- 0.2

WStoVE_angle <- (355.5 * (360/365) + 180) - 360 + 45
x_WStoVE <- -0.3; y_WStoVE <- -0.2

# main loop

d <- 2 # for debugging
for (d in (1:73)) {

  # png output
  png_name <- sprintf("%03d", d)
  png_file <- paste(png_path, "orb_0ka_", png_name, ".png", sep="")
  png(png_file, width = 2250, height = 2250, units = "px", pointsize = 12, res=300)
  
  jj <- as.integer((74 / delta_1deg) / 2) + d - 1
  if (jj > 73) jj <- jj - 73
  
  jday <- as.integer((p_angle[jj])*(365.0/360.0))
  
  # calculate orbital speed
  mu <- 1.327e20; a_actual <- 1.496e11
  orbital_speed <- sqrt(mu *((2.0/(rv[d]*a_actual)) - (1.0/a_actual)))
  
  # set up for plotting
  oldpar <- par(mar=c(0,0,0,0))
  plot(0, 0, xlim = xlim, ylim = ylim, type = "n", asp = 1,
       axes=FALSE, ann = FALSE, xaxt="n", yaxt="n")

  # draw equinox and solstice lines
  segments(0, 0, x_ve, y_ve, lwd=2, col="gray80")
  text(x_ve, y_ve, "VE", pos=1, cex=2, col="gray60")
  segments(0, 0, x_ss, y_ss, lwd=2, col="gray80")
  text(x_ss, y_ss, "SS", pos=4, cex=2, col="gray60")
  segments(0, 0, x_ae, y_ae, lwd=2, col="gray80")
  text(x_ae, y_ae, "AE", pos=3, cex=2, col="gray60")
  segments(0, 0, x_ws, y_ws, lwd=2, col="gray80")
  text(x_ws, y_ws, "WS", pos=2, cex=2, col="gray60")
  
  # labels
  dtxt <- paste(sprintf("%.2f", monlen$VEtoSS[i]), "d")
  text(x_VEtoSS, y_VEtoSS, dtxt, pos=NULL, cex=, col="gray60")
  dtxt <- paste(sprintf("%.2f", monlen$SStoAE[i]), "d")
  text(x_SStoAE, y_SStoAE, dtxt, pos=NULL, cex=1, col="gray60")
  dtxt <- paste(sprintf("%.2f", monlen$AEtoWS[i]), "d")
  text(x_AEtoWS, y_AEtoWS, dtxt, pos=NULL, cex=1, col="gray60")
  dtxt <- paste(sprintf("%.2f", monlen$WStoVE[i]), "d")
  text(x_WStoVE, y_WStoVE, dtxt, pos=NULL, cex=1, col="gray60")
  
  # # plot points along auxillary circle
  # points(x, y, pch=16, cex=0.5, col="black")
  # text(x, y, as.character(as.integer(angle_deg)), cex=0.5, col="gray70")
  # segments(2, 0, -2, 0, col="gray80")
  
  # # points along elliptical orbit
  # points(xe, ye, pch = 16, cex=0.5, col="red")
  # text(xe, ye, as.character(as.integer(e_angle)), cex=0.8, col="black")
  
  # Earth positions at 5-deg intervals-- gray points
  points(xp, yp, pch=16, cex = 1.5, col = "gray70")
  
  # Earth
  points(xp[jj], yp[jj], pch=16, cex = 3.0, col = "lightblue")
  points(xp[jj], yp[jj], pch=1, cex = 3.0, lwd = 0.5, col = "darkblue")
 
  # Sun
  points(0, 0, pch = 16, col = "yellow3", cex=4)
  points(0, 0, pch = 16, col = "yellow1", cex=3.75)
  points(0, 0, pch = 16, col = "yellow", cex=2)
  # points(0, 0, pch = "+", col = "black")
  
  # perihelion
  points(x_perih, y_perih, pch=16, cex=1.5, col="magenta")
  
  # labels
  age <- sprintf("%3d", -1 * orb_params$YearBP[i] / 1000)
  day_char <- sprintf("%03d", jday)
  title <- paste("Day ", day_char, "  (0 ka ", day_names$day[jday], ")",sep="")
  text((-1*lim), (-1*lim)+.1, title, cex=2, pos=4)
  
  orbital_speed_char <- paste("Orbital speed = ", formatC(orbital_speed, format="f", 
    big.mark=",", digits=0), " m/s", sep="")
  text((-1*lim), (-1*lim)+.3, orbital_speed_char, cex=1.5, pos=4)

  # parameter labels
  xpos <- 1.0; ypos <- -1.15
  ecc_char <- paste("=", sprintf("%.5f", orb_params$Eccen[i]), "(x 10)")
  ecc_title <- expression(italic("e"))
  text(xpos, ypos, ecc_title)
  text(xpos, ypos, ecc_char, pos=4)
  
  obl_char <- paste("=", sprintf("%.2f", orb_params$Obliq_deg[i]), "°")
  obl_title <- expression(epsilon)
  text(xpos, ypos-0.1, obl_title)
  text(xpos, ypos-0.1, obl_char, pos=4)
  
  perih_char <- paste("=", sprintf("%.2f", orb_params$Perih_deg[i]), "°")
  perih_title <- expression(omega)
  text(xpos, ypos-0.2, perih_title)
  text(xpos, ypos-0.2, perih_char, pos=4)
  
  points(xpos, ypos-0.3, pch=16, cex=1.5, col="magenta")
  text(xpos, ypos-0.3, " Perihelion", pos=4)
  
  delta_char <- paste(" =", sprintf("%2d", delta_5deg), "°")
  delta_title <- expression(paste(Delta, nu))
  text(xpos, ypos-0.4, delta_title)
  text(xpos, ypos-0.4, delta_char, pos=4)

dev.off()
}

# ImageMagick 7 command line
# magick -delay 10 @orb_0ka_files.txt -loop 0 -layers Optimize orb_0ka.gif 

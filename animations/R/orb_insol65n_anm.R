# orb_insol65n_anm.R
# insolation long-term mean differences, 65N, 150 ka to 0 ka

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
angle_smdeg <- seq(0, 365, by=delta_sm) * (360/365) 
angle_sm <- radians(angle_smdeg)
npts <- length(angle_smdeg)

# points along the circle
radius <- 1.0
x <- radius * cos(angle5)
y <- radius * sin(angle5)

# eccentricity exaggeration
ecc_factor <- 10

# get orbital parameters
path_prefix <- "/Users/bartlein/" # "e:/" #  
orb_path <- "Projects/Calendar/animations/data/input/"
orb_file <- "orb_elt_150ka_1kyr.csv"
orb_params <- read.csv(paste(path_prefix, orb_path, orb_file, sep=""))
names(orb_params)
nt <- length(orb_params$YearBP)

# get month-length data
monlen_path <- "Projects/Calendar/animations/data/input/"
monlen_file <- "PMIP4_test_cal_noleap_rmonlen.csv"
monlen <- read.csv(paste(path_prefix, monlen_path, monlen_file, sep=""))
names(monlen)

# get insolation data
insol_path <- "Projects/Calendar/animations/data/input/"
insol_file <- "65n150k1k.csv"
insol <- read.csv(paste(path_prefix, insol_path, insol_file, sep=""))
names(insol)

# output .pngs path
png_path <- paste(path_prefix, "Projects/Calendar/animations/data/orb_insol65n_anm/", sep="")

# plot limits
lim <- 1.6
xlim <- c(-1 * lim, lim); ylim <- c((-1 * lim) - 0.15, lim - 0.15)

# color numbers
mcv <- 255
#insol_anm_clr <- c("#000066", "#0000FF", "#00CCFF", "#99FFFF", "#FFFFCC", "#FFCC00", "#FF6633", "#CC3300", "#CCCCCC")
insol_anm_clr <- c("#260053", "#8400C0", "#BE44FF", "#BF80FF", "#FFFFC1", "#FEC30B", "#FC4F09", "#BE2005", "#CCCCCC")
insol_anm_breaks <- c(-400, -50.0, -20.0, -10.0, 0.0, 10.0, 20.0, 50.0, 400)
month_name <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# main loop

i <- 151 - 6
for (i in (1:nt)) {
  
  # png output
  png_name <- sprintf("%03d", -1 * orb_params$YearBP[i] / 1000)
  png_file <- paste(png_path, "orb_insol65n_anm_", png_name, "_ka.png", sep="")
  png(png_file, width = 2250, height = 2250, units = "px", pointsize = 12, res=300)
  
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
  
  # perihelion x and y
  theta_perih <- pi
  r_perih =  (a * (1.0 - ecc^2) / (1.0 + ecc * cos(theta_perih))) + 0.06
  x_perih <- r_perih * cos(theta_perih) + 2 * ae
  y_perih <- r_perih * sin(theta_perih)
  rxy <- rotate(x_perih, y_perih, perih_angle)
  x_perih <- rxy[[1]]; y_perih <- rxy[[2]]

  # equinox and solstice positions
  ve_angle <- (80.5 * (360/365) + 180) 
  eqnx <- xyr_pos(ve_angle, xe, ye, e_angle, 0.06) 
  r_ve <- eqnx[[1]]; x_ve <- eqnx[[2]]; y_ve <- eqnx[[3]]; idx_ve <- eqnx[[4]]
  
  ss_angle <- (173.0 * (360/365) + 180)
  eqnx <- xyr_pos(ss_angle, xe, ye, e_angle, 0.06) 
  r_ss <- eqnx[[1]]; x_ss <- eqnx[[2]]; y_ss <- eqnx[[3]]; idx_ss <- eqnx[[4]]
  
  ae_angle <- (263.0 * (360/365) + 180) - 360
  eqnx <- xyr_pos(ae_angle, xe, ye, e_angle, 0.06) 
  r_ae <- eqnx[[1]]; x_ae <- eqnx[[2]]; y_ae <- eqnx[[3]]; idx_ae <- eqnx[[4]]
  
  ws_angle <- (355.5 * (360/365) + 180) - 360
  eqnx <- xyr_pos(ws_angle, xe, ye, e_angle, 0.06) 
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
  
  oldpar <- par(mar=c(0,0,0,0))
  plot(0, 0, xlim = xlim, ylim = ylim, type = "n", asp = 1, 
       axes=FALSE, ann = FALSE, xaxt="n", yaxt="n")

  segments(0, 0, x_ve, y_ve, lwd=2, col="gray80")
  text(x_ve, y_ve, "VE", pos=1, cex=1.5, col="gray60")
  segments(0, 0, x_ss, y_ss, lwd=2, col="gray80")
  text(x_ss, y_ss, "SS", pos=4, cex=1.5, col="gray60")
  segments(0, 0, x_ae, y_ae, lwd=2, col="gray80")
  text(x_ae, y_ae, "AE", pos=3, cex=1.5, col="gray60")
  segments(0, 0, x_ws, y_ws, lwd=2, col="gray80")
  text(x_ws, y_ws, "WS", pos=2, cex=1.5, col="gray60")
  
  dtxt <- paste(sprintf("%.2f", monlen$VEtoSS[i]), "d")
  text(x_VEtoSS, y_VEtoSS, dtxt, pos=NULL, cex=1, col="gray60")
  dtxt <- paste(sprintf("%.2f", monlen$SStoAE[i]), "d")
  text(x_SStoAE, y_SStoAE, dtxt, pos=NULL, cex=1, col="gray60")
  dtxt <- paste(sprintf("%.2f", monlen$AEtoWS[i]), "d")
  text(x_AEtoWS, y_AEtoWS, dtxt, pos=NULL, cex=1, col="gray60")
  dtxt <- paste(sprintf("%.2f", monlen$WStoVE[i]), "d")
  text(x_WStoVE, y_WStoVE, dtxt, pos=NULL, cex=1, col="gray60")
  
  # month beginning and end and SS difference anomalies
  
  for (m in (1:12)) {
    insol_anm <- insol[i, (2 + m)] -  insol[151, (2 + m)]
    colornum <- cut(insol_anm, insol_anm_breaks, right=TRUE)
    if (i == 151) colornum = 9
    
    start_angle <- ( monlen[i, (26 + m)] * (360/365) + 180) 
    if (start_angle > 360.0) start_angle <- start_angle - 360.0
    eqnx_start_inner <- xyr_pos(start_angle, xe, ye, e_angle, -0.03) 
    r_start_inner <- eqnx_start_inner[[1]]; x_start_inner <- eqnx_start_inner[[2]]
    y_start_inner <- eqnx_start_inner[[3]]; idx_start_inner <- eqnx_start_inner[[4]]
    eqnx_start_outer <- xyr_pos(start_angle, xe, ye, e_angle, 0.03) 
    r_start_outer <- eqnx_start_outer[[1]]; x_start_outer <- eqnx_start_outer[[2]]
    y_start_outer <- eqnx_start_outer[[3]]; idx_start_outer <- eqnx_start_outer[[4]]
    
    end_angle <- (monlen[i, (38 + m)] * (360/365) + 180) 
    if (end_angle > 360.0) end_angle <- end_angle - 360.0
    eqnx_end_inner <- xyr_pos(end_angle, xe, ye, e_angle, -0.03) 
    r_end_inner <- eqnx_end_inner[[1]]; x_end_inner <- eqnx_end_inner[[2]]
    y_end_inner <- eqnx_end_inner[[3]]; idx_end_inner <- eqnx_end_inner[[4]]
    eqnx_end_outer <- xyr_pos(end_angle, xe, ye, e_angle, 0.03) 
    r_end_outer <- eqnx_end_outer[[1]]; x_end_outer <- eqnx_end_outer[[2]]
    y_end_outer <- eqnx_end_outer[[3]]; idx_end_outer <- eqnx_end_outer[[4]]
    
    mid_angle <- (monlen[i, (14 + m)] * (360/365) + 180) 
    if (mid_angle > 360.0) mid_angle <- mid_angle - 360.0
    eqnx_mid <- xyr_pos(mid_angle, xe, ye, e_angle, -0.125) 
    r_mid <- eqnx_mid[[1]]; x_mid <- eqnx_mid[[2]]
    y_mid <- eqnx_mid[[3]]; idx_mid <- eqnx_mid[[4]]
    
    if (idx_end_inner < idx_start_inner) {
      idx_start <- eqnx_start_inner[[4]]
      for (ii in (idx_start_inner:npts)) {
        seg_angle <- (e_angle[ii] )#* (360/365) + 180) 
        eqnx_seg_inner <- xyr_pos(seg_angle, xe, ye, e_angle, -0.03) 
        r_seg_inner <- eqnx_seg_inner[[1]]; x_seg_inner <- eqnx_seg_inner[[2]]
        y_seg_inner <- eqnx_seg_inner[[3]]; idx_seg_inner <- eqnx_seg_inner[[4]]
        eqnx_seg_outer <- xyr_pos(seg_angle, xe, ye, e_angle, 0.03) 
        r_seg_outer <- eqnx_seg_outer[[1]]; x_seg_outer <- eqnx_seg_outer[[2]]
        y_seg_outer <- eqnx_seg_outer[[3]]; idx_seg_outer <- eqnx_seg_outer[[4]]
        segments(x_seg_inner, y_seg_inner, x_seg_outer, y_seg_outer, lwd=2, col=insol_anm_clr[colornum])
      }
      for (ii in (1:idx_end_inner)) {
        seg_angle <- (e_angle[ii] )#* (360/365) + 180) 
        eqnx_seg_inner <- xyr_pos(seg_angle, xe, ye, e_angle, -0.03) 
        r_seg_inner <- eqnx_seg_inner[[1]]; x_seg_inner <- eqnx_seg_inner[[2]]
        y_seg_inner <- eqnx_seg_inner[[3]]; idx_seg_inner <- eqnx_seg_inner[[4]]
        eqnx_seg_outer <- xyr_pos(seg_angle, xe, ye, e_angle, 0.03) 
        r_seg_outer <- eqnx_seg_outer[[1]]; x_seg_outer <- eqnx_seg_outer[[2]]
        y_seg_outer <- eqnx_seg_outer[[3]]; idx_seg_outer <- eqnx_seg_outer[[4]]
        segments(x_seg_inner, y_seg_inner, x_seg_outer, y_seg_outer, lwd=2, col=insol_anm_clr[colornum])
      }
    } else {
      for (ii in (idx_start_inner:idx_end_inner)) {
        seg_angle <- (e_angle[ii] )#* (360/365) + 180) 
        eqnx_seg_inner <- xyr_pos(seg_angle, xe, ye, e_angle, -0.03) 
        r_seg_inner <- eqnx_seg_inner[[1]]; x_seg_inner <- eqnx_seg_inner[[2]]
        y_seg_inner <- eqnx_seg_inner[[3]]; idx_seg_inner <- eqnx_seg_inner[[4]]
        eqnx_seg_outer <- xyr_pos(seg_angle, xe, ye, e_angle, 0.03) 
        r_seg_outer <- eqnx_seg_outer[[1]]; x_seg_outer <- eqnx_seg_outer[[2]]
        y_seg_outer <- eqnx_seg_outer[[3]]; idx_seg_outer <- eqnx_seg_outer[[4]]
        segments(x_seg_inner, y_seg_inner, x_seg_outer, y_seg_outer, lwd=2, col=insol_anm_clr[colornum])
      }
    }
    segments(x_start_inner, y_start_inner, x_start_outer, y_start_outer, lwd=3, col="black")
    segments(x_end_inner, y_end_inner, x_end_outer, y_end_outer, lwd=3, col="black")
    
    text(x_mid, y_mid, month_name[m], pos=NULL, cex=0.9, col="gray60")
  }
  
  # Sun
  points(0, 0, pch = 16, col = "yellow3", cex=4)
  points(0, 0, pch = 16, col = "yellow1", cex=3.75)
  points(0, 0, pch = 16, col = "yellow", cex=2)
  # points(0, 0, pch = "+", col = "black")  
  
  # perihelion
  points(x_perih, y_perih, pch=16, cex=1.5, col="magenta")

  # labels
  age <- sprintf("%3d", -1 * orb_params$YearBP[i] / 1000)
  title <- paste(age, "ka")
  text((-1*lim)+.6, (-1*lim)+.2, title, cex=2, pos=2)
  
  xpos <- -1.45 ; ypos <- -1.82
  ecc_char <- paste("= ", sprintf("%.5f", orb_params$Eccen[i]), " (x ", as.integer(ecc_factor),  ")", sep="")
  ecc_title <- expression(italic("e"))
  text(xpos, ypos, ecc_title)
  text(xpos, ypos, ecc_char, pos=4)
  
  obl_char <- paste("=", sprintf("%.2f", orb_params$Obliq_deg[i]), "°")
  obl_title <- expression(epsilon)
  text(xpos+0.9, ypos, obl_title)
  text(xpos+0.9, ypos, obl_char, pos=4)
  
  perih_char <- paste("=", sprintf("%.2f", orb_params$Perih_deg[i]), "°")
  perih_title <- expression(omega)
  text(xpos+1.5, ypos, perih_title)
  text(xpos+1.5, ypos, perih_char, pos=4)
  
  points(xpos+2.1, ypos, pch=16, cex=1.5, col="magenta")
  text(xpos+2.1, ypos, " Perihelion", pos=4)
  
  # scale
  xorig <- -1.0; yorig <- -1.67; dx <- 0.2; dy <- 0.07
  xleft <- xorig + seq(0, dx*7, by=dx)
  xright <- xleft + dx
  ybottom <- rep(yorig, 8)
  ytop <- ybottom + dy
  rect(xleft, ybottom, xright, ytop, col=insol_anm_clr, lwd=1.5)
  expr <- expression(paste("65°N Mid-month Insolation Differences (", Wm^-2, ")", sep=""))
  text(xleft[5], ytop + dy, expr)
  for (l in (1:7)) {
    text(xleft[1+l], ybottom - dy/2, insol_anm_breaks[1+l], cex=0.8)
  }
  
  dev.off()
  }
  
  # ImageMagick 7 command line
  # magick -delay 50 @orb_insol65n_anm_files.txt -loop 0 -layers Optimize orb_insol65n_anm.gif
  

u2 <-expression(Wm^-2*X)

u2 <-expression(Wm^-2 ~ ~X)
u2 <- expression(paste("Wm^2", "~ ~x", sep-=""))
plot(NULL, NULL, xlim=c(-1,1), ylim=c(-1,1))
text(0, 0, u2)


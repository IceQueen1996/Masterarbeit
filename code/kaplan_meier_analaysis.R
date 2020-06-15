library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(survival)
library(survminer)

# read AE data
ae <- read.csv("C:\\Users\\Elly\\Dropbox\\Biotronik\\data\\200603 C1503 BIOSOLVE-IV - csv datasets\\AE.csv")
enrol <- read.csv("C:\\Users\\Elly\\Dropbox\\Biotronik\\data\\200603 C1503 BIOSOLVE-IV - csv datasets\\ENROL.csv")

# create needed functions
descrption_table <- function(patient_count, survFit){
  time <- c(0,180,365,540,730)
  # get indices
  tmp <- if (length(which(survFit$time < 1)) != 0) max(which(survFit$time < 1)) else 0
  tmp2 <- if (length(which(survFit$time < 181)) != 0) max(which(survFit$time < 181)) else 0
  tmp3 <- if (length(which(survFit$time < 366)) != 0) max(which(survFit$time < 366)) else 0
  tmp4 <- if (length(which(survFit$time < 541)) != 0) max(which(survFit$time < 541)) else 0
  tmp5 <- if (length(which(survFit$time < 731)) != 0) max(which(survFit$time < 731)) else 0

  patients_at_risk = c(patient_count, survFit$n.risk[tmp2],
                       survFit$n.risk[tmp3],
                       survFit$n.risk[tmp4],
                       survFit$n.risk[tmp5])
  events <- c(sum(survFit$n.event[1:tmp]),
              sum(survFit$n.event[1:tmp2]),
              sum(survFit$n.event[1:tmp3]),
              sum(survFit$n.event[1:tmp4]),
              sum(survFit$n.event[1:tmp5]))
  anteil <- 100*(events/patient_count)
  lower_bound <- round((1-c(if (length(survFit$upper[tmp]) != 0) survFit$upper[tmp] else 0,
                            if (length(survFit$upper[tmp2]) != 0) survFit$upper[tmp2] else 0,
                            if (length(survFit$upper[tmp3]) != 0) survFit$upper[tmp3] else 0,
                            if (length(survFit$upper[tmp4]) != 0) survFit$upper[tmp4] else 0,
                            if (length(survFit$upper[tmp5]) != 0) survFit$upper[tmp5] else 0))*100,2)

  upper_bound <- round((1-c(if (length(survFit$lower[tmp]) != 0) survFit$lower[tmp] else 0,
                            if (length(survFit$lower[tmp2]) != 0) survFit$lower[tmp2] else 0,
                            if (length(survFit$lower[tmp3]) != 0) survFit$lower[tmp3] else 0,
                            if (length(survFit$lower[tmp4]) != 0) survFit$lower[tmp4] else 0,
                            if (length(survFit$lower[tmp5]) != 0) survFit$lower[tmp5] else 0))*100,2)
  descr <- as.data.frame(t(data.frame(patients_at_risk, events, paste0(round(100*(events/patient_count),2)," [",lower_bound, ", ", upper_bound, "]", sep=""))))
  colnames(descr) <- c(paste0("day ", time, sep=""))
  rownames(descr) <- c("Patients at risk", "n", "Failure (estimate)")
  return(descr)
}

render_survPlot <- function(patient_count, survFit, ylab, ylim, yPos_of_line_descr, outPath){
  jpeg(file=outPath, width=746, height=412)
  
  # get description dataframe
  description <- descrption_table(patient_count, survFit = survFit)

  tlr_plot <- ggsurvplot(survFit, fun = "event", conf.int = T, risk.table = FALSE,
                         palette = "blue", legend = "none", xlab = "Time to Event [days]", ylab = ylab,
                         surv.scale="percent", ylim = ylim, xlim = c(1,731))

  tlr_plot <- tlr_plot$plot +
    scale_x_discrete(breaks = c(0,180,365,540,730), labels = c(0,180,365,540,730), limits=c(0:731)) +
    annotate("text", x=610, y=yPos_of_line_descr, size = 4,
             label= paste(
               round((sum(survFit$n.event[1:max(which(survFit$time < 730))])/patient_count)*100,2),
               " [95%CI: ",
               round((1-survFit$upper[max(which(survFit$time < 730))])*100,2),
               ", ",
               round((1-survFit$lower[max(which(survFit$time < 730))])*100, 2),
               "]", sep=""))

  # Add table below the graph
  tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),
                       base_size = 10,
                       padding = unit(c(2, 4), "mm"))
  grid.arrange(
    tlr_plot,
    tableGrob(description, theme=tt),
    nrow=2, as.table = TRUE, heights = c(2, 1),
    vp=viewport(width=0.95, height=1))
  
  dev.off()
}

############################### TV_MI ################################
# filter for death and TV_MI
idx <- ae %>% filter(((aemiehdef == "Yes" | aemiscaidef == "Yes") & aemibytv == "Yes") | aedth == "Yes")

# set up variables for surv analysis
TV_MI <- seq(0,0,length.out=length(enrol$subjid))
TV_MI[which(idx$aemibytv == "Yes")] <- 1

aedays <- c(idx$aedays, seq(max(idx$aedays),max(idx$aedays),length.out=length(enrol$subjid)-length(idx$aedays)))

# survival analysis for TV_MI
survFit_tv_mi <- survfit(Surv(aedays, TV_MI) ~ 1, data = data.frame(aedays, TV_MI))

############################### clinically driven TLR ################################
# filter for death and TLR
idx <- ae %>% filter(aedth == "Yes" |
                       (aerevasc == "Target lesion revascularization" &
                          aerevasccld == "Yes"))

# set up variables for surv analysis
TLR <- seq(0,0,length.out=length(enrol$subjid))
TLR[which(idx$aerevasccld == "Yes")] <- 1

aedays <- c(idx$aedays, seq(max(idx$aedays),max(idx$aedays),length.out=length(enrol$subjid)-length(idx$aedays)))

# survival analysis for TLR
survFit_TLR <- survfit(Surv(aedays, TLR) ~ 1, data = data.frame(aedays, TLR))

############################### CABG ################################
# filter for death and CABG
idx <- ae %>% filter(aedth == "Yes" |
                       aecabg == "Yes")

# set up variables for surv analysis
cabg <- seq(0,0,length.out=length(enrol$subjid))
cabg[which(idx$aecabg == "Yes")] <- 1

aedays <- c(idx$aedays, seq(max(idx$aedays),max(idx$aedays),length.out=length(enrol$subjid)-length(idx$aedays)))

# survival analysis for cardiac death
survFit_cabg <- survfit(Surv(aedays, cabg) ~ 1, data.frame(aedays, cabg))


############################### cardiac death ################################
# filter for death and cardiac death
idx <- ae %>% filter(aedth == "Yes")

# set up variables for surv analysis
cardiacDeath <- seq(0,0,length.out=length(enrol$subjid))
cardiacDeath[which(idx$aedthcs == "Cardiac")] <- 1

aedays <- c(idx$aedays, seq(max(idx$aedays),max(idx$aedays),length.out=length(enrol$subjid)-length(idx$aedays)))

# survival analysis for cardiac death
survFit_carDeath <- survfit(Surv(aedays, cardiacDeath) ~ 1, data.frame(aedays, cardiacDeath))



############################### target lesion failure (TLF) ################################
# filter for death and TLF
idx <- ae %>% filter(aedth == "Yes" |
                       ((aemiehdef == "Yes" | aemiscaidef == "Yes") & aemibytv == "Yes") |
                       (aerevasc == "Target lesion revascularization" &
                                         aerevasccld == "Yes") |
                       aecabg == "Yes")

# set up variables for surv analysis
TLF <- seq(0,0,length.out=length(enrol$subjid))
TLF[which(idx$aedthcs == "Cardiac")] <- 1
TLF[which(idx$aemibytv == "Yes")] <- 1
TLF[which(idx$aerevasccld == "Yes")] <- 1
TLF[which(idx$aecabg == "Yes")] <-1
sum(TLF)

aedays <- c(idx$aedays, seq(max(idx$aedays),max(idx$aedays),length.out=length(enrol$subjid)-length(idx$aedays)))

# survival analysis for TLF
survFit_TLF<- survfit(Surv(aedays, TLF) ~ 1, data = data.frame(aedays, TLF))

############################## Plot surival data ###########################################
# plot TLR
render_survPlot(patient_count = length(enrol$subjid), survFit = survFit_TLR, ylab = "TLR Rate", ylim = c(0,0.15), yPos_of_line_descr = 0.07, outPath = "C:\\Users\\Elly\\Uni\\Masterarbeit\\plots\\TLR_KMP.jpeg")

# plot cardiac death
render_survPlot(patient_count = length(enrol$subjid), survFit = survFit_carDeath, ylab = "Cardiac Death", ylim = c(0,0.15), yPos_of_line_descr = 0.015, outPath = "C:\\Users\\Elly\\Uni\\Masterarbeit\\plots\\CD_KMP.jpeg")

# plot TV_MI
render_survPlot(patient_count = length(enrol$subjid), survFit = survFit_tv_mi, ylab = "TV-MI Rate", ylim = c(0,0.15), yPos_of_line_descr = 0.025, outPath = "C:\\Users\\Elly\\Uni\\Masterarbeit\\plots\\TV_MI_KMP.jpeg")

# plot CABG
render_survPlot(patient_count = length(enrol$subjid), survFit = survFit_cabg, ylab = "CABG Rate", ylim = c(0,0.15), yPos_of_line_descr = 0.025, outPath = "C:\\Users\\Elly\\Uni\\Masterarbeit\\plots\\CABG_KMP.jpeg")


# plot TLF
render_survPlot(patient_count = length(enrol$subjid), survFit = survFit_TLF, ylab = "TLF Rate", ylim = c(0,0.15), yPos_of_line_descr = 0.085, outPath = "C:\\Users\\Elly\\Uni\\Masterarbeit\\plots\\TLF_KMP.jpeg")



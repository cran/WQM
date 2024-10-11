## ----include = FALSE----------------------------------------------------------
#------------------------------------------------------------------------------#
# global parameters----
op <- par() # save original par

# WBC parameters
wavelet <- "morlet"

dt <- 1
dj <- 1

num.sim <- 5
theta <- switch(1, 0.1, 0.2, 0.5, 1) # mm/h or mm/d
QM <- switch(3, "MBCr","MRS","QDM") # MBC include MBCp, MBCr, MBCn
variable <- switch(2, "sine", "prep")

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", 
  
  out.width = "85%", fig.width = 9, fig.height = 9, dpi=300, 
  fig.align = "center", fig.pos = "h!"
)

## ----setup--------------------------------------------------------------------
library(WQM)

library(MBC)
library(WaveletComp)

library(tidyr)
library(dplyr)
library(data.table)
library(scales)
library(ggplot2)

## ----dat----------------------------------------------------------------------
data(NWP.rain)
summary(NWP.rain)

station.no <- levels(NWP.rain$Station)

Ind.samp <- c(34,75)[1]; Ind.samp
station.samp <- station.no[Ind.samp]; station.samp
station.samp

NWP.rain.h <- NWP.rain %>% na.omit() %>% subset(Station==station.samp)
summary(NWP.rain.h)

#obs overview
plot(NWP.rain.h$obs, type="l",col=1, ylab="P", xlab=NA)
lines(NWP.rain.h$mod, type="l", col=2)

## ----cwt----------------------------------------------------------------------
folds <- cut(seq(1,nrow(NWP.rain.h)),breaks=2,labels=FALSE)
subset <- which(folds==1, arr.ind=TRUE)

variable <- "prep"
PRand <- TRUE
seed <- 2021-11-28


###===============================###===============================###
if(TRUE){
  data <- list()
  l <- 1
  data[[l]] <- NWP.rain.h
  #cat("Station:", l)

  ### list for storing results
  data_sim <- list()
  data[[l]]$index <- as.numeric(format(data[[l]]$Date,format='%j'))

  summary(data[[l]]$obs[subset])
  wt_o <- t(WaveletComp::WaveletTransform(x=data[[l]]$obs[subset],dt=dt,dj=dj)$Wave)
  wt_m <- t(WaveletComp::WaveletTransform(x=data[[l]]$mod[subset],dt=dt,dj=dj)$Wave)
  wt_p <- t(WaveletComp::WaveletTransform(x=data[[l]]$mod[-subset],dt=dt,dj=dj)$Wave)

  if(l==1) print(paste0("No of decomposition level: ", ncol(wt_o)))

  ### return CWT coefficients as a complex matrix with rows and columns representing times and scales, respectively.
  wt_o_mat <- as.matrix(wt_o)
  real.o <- Re(wt_o_mat)
  modulus.o <- Mod(wt_o_mat)       ### derive modulus of complex numbers (radius)
  phases.o <- Arg(wt_o_mat)        ### extract phases (argument)

  apply(modulus.o, 2, var)

  wt_m_mat <- as.matrix(wt_m)
  real.m <- Re(wt_m_mat)
  modulus.m <- Mod(wt_m_mat)
  phases.m <- Arg(wt_m_mat)

  wt_p_mat <- as.matrix(wt_p)
  real.p <- Re(wt_p_mat)
  modulus.p <- Mod(wt_p_mat)
  phases.p <- Arg(wt_p_mat)
}


## ----fig----------------------------------------------------------------------
# CWT plot----
# selected decomposition level to plot
level.samp <- c(seq(1,ncol(modulus.o), by=1));level.samp

df.cwt <- data.frame(No=data[[l]]$Date[subset],
                     obs=data[[l]]$obs[subset],modulus.o) %>%
  gather(Group, Value,2:(ncol(modulus.o)+2))

df.cwt_sub <- subset(df.cwt, Group %in% c("obs",paste0("X",level.samp)))
summary(factor(df.cwt_sub$Group))

##observations----
p0 <- ggplot(subset(df.cwt_sub,Group=="obs"), aes(x=No, y=Value)) +
  geom_line()+
  #facet_grid(Group~., switch="both") +
  #facet_wrap(Group~., strip.position = "left", ncol=1, labeller = label_parsed) +
  scale_y_continuous(breaks = pretty_breaks(n = 4)) +
  scale_x_datetime(date_breaks="3 month", date_labels ="%Y-%m", expand = c(0,0)) +
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  labs(y="Precipitation(mm/h)") +
  theme_bw() +
  theme(text = element_text(size = 16),
        plot.margin = unit(c(0.5,1,0.5, 1), "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),

        strip.text.y = element_text(size=16),
        # axis.text.y = element_text(angle = 90, hjust=0.5, size=12),
        # axis.title.x = element_text(size=16),
        # axis.title.y = element_text(size=16),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(),
        #axis.title.y = element_blank()
  )
p0

##amplitudes of observations----
df.cwt_sub$Group <- factor(df.cwt_sub$Group, labels=c("obs",paste('italic(s)[',level.samp-1,']',sep = '')))

p1 <- ggplot(subset(df.cwt_sub,Group!="obs"), aes(x=No, y=Value)) +
      geom_line()+
      #facet_grid(Group~., switch="both") +
      facet_wrap(Group~., strip.position = "left", ncol=1, labeller = label_parsed) +
      scale_y_continuous(breaks = pretty_breaks(n = 4)) +
      scale_x_datetime(date_breaks="3 month", date_labels ="%Y-%m",expand = c(0,0)) +
      #scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
      labs(y="Amplitudes", color="temp") +
      theme_bw() +
      theme(text = element_text(size = 16),
            plot.margin = unit(c(0.2,0.4,0.5, 0), "cm"),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),

            strip.text.y = element_text(size=16),
            # axis.text.y = element_text(angle = 90, hjust=0.5, size=12),
            # axis.title.x = element_text(size=16),
            # axis.title.y = element_text(size=16),
            strip.background = element_blank(),
            strip.placement = "outside",
            axis.title.x = element_blank(),
            axis.title.y = element_text(margin = margin(t = 0, r = -2, b = 0, l = 0))
            #axis.title.y = element_blank()
      )
p1 %>% print()

##phases of observations----
df.cwt <- data.frame(No=data[[l]]$Date[subset],
                     obs=data[[l]]$obs[subset],phases.o) %>%
  gather(Group, Value,2:(ncol(modulus.o)+2))

df.cwt_sub <- subset(df.cwt, Group %in% c("obs",paste0("X",level.samp)))
summary(factor(df.cwt_sub$Group))

df.cwt_sub$Group <- factor(df.cwt_sub$Group, labels=c("obs",paste('italic(s)[',level.samp-1,']',sep = '')))

p2 <- ggplot(subset(df.cwt_sub,Group!="obs"), aes(x=No, y=Value)) +
      geom_line()+
      #facet_grid(Group~., switch="both") +
      facet_wrap(Group~., strip.position = "left", ncol=1, labeller = label_parsed) +
      scale_y_continuous(limits=c(-pi,pi),breaks = pretty_breaks(n = 4)) +
      #scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
      scale_x_datetime(date_breaks="3 month", date_labels ="%Y-%m",expand = c(0,0)) +
      labs(y="Phases",color="temp") +
      theme_bw() +
      theme(text = element_text(size = 16),
            plot.margin = unit(c(0.2,0.4,0.5, 0), "cm"),
            panel.border = element_rect(color = 'black'),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),

            strip.text.y = element_text(size=16),

            # axis.text.y = element_text(angle = 90, hjust=0.5, size=12),
            # axis.title.x = element_text(size=16),
            # axis.title.y = element_text(size=16),
            strip.background = element_blank(),
            strip.placement = "outside",
            axis.title.x = element_blank(),
            axis.title.y = element_text(margin = margin(t = 0, r = -2, b = 0, l = 0))
            #axis.title.y = element_blank()
      )
p2 %>% print()


## -----------------------------------------------------------------------------
if(QM %like% "MBC"){
  modulus.tmp <- do.call(QM,list(o.c=modulus.o, m.c=modulus.m,
                                 m.p=modulus.p, ratio.seq=rep(TRUE, ncol(modulus.m)), #trace=TRUE,
                                 silent=TRUE,n.tau=100
                                 ))
  modulus.bcc <- modulus.tmp$mhat.c
  modulus.bcf <- modulus.tmp$mhat.p
} else if(QM=="MRS") {
  modulus.tmp <- do.call(QM,list(o.c=modulus.o, m.c=modulus.m, m.p=modulus.p))
  modulus.bcc <- modulus.tmp$mhat.c
  modulus.bcf <- modulus.tmp$mhat.p
} else if(QM=="QDM") {
  #cat("QDM with ratio=T \n")
  modulus.tmp <- lapply(1:ncol(modulus.o), function(i)
    QDM(o.c=modulus.o[,i], m.c=modulus.m[,i], m.p=modulus.p[,i], ratio=FALSE))
  modulus.bcc <- sapply(modulus.tmp, function(ls) ls$mhat.c)
  modulus.bcf <- sapply(modulus.tmp, function(ls) ls$mhat.p)
}

sum(modulus.bcc<0)
sum(modulus.bcf<0)

modulus.bcf[modulus.bcf<0] <-0
modulus.bcc[modulus.bcc<0] <-0

phases.tmp <- lapply(1:ncol(phases.o), function(i)
  QDM(o.c=phases.o[,i], m.c=phases.m[,i], m.p=phases.p[,i]))
phases.bcc <- sapply(phases.tmp, function(ls) ls$mhat.c)
phases.bcf <- sapply(phases.tmp, function(ls) ls$mhat.p)

## modulus----
df.modulus <- rbind(data.frame(mod="obs",no=data[[l]]$Date[subset],modulus.o),
                      data.frame(mod="cur",no=data[[l]]$Date[subset],modulus.m),
                      data.frame(mod="bcc",no=data[[l]]$Date[subset],modulus.bcc)) %>%
    gather(lev, amp, 3:((ncol(modulus.o)+2))) #%>% mutate(amp=as.numeric(amp), subset=as.numeric(subset))

summary(df.modulus$mod)
df.modulus$mod <- factor(df.modulus$mod, levels = c("obs","cur","bcc"),
                         labels = c("obs"="Observed","mod"="Raw","QM"))

df.modulus_sub <- subset(df.modulus, lev %in% c(paste0("X",level.samp)))
summary(factor(df.modulus_sub$lev))

df.modulus_sub$lev <- factor(df.modulus_sub$lev, labels=c(paste('italic(s)[',level.samp-1,']',sep = '')))

p3 <- ggplot(df.modulus_sub, aes(x=no, y=amp, color=mod)) +
      geom_line(linewidth=0.5)+
      #facet_grid(Group~., switch="both") +
      facet_wrap(lev~., strip.position = "left", ncol=1, labeller = label_parsed) +
      scale_color_manual(values=c("black","red","blue"))+
      #scale_y_continuous(limits=c(-pi,pi),breaks = scales::pretty_breaks(n = 3)) +
      #scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
      scale_x_datetime(date_breaks="3 month", date_labels ="%Y-%m", expand = c(0,0)) +
      labs(color=NULL, x="Date") +
      theme_bw() +
      theme(text = element_text(size = 12),
            plot.margin = unit(c(0.2,0.1,0.5, 0), "cm"),
            panel.border = element_rect(color = 'black'),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),

            strip.text.y = element_text(size=16),

            # axis.text.y = element_text(angle = 90, hjust=0.5, size=12),
            # axis.title.x = element_text(size=16),
            # axis.title.y = element_text(size=16),
            #legend.position = 'none',
            legend.position = c(0.68,0.98),
            legend.direction = "horizontal",
            legend.background = element_rect(fill = "transparent"),

            strip.background = element_blank(),
            strip.placement = "outside",
            #axis.title.x = element_blank(),
            axis.title.y = element_blank()
      )

p3

summary(df.modulus_sub)

p4 <- ggplot(df.modulus_sub) +
        stat_ecdf(aes(x=amp, color=mod, size=mod), geom = "step",pad = TRUE) +
        facet_wrap(lev~., strip.position = "left", ncol=1, labeller = label_parsed)+

        scale_x_continuous(trans="log2",limits = c(1/2^10,16), breaks=c(0.001,0.05,1,6)) +
        #scale_x_continuous(trans="sqrt",limits = c(0,4), breaks=c(0,1,2,4)) +

        scale_color_manual(values=c("black","red","blue"))+
        scale_size_manual(values=c(1,0.5,0.5))+
        #guides(color=guide_legend(override.aes = list(size=5))) +
        labs(color=NULL, size=NULL, x="log2") +
        theme_bw() +
        theme(text = element_text(size = 12),
              plot.margin = unit(c(0.2,0.1,0.5, 0), "cm"),
              panel.border = element_rect(color = 'black'),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),

              strip.text.y = element_text(size=16),
              # axis.text.y = element_text(angle = 90, hjust=0.5, size=12),
              # axis.title.x = element_text(size=16),
              # axis.title.y = element_text(size=16),
              legend.position = 'none',

              strip.background = element_blank(),
              strip.placement = "outside",
              #axis.title.x = element_blank(),
              axis.title.y = element_blank()
        )
p4


## ----Prand--------------------------------------------------------------------
# Phase random - validation same for calibration
noise_mat <- list()
for (r in 1:num.sim) {

  data.obs <- as.vector(sapply(data, function(ls) ls$obs[subset]))
  ts_wn <- sample(data.obs, size=length(data[[1]]$obs[-subset]), replace=T)

  wt_noise <- t(WaveletComp::WaveletTransform(x=ts_wn,dt=dt,dj=dj)$Wave)

  noise_mat[[r]] <- as.matrix(wt_noise)
}

phases.rand <- lapply(1:num.sim, function(r)  {
  tmp <- Arg(noise_mat[[r]])

  ord.phases <- apply(phases.p, 2, order)
  #ord.phases <- apply(phases.p, 2, order)

  tmp.rank <- apply(tmp, 2, sort)
  tmp.n <- sapply(1:ncol(tmp), function(ii) {tmp[ord.phases[,ii],ii] <- tmp.rank[,ii];
  return(tmp[,ii])})

  return(tmp.n)
})

## bcf----
mat_new_val <- matrix(complex(modulus=modulus.bcf,argument=phases.p),ncol=ncol(phases.p))
rec_val <- fun_icwt(x.wave=mat_new_val,dt=dt,dj=dj)
if(variable=="prep") rec_val[rec_val<=theta] <- 0

### apply wavelet reconstruction to randomized signal----
mat_val_r <- prsim(modulus.bcf, phases.p, noise_mat)

data_sim_val <- sapply(1:num.sim, function(r) fun_icwt(x.wave=mat_val_r[[r]], dt=dt, dj=dj))
if(variable=="prep") data_sim_val[data_sim_val<=theta] <- 0
colnames(data_sim_val) <- paste0("r",seq(1:num.sim))

data.val <- data[[l]][-subset,]
data.val$bcf <-  rec_val
data.val <- data.frame(data.val, data_sim_val)

summary(data.val)

## ----qm-----------------------------------------------------------------------
# QM approach-------------------------------------------------------------------
for(k_fold in 1){

  data.o <- data[[l]]$obs[subset];
  data.m <- data[[l]]$mod[subset];
  data.p <- data[[l]]$mod[-subset];

  data.tmp <- QDM(o.c=data.o, m.c=data.m, m.p=data.p, ratio=TRUE)#, ratio.max = 1, trace=theta))
  data.bcc <- data.tmp$mhat.c
  data.bcf <- data.tmp$mhat.p

}

## ----comp---------------------------------------------------------------------
df_fut <- cbind(data.val[,c('Date',"obs","mod","bcf",paste0("r",1:num.sim))],qm=data.bcf) %>%
  gather(Group, Value,2:(4+num.sim+1))
summary(factor(df_fut$Group))

summary(df_fut)
p5 <- ggplot(data=df_fut,aes(x=Value,color = Group, size=Group)) +
  stat_ecdf(aes(color = Group, size=Group), geom = "step", pad = FALSE) +
  #stat_density(geom='line', position='identity') +
  #facet_grid(W~Station, scales = "free") + #, labeller = labeller(Grid = Predictor.labs)) +
  # scale_color_manual(values=c("black","red","green","blue"))+#,
  # scale_size_manual(values=c(1,1,0.5,0.5))+
  scale_color_manual(values = c("Black",alpha("Red",0.6),alpha("green",0.6),
                                alpha("blue",0.6), rep('grey',5))) +
  scale_size_manual(values=c(1,1,0.5,0.5, rep(0.1,5))) +

  #labels=c("obs","mod.c","mod.bcc")) +

  scale_x_continuous(limits=c(0, 20), breaks = pretty_breaks(n = 3)) +
  #labs(color=paste0(variable.ncep[i], "-level: ",level[j])) +
  guides(color=guide_legend(override.aes = list(alpha=1)),
         size=guide_legend(override.aes = list(size=1))) +
  theme_bw() +
  theme(text = element_text(size = 16),
        plot.margin = unit(c(1,1,0.5, 1), "cm"),
        # panel.grid.minor = element_blank(),
        # panel.grid.major = element_blank(),
        # axis.text.y = element_text(angle = 90, hjust=0.5, size=12),
        # axis.title.x = element_text(size=16),
        # axis.title.y = element_text(size=16),
        legend.title=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
  )

p5


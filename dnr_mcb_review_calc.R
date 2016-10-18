library(ggplot2, quietly = TRUE)
library(stringr, quietly = TRUE)
#library(plyr, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(data.table)
library(base)

setwd("D:\\Cloud\\OneDrive\\!DNR\\mcbm_review")
rrca.path <- "E:\\dnr\\RRCA\\MCB_Review\\RRCA_RUNS\\"
mcb.path <- "E:\\dnr\\RRCA\\MCB_Review\\MCB_ModelFiles_ModelsOnly20160714\\ModelFiles20160714\\MCB_Model\\working_TR\\"
acre2foot = 43560.0


# groundwater level analysis -----------------------------------------------------------------------
read.gwh = function (file) {
  gwh = readLines(file)
  gwh[1] = str_sub(gwh[1], str_locate(gwh[1], "TIME")[1], -1)
  sim = read.table(text = gwh, header = TRUE)
  return(sim)
}

RRCA.readhyd <- function(filename){
  to.read = file(filename, "rb")
  
  num_before = readBin(to.read, integer(), n = 1, size = 4)
  NUMH = readBin(to.read, integer(), n = 1, size = 4)
  ITMUNI = readBin(to.read, integer(), n = 1, size = 4)
  num_after = readBin(to.read, integer(), n = 1, size = 4)
  
  num_before = readBin(to.read, integer(), n = 1, size = 4)
  data <- read.csv( text = paste0(readChar(to.read, c(4,rep(20,NUMH))), collapse = ",") )
  num_after = readBin(to.read, integer(), n = 1, size = 4)
  
  m = 0
  num_before = readBin(to.read, integer(), n = 1, size = 4)
  l <- readBin(to.read, numeric(), n = NUMH + 1, size = 8)
  num_after = readBin(to.read, integer(), n = 1, size = 4)
  while(length(l)>0) {
    m = m + 1
    data[m,] <- l
    num_before = readBin(to.read, integer(), n = 1, size = 4)
    l <- readBin(to.read, numeric(), n = NUMH + 1, size = 8)
    num_after = readBin(to.read, integer(), n = 1, size = 4)
  }
  
  
  return(data) 
}


### read observation ###
obs <- read.csv("ob.csv", stringsAsFactors = FALSE)
obs$Well <- str_replace_all(str_replace_all(str_replace_all(obs$Well, " ", ""), "-", ""), "#", "")
idx_long_name <- str_length(obs$Well) > 6
obs$Well[idx_long_name] <- str_c(str_sub(obs$Well[idx_long_name], 1, 1), str_sub(obs$Well[idx_long_name], -5, -1))
obs$Date = as.Date(obs$Date, "%m/%e/%Y")

#### mcb gw results #### 
mcb.sim = read.gwh(paste0(mcb.path,"mcb.hydro.gwh"))
head(mcb.sim)
sim.names <- names(mcb.sim) 
sim.names <- str_replace(sim.names, "HDC002", "") 
names(mcb.sim) <- sim.names

mcb.sim$TIME = as.Date("1949-12-31") + mcb.sim$TIME
idx = match(obs$Date, mcb.sim$TIME)

obs$mcb_sim = 0.0
for (i in 1:nrow(obs)) {
  obs$mcb_sim[i] = mcb.sim[idx[i],obs$Well[i]]
}

#### rrca gw results ####
rrca.sim = RRCA.readhyd(paste0(rrca.path,"2000.hyd.result"))
rrca.sim$TIME = as.Date("1950-12-31") + rrca.sim$TIME / 86400
for (i in 2001:2013) {
  gwh = RRCA.readhyd(paste0(rrca.path,i, ".hyd.result"))
  gwh = gwh[-1, ]
  gwh$TIME = as.Date(paste0(i-1,"-12-31")) + gwh$TIME / 86400
  rrca.sim = rbind(rrca.sim, gwh)
}

sim.names <- names(rrca.sim) 
sim.names <- str_replace(sim.names, "HDC001", "") 
names(rrca.sim) <- sim.names
tail(rrca.sim, 30)

idx = match(obs$Date, rrca.sim$TIME)
obs$Date[is.na(idx)]
obs$rrca_sim = 0.0
for (i in 1:nrow(obs)) {
  obs$rrca_sim[i] = rrca.sim[idx[i],obs$Well[i]]
}

#### stats of gw results ####


by_obswel <- group_by(obs, Well)
obs.summary <- arrange(
  summarise(by_obswel,
            N = n(),
            R2_MCB =  cor(Ob, mcb_sim)^2,
            AME_MCB = sum(abs(Ob-mcb_sim))/N, 
            R2_RRCA =  cor(Ob, rrca_sim)^2,
            AME_RRCA = sum(abs(Ob-rrca_sim))/N), 
  Well)

write.csv(file = "obs_summary.csv",obs.summary)

#### plots of gw results ####
gwlevel.mcb <- mutate(melt(mcb.sim, id.vars="TIME", variable.name="Well"), Model="MCB")
gwlevel.rrca <- mutate(melt(rrca.sim, id.vars="TIME", variable.name="Well"), Model="RRCA")
gwlevel.obs <- mutate(select(obs, TIME=Date, Well, value=Ob), Model="OBS")
gwlevel <- rbind(gwlevel.obs, gwlevel.mcb, gwlevel.rrca) 

g <- ggplot() + #data=gwlevel, aes(x=TIME)
  geom_point(data=filter(gwlevel, Model=="OBS"), aes(x=TIME,y=value, color=Model), size=0.3) + 
  geom_line(data=filter(gwlevel, Model=="MCB"), aes(x=TIME,y=value, color=Model)) +
  geom_line(data=filter(gwlevel, Model=="RRCA"), aes(x=TIME,y=value, color=Model)) +
  ylab("Groundwater Level, (ft)") + xlab("Time") + 
  facet_wrap(~Well, ncol=5, scales="free_y") + theme(legend.position="none")
scale=1.5
ggsave("gw_level.png",g,dpi=600,width=6*scale, height=8*scale,units="in")


# plot one by one

# budget analysis -----------------------------------------------------------------------

readBud <- function(fname, itype, zn){
  #read basline first
  bs <- fread(fname, sep = ",", stringsAsFactors = FALSE, data.table = FALSE, check.names = TRUE)
  #head(bs)
  #print(fname)
  
  ntype <- match("From Other Zones", names(bs)) - match("ZONE", names(bs)) - 1
  if ( missing(itype)) itype <- 1:ntype
  itype <- itype[itype >=1 & itype <= ntype]
  
  if ( ! missing(zn) ) bs$ZONE <- zn[bs$ZONE]
  
  # remove zone number is zero
  bs <- filter(bs, ZONE!=0)
  
  #get the type you want
  f<- bs[,4+itype] - bs[,4+itype+ntype+2]
  
  #write.csv(bs,"text.csv")
  #aggregate by zone and time step
  f.out <- aggregate(f, list(ZONE = bs$ZONE, STEP = bs$STEP, PERIOD = bs$PERIOD, TIME = bs$TOTIM), sum)
  
  #f.out$Step <- 1:nrow(f.out)
  return(f.out)
}

numberOfDays <- function(date) {
  m <- format(date, format="%m")
  
  while (format(date, format="%m") == m) {
    date <- date + 1
  }
  
  return(as.integer(format(date - 1, format="%d")))
}

ndays <- function(d) {
  last_days <- 28:31
  rev(last_days[which(!is.na( 
    as.Date( paste( substr(d, 1, 8), 
                    last_days, sep = ''), 
             '%Y-%m-%d')))])[1]
}

#system.time(sapply(rep(rdate,10), FUN = ndays)) 
#user  system elapsed 
#1.03    0.11    1.14 
#system.time(sapply(rep(rdate,10), FUN = numberOfDays))
#user  system elapsed 
#8.15    0.01    8.16



#### rrca budget ####

rrca.zon = scan(paste0(rrca.path, "RRCA.zon"), skip=2)
rrca.zon.sum = summarise(group_by(data.frame(zone = rrca.zon), zone), n=n())
rrca.zon.n = sum(rrca.zon.sum$n[rrca.zon.sum$zone %in% c(1,3)])
rrca.area = rrca.zon.n * 5280 ^ 2

rzone = c(1, 0, 1) 
rdate = seq(as.Date("1951-01-15"),as.Date("2013-12-15"),by="month")
rnday = sapply(rdate, FUN = ndays)

rrca.flow = readBud(fname = paste0(rrca.path, "baseline\\2000.2.csv"), zn = rzone)
nper = 600
for (i in 2001:2013) {
  flow = readBud(fname = paste0(rrca.path, "baseline\\", i, ".2.csv"), zn = rzone)
  flow$PERIOD = flow$PERIOD + nper
  rrca.flow = rbind(rrca.flow, flow)
  nper = nper + 12
}
tail(rrca.flow$PERIOD, 100)
rrca.flow.mon = aggregate(rrca.flow[,5:11], list(MONTH=rrca.flow$PERIOD), mean)
rrca.flow.mon$MONTH = rdate
rrca.flow.mon[,2:8] = rrca.flow.mon[,2:8] * rnday * 86400 / rrca.area * 12 # inch / month
rrca.flow.year = aggregate(rrca.flow.mon[,2:8], list(YEAR=year(rdate)), sum) 

#### mcb budget ####

mcb.zon = scan(paste0(mcb.path, "MedicineCreek.zone"), skip=262)
mcb.zon.sum = summarise(group_by(data.frame(zone = mcb.zon), zone), n=n())
mcb.zon.n = sum(mcb.zon.sum$n[mcb.zon.sum$zone %in% c(1)])
mcb.area = mcb.zon.n * 1320 ^ 2

mzone = c(1, 0)
mdate = seq(as.Date("1950-01-15"),as.Date("2013-12-15"),by="month")
mnday = sapply(mdate, FUN = ndays)

mcb.flow = readBud(fname = paste0(mcb.path, "baseline\\mcb.2.csv"), zn = rzone)
mcb.flow.mon = aggregate(mcb.flow[,5:14], list(MONTH=mcb.flow$PERIOD), mean)
mcb.flow.mon$MONTH = mdate
mcb.flow.mon[,2:11] = mcb.flow.mon[,2:11] * mnday / mcb.area * 12 # inch / month  # /acre2foot # / mcb.area * 12 # inch / month
mcb.flow.year = aggregate(mcb.flow.mon[,2:11], list(YEAR=year(mdate)), sum) 


#### budget plots ####
pn = function(x) if (x >= 0) return("Flow In To Aquifer") else return("Flow Out From Aquifer")

modflow_type = function(x) {
  if (x == "RESERV. LEAKAGE") x = "BASEFLOW"
  if (x == "RIVER LEAKAGE") x = "BASEFLOW"
  if (x == "STREAM LEAKAGE") x = "BASEFLOW"
  if (x == "STREAM.LEAKAGE") x = "BASEFLOW"
  if (x == "MNW2") x = "PUMPING"
  if (x == "WELLS") x = "PUMPING"
  return(x)
}

mcb.flow.mon$STREAM.LEAKAGE = mcb.flow.mon$`STREAM LEAKAGE` + mcb.flow.mon$`RESERV. LEAKAGE` + mcb.flow.mon$`RIVER LEAKAGE`
mflow = melt(select(mcb.flow.mon, -DRAINS, -`HEAD DEP BOUNDS`, -`RESERV. LEAKAGE`, -`RIVER LEAKAGE`,-`STREAM LEAKAGE`), id.vars="MONTH")
mflow$InOut = sapply(mflow$value, pn)
mflow$type = sapply(as.character(mflow$variable) , modflow_type)
mflow$model = "MCB"

#ggplot(filter(rflow, type=="PUMPING"), aes(x = MONTH, y = value)) + geom_line(color=4)
rflow = melt(select(rrca.flow.mon, -DRAINS), id.vars="MONTH")
rflow$InOut = sapply(rflow$value, pn)
rflow$type = sapply(as.character(rflow$variable), modflow_type)
rflow$model = "RRCA"

g <- ggplot(filter(rbind(mflow, rflow), value != 0), aes(x = MONTH, y = value))  +
  geom_line(color=4) +
  xlab("Time") +
  ylab("Net Flow Rate (inch/month)") +
  facet_grid(type ~ model, scales="free_y")  + theme(legend.position="none")
scale=1.5
ggsave("monthly_budget.png",g,dpi=600,width=6*scale, height=7*scale,units="in")



### cumulative ###
mcb.flow.cum = mcb.flow.mon[-(1:12),]
rrca.flow.cum = rrca.flow.mon
for (i in 2:11) mcb.flow.cum[,i] = cumsum(mcb.flow.cum[,i])
for (i in 2:8) rrca.flow.cum[,i] = cumsum(rrca.flow.cum[,i])

mcb.flow.cum$STREAM.LEAKAGE = mcb.flow.cum$`STREAM LEAKAGE` + mcb.flow.cum$`RESERV. LEAKAGE` + mcb.flow.cum$`RIVER LEAKAGE`
mflow = melt(select(mcb.flow.cum, -DRAINS, -`HEAD DEP BOUNDS`, -`RESERV. LEAKAGE`, -`RIVER LEAKAGE`,-`STREAM LEAKAGE`), id.vars="MONTH")
mflow$InOut = sapply(mflow$value, pn)
mflow$type = sapply(as.character(mflow$variable) , modflow_type)
mflow$model = "MCB"


rflow = melt(select(rrca.flow.cum, -DRAINS), id.vars="MONTH")
rflow$InOut = sapply(rflow$value, pn)
rflow$type = sapply(as.character(rflow$variable), modflow_type)
rflow$model = "RRCA"

g <- ggplot(filter(rbind(mflow, rflow), value != 0), aes(x = MONTH, y = value, color = model))  +
  geom_line() +
  xlab("Time") +
  ylab("Cumulative Flow (inch)") +
  facet_wrap(~type, ncol=1, scales="free_y") 
scale=1.5
ggsave("cum_budget.png",g,dpi=600,width=4*scale, height=5*scale,units="in")


### annual flow ###

flow.annual = filter(rbind(mflow, rflow), MONTH==as.Date("2013-12-15"))
flow.annual = filter(flow.annual, value != 0)
flow.annual$value = flow.annual$value / (2013 - 1951 + 1)
flow.annual$value[flow.annual$type %in% c("BASEFLOW", "ET", "PUMPING")] = - flow.annual$value[flow.annual$type %in% c("BASEFLOW", "ET", "PUMPING")]
g = ggplot(flow.annual, aes(x=type,y=value,fill=model)) +
  geom_bar(stat="identity",position="dodge") +
  xlab("Flow Term") +
  ylab("Mean Annual Flow (inch/year)")
scale=1
ggsave("annual_budget.png",g,dpi=600,width=6*scale, height=4*scale,units="in")


# depletion test ----------------------------------------------------------

### mcb model ###


mcb.dpump = data.frame(MONTH = mcb.flow.mon$MONTH)

loc = "Loc2"


### rrca model ###
paste0(rrca.path, loc, "\\2000.2.csv")

read.mcb.bud = function(fname){
  mcb.dflow = readBud(fname = fname, zn = mzone)
  mcb.dflow.mon = aggregate(mcb.dflow[,-(1:4)], list(MONTH=mcb.dflow$PERIOD), mean)
  mcb.dflow.mon$MONTH = mdate
  
  mcb.dflow.mon$STREAM.LEAKAGE = mcb.dflow.mon$`STREAM LEAKAGE` + mcb.dflow.mon$`RESERV. LEAKAGE` + mcb.dflow.mon$`RIVER LEAKAGE`
  
  mflow = melt(select(mcb.dflow.mon, -`STREAM LEAKAGE`, -`RESERV. LEAKAGE`, -`RIVER LEAKAGE`), id.vars="MONTH", variable.name = "Flow Term")
  mflow$value = mflow$value / acre2foot 
  mflow$Model = "MCB"
  return(mflow)
}

read.rrca.bud = function(basedir){
  if (!(endsWith(basedir, "\\") | endsWith(basedir, "/"))) basedir = paste0(basedir,"/")
  rrca.dflow = readBud(fname = paste0(basedir, "2000.2.csv"), zn = rzone)
  nper = 600
  for (i in 2001:2013) {
    flow = readBud(fname = paste0(basedir, i, ".2.csv"), zn = rzone)
    flow$PERIOD = flow$PERIOD + nper
    rrca.dflow = rbind(rrca.dflow, flow)
    nper = nper + 12
  }
  
  rrca.dflow.mon = aggregate(rrca.dflow[,-(1:4)], list(MONTH=rrca.dflow$PERIOD), mean)
  rrca.dflow.mon$MONTH = rdate
  #rrca.dflow.mon[,-1] = rrca.dflow.mon[,-1] - rrca.flow.mon[,-1]
  
  rflow = melt(rrca.dflow.mon, id.vars="MONTH", variable.name = "Flow Term")
  rflow$value = rflow$value * 86400 / acre2foot # acre-foot/day
  rflow$Model = "RRCA"
  return(rflow)
}

# read baseline
rflow.bas = read.rrca.bud(basedir = paste0(rrca.path, "baseline"))
mflow.bas = read.mcb.bud(paste0(mcb.path, "baseline.2.csv"))

rflow.1 = read.rrca.bud(paste0(rrca.path, loc))
rflow.1$value = rflow.1$value - rflow.bas$value
rflow.1$InOut = sapply(rflow.1$value, pn)

mflow.1 = read.mcb.bud(paste0(mcb.path, loc,".2.csv"))
mflow.1$value = mflow.1$value - mflow.bas$value
mflow.1$InOut = sapply(mflow.1$value, pn)

g = ggplot(data = filter(rbind(rflow.1, mflow.1), value != 0), aes(x = MONTH, y = value, fill=`Flow Term`, group=`Flow Term`))  +
  geom_area() +
  xlab("Time") +
  ylab("Change in Flow Rate (acre-foot/day)") +
  scale_fill_brewer(palette = "Set1") +
  facet_grid(InOut ~ Model, scales="free_y")
scale=1
ggsave(paste0("depletion_mcb_", loc, ".png"),g,dpi=600,width=7*scale, height=5*scale,units="in")


ggplot(filter(rflow.1, `Flow Term`=="STREAM LEAKAGE"), aes(x=MONTH,y=value)) + geom_line()

g = ggplot(data = rflow, aes(x = MONTH, y = value, fill=variable, group=variable))  +
  geom_area() +
  xlab("Time") +
  ylab("Flow rate (inch/month)") +
  scale_fill_brewer(palette = "Set1") +
  facet_grid(FLOW ~ ., scales="free_y")


rrca.dflow[,-(1:4)] = rrca.dflow1[,-(1:4)] - rrca.dflow[,-(1:4)]
plot(rrca.dflow$WELLS)
rrca.dflow$PERIOD = 1:nrow(rrca.dflow)
rflow = melt(select(rrca.dflow, -TIME,-ZONE,-STEP), id.vars="PERIOD", variable.name = "Flow Term")
rflow$InOut = sapply(rflow$value, pn)
ggplot(data = rflow, aes(x = PERIOD, y = value, fill=`Flow Term`, group=`Flow Term`))  +
  geom_area() +
  scale_fill_brewer(palette = "Set1") +
  facet_grid(InOut ~ ., scales="free_y")
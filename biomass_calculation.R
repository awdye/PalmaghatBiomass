library(dplR)
library(sp)
library(foreign)
library(lattice)
library(ggplot2)
library(sciplot)


setwd("/Users/dmb/Documents/DMB/NE_USA/Palmaghatt/Biomass/R")
pm.field.data<-read.table("PMall_FieldData.csv", header=T,sep=",")
pm.field.data$dbhclass<-cut(pm.field.data$dbh,breaks=c(5,8.999,12.499,24.999,39.999,1000), right=T,
                       labels=c(5,9,12.5,25,40),include.lowest=T)
pm.field.data$convha<-cut(pm.field.data$dbh,breaks=c(5,8.999,12.499,24.999,39.999,1000), right=T,
                     labels=10000/(pi*(c(3,5,8,20,30)^2)),include.lowest=T)
pm.field.data$convha<-as.numeric(as.character(pm.field.data$convha))

bm.eqs<-read.table("biomass_coeff.csv", header=T,sep=",")

dim(unique(subset(pm.field.data, select=tree, site=="PM3")))
unique(subset(ab.pm.all, select=tree, site=="PM3"))
subset(pm.field.data, site=="PM1")


# Function for current biomass
# al.eq is a data.frame with coefficients (a & b) for alometric equations of 
# the form a*DBH^b
# dbh.data is a data.frame with dbh data

setwd("/Users/dmb/Documents/DMB/NE_USA/Palmaghatt/Biomass/RawRW_Data_by_sp")
files_to_read<-list.files()
# files_to_read<-"PM_TulipPoplar.rw"

core.means.all<-NULL
for (i in files_to_read){
  file_read<-read.rwl(i)
  assign(i,file_read)

cores<-data.frame("code"=names(file_read))
cores$recode<- substr(as.character(cores$code),1,5)
  core.means<-data.frame("year"=row.names(file_read))
  core.means<-NULL
  n=1
  file_read$year<-row.names(file_read)
for (j in unique(cores$recode)){
  print(j)
  print(n)
  tree1<-data.frame(file_read[,c(as.character(cores$code[cores$recode==j]),"year")])
  if(dim(tree1)[2]>1){
    tree2<-data.frame(rowMeans(subset(tree1, select=-c(year)),na.rm=T))  
    } else {tree2=tree1}
  tree2<-na.omit(tree2)
    tree2$tree<-j
    tree2$year<-as.numeric(as.character(row.names(tree2)))
      tree2$rACUM<-0
      tree2$rACUM.inv<-0
        for(k in 1:dim(tree2)[1]){
          tree2$rACUM[k]<-sum(tree2[c(1:k),1])
          tree2$rACUM.inv[k-1]<-sum(tree2[c((k):(dim(tree2)[1])),1])
        }
    names(tree2)<-c("meanRW","recode","year","rACUM","rACUM.inv")
    if(n==1)(core.means<-tree2)else(core.means<-rbind(core.means, tree2))
     n=n+1}
 plot(core.means$rACUM.inv~core.means$year, type="p")
  title(i)
  
  core.means.all<-rbind(core.means.all,core.means)
}
core.means.all$site<-substr(as.character(core.means.all$recode),1,3)
dim(unique(subset(core.means.all, select=tree, site=="PM5")))


core.means.all$tree<-as.numeric(as.character(substr(as.character(core.means.all$recode),4,5)))

ab.pm.all<-merge(core.means.all,pm.field.data, by=c("site","tree"))
bm.eqs.sp<-bm.eqs[bm.eqs$eq==1,]
ab.pm.all<-merge(ab.pm.all,bm.eqs.sp,by="species")

#dbhest1==> dbh at the BEGINNING of the growing season
ab.pm.all$dbhest1<-ab.pm.all$dbh-0.2*(ab.pm.all$rACUM.inv+ab.pm.all$meanRW)

#dbhest2==> dbh at the END of the growing season
ab.pm.all$dbhest2<-ab.pm.all$dbh-0.2*(ab.pm.all$rACUM.inv)

ab.pm.all$AB<-ab.pm.all$a*ab.pm.all$dbhest2^ab.pm.all$b

#annAB==>biomass acculumated during the year.
# that is, biomass at the END of the growing season
# minus the biomass at the BEGINNING of the growing season
ab.pm.all$annAB<-ab.pm.all$a*(ab.pm.all$dbhest2^ab.pm.all$b-ab.pm.all$dbhest1^ab.pm.all$b)

#biomass per hectare assuming the plot
#represents the entire hectare
ab.pm.all$annAB.ha<-ab.pm.all$annAB*ab.pm.all$convha


#Calculates total biomass per ha and number of trees
#first we need to define two basic functions
sum.fn<-function(x) sum(x, na.rm=TRUE)
count.fn<-function(x) length(unique(x, na.rm=TRUE))

ab.ha.site<-data.frame(tapply(X=ab.pm.all$annAB,INDEX=list(ab.pm.all$year,ab.pm.all$site),sum.fn))
ab.ha.site$year<-as.numeric(as.character(row.names(ab.ha.site)))
ab.ha.site[is.na(ab.ha.site)]=0

n.trees<-data.frame(tapply(X=ab.pm.all$tree,INDEX=list(ab.pm.all$year,ab.pm.all$site),count.fn))
n.trees$year<-as.numeric(as.character(row.names(n.trees)))
plot(n.trees[,5])
n.trees[is.na(n.trees)]=0



#Calculations END HERE
#From here on it's JUST GRAPHS

p<-ggplot(data=ab.pm.all,aes(x=year, y=AB,color=species))
p+geom_point(size=1)+facet_wrap(~site)

p<-ggplot(data=ab.pm.all,aes(x=year, y=annAB,color=species))
p+geom_point(size=1)+facet_wrap(~site)

p<-ggplot(data=ab.pm.all,aes(x=year, y=annAB.ha,color=species))
p+geom_point(size=1)+facet_wrap(~site)

lineplot.CI(x.factor=year,response=annAB.ha,group=site,
              data=subset(ab.pm.all,site!="PM4"&year>1960))

lineplot.CI(x.factor=year,response=meanRW,group=species,
            data=subset(ab.pm.all,site!="PM4"&site!="PM1"&year>1960))

lineplot.CI(x.factor=year,response=annAB,group=site,
            data=subset(ab.pm.all,site!="PM4"&year>1960))


tot.ab.ha<-lineplot.CI(x.factor=year,response=annAB,group=site,
            data=subset(ab.pm.all,site!="PM4"&year>1960),fun=sum.fn,
             ylab="",xlab="")
axis(1,at=seq(1,60,5),labels=F)
mtext(side=1,"Year",line=2)
mtext(side=2,"Annual above-ground biomass (kg/ha)",line=2)


ggplot(data=subset(ab.pm.all,site!="PM4"), 
       aes(y = annAB, x = year,color=site)) + 
  stat_summary(fun.y = 'sum', fun.ymin = function(x) 0, geom = 'point', 
               aes(fill =site), position = 'dodge',
               fun.ymin = min, fun.ymax = max) 








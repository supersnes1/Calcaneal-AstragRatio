getOneSpecimenAverage <- function(this.specimen, calcmeasure) {
	this.measurements <- c("Length.Calc", "Length..medial.", "Length..lateral.", "Width..between.prox.condyles.", "Width..between.distal.condyles.navicular.facet.", "Length.to.Sustentacular.Process", "Length.to.Calcaneoastragular.facet", "Calcaneoastragular.facet.to.cuboid.facet")
	colMeans(calcmeasure[calcmeasure$Catalog.Number== this.specimen, this.measurements], na.rm=TRUE)
}

getSpecimenAverages <- function(calcmeasure) {
	specimen.vec <- sort(unique(calcmeasure$Catalog.Number))
	n <- t(sapply(specimen.vec, function(x) calcmeasure[calcmeasure$Catalog.Number==x, c("Museum.Acronym", "Catalog.Number", "Genus", "Species")][1,]))
	m <- data.frame(n, t(sapply(specimen.vec, getOneSpecimenAverage, calcmeasure)))
	m
}

plotOneProjectedSpecimen <- function(this.coord, this.taxon) {
	arrows(x0=this.coord[1], x1=pca$x[this.taxon, this.x], y0=this.coord[2], y1=pca$x[this.taxon, this.y], length=0.10, angle=60, col="gray75")
}

plotMultipleProjectedSpecimens <- function(this.taxon) {
	apply(m.proj[m.short$binom==this.taxon,c(this.x, this.y)], 1, plotOneProjectedSpecimen, this.taxon)
}

#plotOneProjectedMeasures<- function(this.coord, this.taxon) {
# arrows(x0=this.coord[1], x1=pca$x[this.taxon, this.x], y0=this.coord[2], y1=pca$x[this.taxon, this.y], length=0.10, angle=60, col="gray75")
#}

#plotMultipleProjectedMeasures <- function(this.taxon) {
#  apply(m.proj[m.short$binom==this.taxon,c(this.x, this.y)], 1, plotOneProjectedMeasures, this.taxon)
#}

dat <- read.csv("~/Dropbox/CalcanealGearRatio/LACM_CalcaneusMeasurements_2018_8_29.csv", stringsAsFactors=FALSE)
dat <- dat[dat$Genus !="Dicerorhinus" & dat$Genus !="Diceros" & dat$Genus !="Hippopotamus",]

#need to parse input data so only photo measures are being utilized
dat <- dat[dat$Measurement.Type.Manual.Photo. == "Photo",]

m <- getSpecimenAverages(dat)
#remove total length of calcaneus to avoid redundancy
m <- m[,!colnames(m) == "Length.Calc"]

m$binom <- apply(m, 1, function(x) paste(x["Genus"], x["Species"], sep="_"))

m.short <- m[,c("binom", names(m)[sapply(m, class)=="numeric"])]
rownames(m.short) <- m$Catalog.Number 
m.short <- m.short[complete.cases(m.short),]
m.short <- cbind(binom=m.short$binom, log(m.short[,names(m)[sapply(m, class)=="numeric"]]))

m.sp <- aggregate(m[,sapply(m, class)=="numeric"], by=list(m$binom), FUN=mean, na.rm=TRUE)
rownames(m.sp) <- m.sp[,1]
m.sp <- m.sp[,-1]
m.sp.short <- log(m.sp[complete.cases(m.sp),])

pca <- prcomp(m.sp.short)

tax <- table(m$binom)
tax <- names(tax)[tax>1]

m.proj <- data.matrix(m.short)[,sapply(m.short, class)=="numeric"]
m.proj <- t(apply(m.proj, 1,  function(x) x - pca$center))
m.proj <- m.proj %*% pca$rotation

#make single plot to track down where specific specimens occur in relation to species average
quartz(height = 8, width = 11)
this.x <- 2
this.y <- 3
plot(pca$x[,this.x], pca$x[,this.y], type="n", xlim = range(m.proj[,this.x]), ylim = range(m.proj[,this.y]), xlab=paste("PC", this.x), ylab=paste("PC", this.y))
abline(h=0, lty=3, col="gray50")
abline(v=0, lty=3, col="gray50")

text(pca$x[,this.x], pca$x[,this.y], labels=rownames(pca$x), cex=0.5)

sapply(tax, plotMultipleProjectedSpecimens)
			 
text(m.proj[m.short$binom %in% tax,this.x], m.proj[m.short$binom %in% tax,this.y], labels=m.short$binom[m.short$binom %in% tax], col="dodgerblue", cex=0.5)
text(m.proj[m.short$binom %in% tax,this.x], m.proj[m.short$binom %in% tax,this.y], labels=rownames(m.short)[m.short$binom %in% tax], col="dodgerblue", cex=0.5, pos=1)

#select for problematic species (e.g. Cervus canadensis)
tax.cat <- c("Alces_alces", "Cervus_canadensis", "Kobus_kob", "Rangifer_tarandus")
text(m.proj[m.short$binom %in% tax.cat,this.x], m.proj[m.short$binom %in% tax.cat,this.y], labels=m.short$binom[m.short$binom %in% tax.cat], col="red", cex=0.5)
text(m.proj[m.short$binom %in% tax.cat,this.x], m.proj[m.short$binom %in% tax.cat,this.y], labels=rownames(m.short)[m.short$binom %in% tax.cat], col="red", cex=0.5, pos=1)

#make group of plots that step through all species with multiuple specimens
quartz(height = 48, width = 11)
par(mfrow=c((length(tax)/2)+1,2))

for(xx in seq(0, length(tax),1))
{
  this.x <- 2
  this.y <- 3
  plot(pca$x[,this.x], pca$x[,this.y], type="n", xlim = range(m.proj[,this.x]), ylim = range(m.proj[,this.y]), xlab=paste("PC", this.x), ylab=paste("PC", this.y))
  abline(h=0, lty=3, col="gray50")
  abline(v=0, lty=3, col="gray50")
  
  text(pca$x[,this.x], pca$x[,this.y], labels=rownames(pca$x), cex=0.5)
  
  sapply(tax, plotMultipleProjectedSpecimens)
  
  text(m.proj[m.short$binom %in% tax,this.x], m.proj[m.short$binom %in% tax,this.y], labels=m.short$binom[m.short$binom %in% tax], col="dodgerblue", cex=0.5)
  text(m.proj[m.short$binom %in% tax,this.x], m.proj[m.short$binom %in% tax,this.y], labels=rownames(m.short)[m.short$binom %in% tax], col="dodgerblue", cex=0.5, pos=1)
  
  if(xx > 0)
  {
    #select for problematic species (e.g. Cervus canadensis)
    text(m.proj[m.short$binom %in% tax[xx],this.x], m.proj[m.short$binom %in% tax[xx],this.y], labels=m.short$binom[m.short$binom %in% tax[xx]], col="red", cex=0.5)
    text(m.proj[m.short$binom %in% tax[xx],this.x], m.proj[m.short$binom %in% tax[xx],this.y], labels=rownames(m.short)[m.short$binom %in% tax[xx]], col="red", cex=0.5, pos=1)
  }
}

#######################################################################################################################################################
######project specimen data into plot to verify that measures were not overly variable  (look to plot manual vs photo measures)
names <-  c("Museum.Acronym", "Catalog.Number", "Genus", "Species", "Element", "Length..medial.", "Length..lateral.", 
            "Width..between.prox.condyles.", "Width..between.distal.condyles.navicular.facet.", 
            "Length.to.Sustentacular.Process", "Length.to.Calcaneoastragular.facet", "Calcaneoastragular.facet.to.cuboid.facet")
names.astrag <- c("Length..medial.", "Length..lateral.", "Width..between.prox.condyles.", "Width..between.distal.condyles.navicular.facet.")
dat.specFixed <- data.frame()
dat.cropped <- dat[,colnames(dat) %in% names]
##remove specimens that do not match those already in the plot
dat.cropped <- dat.cropped[dat.cropped$Catalog.Number %in% rownames(m.short),]
#need to collapse measures down into separate lines so that calcaneus and astragalus are together (ignore manual measures since they are incomplete in some cases)
for(xx in unique(dat.cropped$Catalog.Number))
{
  dat.specimen <- dat.cropped[dat.cropped$Catalog.Number == xx,]
  dat.calc <- dat.specimen[dat.specimen$Element == "calcaneus",]
  dat.astrag<- dat.specimen[dat.specimen$Element == "astragalus",]
  
  #need to check if rows for calcaneus and astragalus are the same number 
  if(!nrow(dat.calc) ==  nrow(dat.astrag))
    {
    print(paste(xx,"has unequal number of calcaneus and astragalar measurements"), sep=" ")
    #should only be 88929 which is astragalus only measures so it can be removed
  }
  dat.calc[,names.astrag] <- dat.astrag[,names.astrag]
  dat.specFixed <- rbind(dat.specFixed, dat.calc)
}

dat.specFixed <- dat.specFixed[complete.cases(dat.specFixed),]


dat.proj <- log(dat.specFixed[,sapply(dat.specFixed, class)=="numeric"])
dat.proj <- t(apply(dat.proj, 1,  function(x) x - pca$center))
dat.proj <- dat.proj %*% pca$rotation  #getting values that are exaggerated
#predict(object=pca, newdata = dat.proj)

dat.proj <- cbind(dat.proj, dat.specFixed$Catalog.Number)
dat.proj <- as.data.frame(dat.proj)
rownames(dat.proj) <- NULL

#n.rname.rows <- table(dat.proj$V1)
#rnames <- vector()
#rnames.list <- vector()
#for(yy in unique(dat.proj$V1))
#{
#  for(kk in seq(1,n.rname.rows[names(n.rname.rows) == yy],1))
#  {
#    rnames[kk] <- paste(yy, paste("_",kk,sep=""),sep="") 
#  }
#  rnames.list <- append(rnames.list, rnames)
#}
#dat.proj[dat.proj$V1 == yy,]$V1 <- rnames
#dat.proj[,!colnames(dat.proj)=="V1"] <- as.numeric(dat.proj[,!colnames(dat.proj)=="V1"])

quartz(height = 8, width = 11)
this.x <- 2
this.y <- 3
plot(pca$x[,this.x], pca$x[,this.y], type="n", xlim = range(m.proj[,this.x]), ylim = range(m.proj[,this.y]), xlab=paste("PC", this.x), ylab=paste("PC", this.y))
abline(h=0, lty=3, col="gray50")
abline(v=0, lty=3, col="gray50")

text(pca$x[,this.x], pca$x[,this.y], labels=rownames(pca$x), cex=0.5)

#sapply(tax, plotMultipleProjectedSpecimens)

text(m.proj[m.short$binom %in% tax,this.x], m.proj[m.short$binom %in% tax,this.y], labels=m.short$binom[m.short$binom %in% tax], col="dodgerblue", cex=0.5)
text(m.proj[m.short$binom %in% tax,this.x], m.proj[m.short$binom %in% tax,this.y], labels=rownames(m.short)[m.short$binom %in% tax], col="dodgerblue", cex=0.5, pos=1)

text(m.proj[!m.short$binom %in% tax,this.x], m.proj[!m.short$binom %in% tax,this.y], labels=rownames(m.short)[!m.short$binom %in% tax], col="dodgerblue", cex=0.5, pos=1)
#plot unique specimens into PCA space
dat.x <- dat.proj[,this.x]
dat.y <- dat.proj[,this.y]
names(dat.x) <- NULL
names(dat.y) <- NULL

text(as.numeric(as.character(dat.x)), as.numeric(as.character(dat.y)), labels=dat.specFixed$Catalog.Number, col="chartreuse4", cex=0.5)

#make group of plots that step through all species with multiuple specimens
quartz(height = 48, width = 11)
par(mfrow=c((length(tax)/2)+1,2))

for(xx in seq(0, length(tax),1))
{
  this.x <- 2
  this.y <- 3
  plot(pca$x[,this.x], pca$x[,this.y], type="n", xlim = range(m.proj[,this.x]), ylim = range(m.proj[,this.y]), xlab=paste("PC", this.x), ylab=paste("PC", this.y))
  abline(h=0, lty=3, col="gray50")
  abline(v=0, lty=3, col="gray50")
  
  text(pca$x[,this.x], pca$x[,this.y], labels=rownames(pca$x), cex=0.5)
  
  sapply(tax, plotMultipleProjectedSpecimens)
  
  text(m.proj[m.short$binom %in% tax,this.x], m.proj[m.short$binom %in% tax,this.y], labels=m.short$binom[m.short$binom %in% tax], col="dodgerblue", cex=0.5)
  text(m.proj[m.short$binom %in% tax,this.x], m.proj[m.short$binom %in% tax,this.y], labels=rownames(m.short)[m.short$binom %in% tax], col="dodgerblue", cex=0.5, pos=1)
  
  if(xx > 0)
  {
    #select for problematic species (e.g. Cervus canadensis)
    text(m.proj[m.short$binom %in% tax[xx],this.x], m.proj[m.short$binom %in% tax[xx],this.y], labels=m.short$binom[m.short$binom %in% tax[xx]], col="red", cex=0.5)
    text(m.proj[m.short$binom %in% tax[xx],this.x], m.proj[m.short$binom %in% tax[xx],this.y], labels=rownames(m.short)[m.short$binom %in% tax[xx]], col="red", cex=0.5, pos=1)
  
    text(dat.proj[,this.x], dat.proj[,this.y], labels=dat.specFixed$Catalog.Number, col="chartreuse4", cex=0.5)
  }
}


######################################################################################################################################
##Check astragalus measures with those of Barr 2014
barr.dat <- read.csv("/Users/emdoughty/Dropbox/CalcanealGearRatio/Data_Barr2014.csv", stringsAsFactors=FALSE)
barr.dat$Taxon <- gsub(" ", "_", barr.dat$Taxon)

#Check if unwanted genera are removed
dat.tax <- matrix(nrow = nrow(barr.dat), ncol = 2)
colnames(dat.tax) <- c("Genus", "Species") 
test <- strsplit(barr.dat$Taxon, "_")
for(xx in seq(1, length(test),1)) { dat.tax[xx,] <- c(test[[xx]][1], test[[xx]][2])}
dat <- cbind(barr.dat, dat.tax)
nrow(barr.dat)
barr.dat <- barr.dat[dat$Genus !="Dicerorhinus" & dat$Genus !="Diceros" & dat$Genus !="Hippopotamus",]
nrow(barr.dat)
barr.dat <- barr.dat[,!colnames(barr.dat) %in% c("Genus", "Species")]
barr.dat <- barr.dat[,c(1,2,4,3,5,6)]

#rerun pca for astragalus only
astrag.sp.short <- log(m.sp[complete.cases(m.sp), c("Length..medial.", "Length..lateral.", "Width..between.prox.condyles.",
                                                    "Width..between.distal.condyles.navicular.facet.")])
astrag.pca <- prcomp(astrag.sp.short)

astrag.short <- m.short[, colnames(m.short) %in% c( "binom","Length..medial.", "Length..lateral.", 
                                      "Width..between.prox.condyles.","Width..between.distal.condyles.navicular.facet.")]

astrag.proj <- astrag.short[,sapply(astrag.short, class)=="numeric"]
astrag.proj <- t(apply(astrag.proj, 1,  function(x) x - astrag.pca$center))
astrag.proj <- astrag.proj %*% astrag.pca$rotation

#prepare for comparison with my data
barr.short <- barr.dat$Taxon
barr.specNum <- barr.dat$individual
barr.dat <- barr.dat[,colnames(barr.dat) %in% c("individual","Taxon", "LML", "MML","WAF", "WAT")]

barr.proj <- log(barr.dat[,sapply(barr.dat, class)=="numeric"])
barr.proj <- t(apply(barr.proj, 1,  function(x) x - astrag.pca$center))
barr.proj <- barr.proj %*% astrag.pca$rotation  #getting values that are exaggerated

barr.tax <- barr.short[barr.short %in% rownames(astrag.sp.short)]#get list of taxa shared with Barr 2014
barr.tax <- which(barr.short %in% rownames(astrag.sp.short))#get list of taxa shared with Barr 2014


#plot comaprison
quartz(height = 8, width = 11)
this.x <- 2
this.y <- 3
plot(astrag.pca$x[,this.x], astrag.pca$x[,this.y], type="n",xlim = range(c(astrag.proj[,this.x], astrag.pca$x[,this.x])), 
     ylim = range(c(astrag.proj[,this.y], astrag.pca$x[,this.y])), 
     xlab=paste("PC", this.x), ylab=paste("PC", this.y))
abline(h=0, lty=3, col="gray50")
abline(v=0, lty=3, col="gray50")

text(astrag.pca$x[,this.x], astrag.pca$x[,this.y], labels=rownames(astrag.pca$x), cex=0.5)

#sapply(tax, plotMultipleProjectedSpecimens)

text(astrag.proj[astrag.short$binom %in% tax,this.x], astrag.proj[astrag.short$binom %in% tax,this.y], labels=astrag.short$binom[astrag.short$binom %in% tax], col="dodgerblue", cex=0.5)
text(astrag.proj[astrag.short$binom %in% tax,this.x], astrag.proj[astrag.short$binom %in% tax,this.y], labels=rownames(astrag.short)[astrag.short$binom %in% tax], col="dodgerblue", cex=0.5, pos=1)

text(barr.proj[barr.tax,this.x], barr.proj[barr.tax,this.y], labels=barr.short[barr.tax], col="orange", cex=0.5)
text(barr.proj[barr.tax,this.x], barr.proj[barr.tax,this.y], labels=barr.specNum[barr.tax], col="orange", cex=0.5, pos=1)

#make group of plots that step through all species with multiuple specimens
quartz(height = 48, width = 11)
barr.uni.sp <- unique(barr.short[barr.tax])
par(mfrow=c((length(barr.uni.sp)/2)+1,2))

for(xx in seq(0, length(barr.uni.sp),1))
{
  this.x <- 2
  this.y <- 3
  plot(astrag.pca$x[,this.x], astrag.pca$x[,this.y], type="n",xlim = range(c(astrag.proj[,this.x], astrag.pca$x[,this.x], barr.proj[,this.x])), 
       ylim = range(c(astrag.proj[,this.y], astrag.pca$x[,this.y], barr.proj[,this.y])), 
       xlab=paste("PC", this.x), ylab=paste("PC", this.y))
  abline(h=0, lty=3, col="gray50")
  abline(v=0, lty=3, col="gray50")
  
  text(astrag.pca$x[,this.x], astrag.pca$x[,this.y], labels=rownames(astrag.pca$x), cex=0.5)
  
  #sapply(tax, plotMultipleProjectedSpecimens)

  text(astrag.proj[astrag.short$binom %in% tax,this.x], astrag.proj[astrag.short$binom %in% tax,this.y], labels=astrag.short$binom[astrag.short$binom %in% tax], col="dodgerblue", cex=0.5)
  text(astrag.proj[astrag.short$binom %in% tax,this.x], astrag.proj[astrag.short$binom %in% tax,this.y], labels=rownames(astrag.short)[astrag.short$binom %in% tax], col="dodgerblue", cex=0.5, pos=1)
  
  text(barr.proj[barr.tax,this.x], barr.proj[barr.tax,this.y], labels=barr.short[barr.tax], col="darkseagreen", cex=0.5)
  text(barr.proj[barr.tax,this.x], barr.proj[barr.tax,this.y], labels=barr.specNum[barr.tax], col="darkseagreen", cex=0.5, pos=1)
  
  if(xx > 0)
  {
    #select for problematic species (e.g. Cervus canadensis)
    text(astrag.pca$x[rownames(astrag.sp.short) %in% barr.uni.sp[xx],this.x], astrag.pca$x[rownames(astrag.sp.short) %in% barr.uni.sp[xx],this.y], 
         labels=rownames(astrag.pca$x)[xx],
         cex=0.5, col="darkblue")
   
    text(astrag.proj[astrag.short$binom %in% barr.uni.sp[xx],this.x], astrag.proj[astrag.short$binom %in% barr.uni.sp[xx],this.y],
         labels=astrag.short$binom[astrag.short$binom %in% barr.uni.sp[xx]], col="red", cex=0.5)
    text(astrag.proj[astrag.short$binom %in% barr.uni.sp[xx],this.x], astrag.proj[astrag.short$binom %in% barr.uni.sp[xx],this.y], 
         labels=rownames(astrag.short)[astrag.short$binom %in% barr.uni.sp[xx]], col="red", cex=0.5, pos=1)
    
    text(barr.proj[barr.short %in% barr.uni.sp[xx],this.x], barr.proj[barr.short %in% barr.uni.sp[xx],this.y], 
         labels=barr.short[barr.short %in% barr.uni.sp[xx]], col="deeppink", cex=0.5)
    text(barr.proj[barr.short %in% barr.uni.sp[xx],this.x], barr.proj[barr.short %in% barr.uni.sp[xx],this.y], 
         labels=barr.specNum[barr.short %in% barr.uni.sp[xx]], col="deeppink", cex=0.5, pos=1)
  }
}


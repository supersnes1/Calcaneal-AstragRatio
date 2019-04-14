#############################################################################################
###Data Set-up
#############################################################################################
setwd("~/Dropbox/CalcanealGearRatio")
calcmeasure <- read.csv("~/Dropbox/CalcanealGearRatio/LACM_CalcaneusMeasurements_2018_8_29.csv", stringsAsFactors=FALSE)
calcmeasure <- calcmeasure[calcmeasure["Measurement.Type.Manual.Photo."] == "Photo",1:14]
calcSpecies <- unique(paste(calcmeasure$Genus, calcmeasure$Species,sep="_"))

MOM.Mat <- read.csv("~/Dropbox/CalcanealGearRatio/MOM_v4.1.csv")
#remove empty rows
MOM.Mat <- MOM.Mat[!MOM.Mat$Species == "",]
#aggregate to get single species entries
MOM.Mat[,"TaxonomicName"] <- paste(MOM.Mat[,"Genus"],MOM.Mat[,"Species"], sep="_")
#MOM.Mat[,"TaxonomicName"] <- gsub(pattern = " ", replacement =  "_", x = MOM.Mat[,"TaxonomicName"],)

MOM.Mat[,"CombinedMass(kg)"] <- MOM.Mat[,"Combined.Mass..g."]/1000
MOM.Mat[,"logMass.Kg"] <- log10(MOM.Mat[,"CombinedMass(kg)"])
MOM.Mat <- MOM.Mat[,c("Order","FAMILY","Genus","Species","TaxonomicName","CombinedMass(kg)","logMass.Kg")]
#combine same species if multiple entires(e.g. Muntiacus_reevesi has separate entries for mainland Asia and Insular)
MOM.Mat <- aggregate(MOM.Mat[,-c(1:5)], by = list(MOM.Mat$Order, MOM.Mat$FAMILY,MOM.Mat$Genus,
                                                  MOM.Mat$Species, MOM.Mat$TaxonomicName), mean)
colnames(MOM.Mat)[1:5] <- c("Order","FAMILY","Genus","Species","TaxonomicName")

ecoMat <- read.csv("~/Dropbox/CalcanealGearRatio/ExtantUngulate_BM_Eco_Mendoza_etal2007.csv")
ecoMat["TaxonomicName"] <- paste(ecoMat[,"Genus"],paste("_",ecoMat[,"Species"],sep=""),sep="")

#############################################################################################
###Analysis
#############################################################################################
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

getSpeciesAverages <- function(measureMat) {
  rownames(measureMat) <- NULL
  species.vec <- sort(unique(measureMat$taxon))
  this.measurements <- c("Length.Calc", "Length..medial.", "Length..lateral.", "Width..between.prox.condyles.", "Width..between.distal.condyles.navicular.facet.", "Length.to.Sustentacular.Process", "Length.to.Calcaneoastragular.facet", "Calcaneoastragular.facet.to.cuboid.facet")
  tst <- species.vec[1]
  this.species <- tst
  spMat <-matrix(nrow=length(species.vec),ncol = length(this.measurements))
  for(jj in seq(1, length(species.vec),1))
  {
    this.species <- species.vec[jj]
    spMat[jj,] <- colMeans(measureMat[measureMat$taxon == this.species, this.measurements], na.rm=TRUE)
  }
  colnames(spMat) <- this.measurements
  rownames(spMat) <- species.vec
  spMat
}

m <- getSpecimenAverages(calcmeasure)
m$taxon <- paste(m$Genus, tolower(m$Species), sep="_")
#m$R2 <- m$Length.to.Calcaneoastragular.facet/m$Calcaneoastragular.facet.to.cuboid.facet #this is based on the figure caption
m$R2 <- m$Calcaneoastragular.facet.to.cuboid.facet/m$Length.to.Calcaneoastragular.facet
#m$R2  <- m$Length.to.Calcaneoastragular.facet/m$Length.Calc #this is based on the designation listed in the figure
m$R3 <- m$Length.Calc/m$Length.to.Sustentacular.Process
m$A1 <- rowMeans(m[,c("Length..medial.", "Length..lateral.")], na.rm=TRUE)/rowMeans(m[,c("Width..between.prox.condyles.", "Width..between.distal.condyles.navicular.facet.")], na.rm=TRUE)
m$habitat <- ecoMat$Habitat[match(m$taxon, ecoMat$TaxonomicName)]
m$BM <- MOM.Mat$logMass.Kg[match(m$taxon, MOM.Mat$TaxonomicName)]

#family
m$Family <- MOM.Mat$FAMILY[match(m$taxon, MOM.Mat$TaxonomicName)]
#order
#family
m$Order <- MOM.Mat$Order[match(m$taxon, MOM.Mat$TaxonomicName)]

##remove Perissodactyls (3 species)
m <- m[!m$Order == "Perissodactyla",]

##species level
m.sp <- as.data.frame(getSpeciesAverages(m))
m.sp$taxon <- row.names(m.sp)
#m$R2 <- m$Length.to.Calcaneoastragular.facet/m$Calcaneoastragular.facet.to.cuboid.facet #this is based on the figure caption
m.sp$R2 <- m.sp$Calcaneoastragular.facet.to.cuboid.facet/m.sp$Length.to.Calcaneoastragular.facet
#m$R2  <- m$Length.to.Calcaneoastragular.facet/m$Length.Calc #this is based on the designation listed in the figure
m.sp$R3 <- m.sp$Length.Calc/m.sp$Length.to.Sustentacular.Process
m.sp$A1 <- rowMeans(m.sp[,c("Length..medial.", "Length..lateral.")], na.rm=TRUE)/rowMeans(m.sp[,c("Width..between.prox.condyles.", "Width..between.distal.condyles.navicular.facet.")], na.rm=TRUE)
m.sp$habitat <- ecoMat$Habitat[match(m.sp$taxon, ecoMat$TaxonomicName)]
m.sp$BM <- MOM.Mat$logMass.Kg[match(m.sp$taxon, MOM.Mat$TaxonomicName)]
m.sp$Family <- MOM.Mat$FAMILY[match(m.sp$taxon, MOM.Mat$TaxonomicName)]
#reorder
m.sp <- m.sp[,c(9,13,1:8,14, 10:12)]
m.sp <- m.sp[complete.cases(m.sp),] #68 species
m <- m.sp

habitat.colors <- c("dodgerblue", "chartreuse", "goldenrod1")
#### Polly 2010 Figure 4c
	# habitat.colors <- rainbow(length(unique(m$habitat)))
	plot(m$R2, m$R3, xlim=c(0,1), ylim=c(1, 1.5), pch=21, cex=1, bg=habitat.colors[m$habitat])
	# plot(m$R2, m$R3, pch=21, cex=1, bg=habitat.colors[m$habitat])
	text(m$R2, m$R3, labels=m$taxon, pos=3, cex=0.3)
	abline(lm(m$R3 ~ m$R2))

### Polly 2010 Figure 6
	plot(sort(m$R3), pch=21, bg=habitat.colors[m$habitat], type="n", ylim=c(1,1.5))
	rect(xleft=-10, ybottom=mean(m$R3, na.rm=TRUE)-sd(m$R3, na.rm=TRUE), xright=150, ytop=mean(m$R3, na.rm=TRUE)+sd(m$R3, na.rm=TRUE), col=adjustcolor("black", alpha.f=0.25))
	abline(h=mean(m$R3, na.rm=TRUE), lty=2)
	points(sort(m$R3), pch=21, bg=habitat.colors[m$habitat[order(m$R3)]])
	text(seq_len(nrow(m)), sort(m$R3), labels=m$taxon[order(m$R3)], pos=3, cex=0.3)

### Polly 2010 Figure 6 but for R2, instead of R3
	plot(sort(m$R2), pch=21, bg=habitat.colors[m$habitat], type="n")
	rect(xleft=-10, ybottom=mean(m$R2, na.rm=TRUE)-sd(m$R2, na.rm=TRUE), xright=150, ytop=mean(m$R2, na.rm=TRUE)+sd(m$R2, na.rm=TRUE), col=adjustcolor("black", alpha.f=0.25))
	abline(h=mean(m$R2, na.rm=TRUE), lty=2)
	points(sort(m$R2), pch=21, bg=habitat.colors[m$habitat[order(m$R2)]])
	text(seq_len(nrow(m)), sort(m$R2), labels=m$taxon[order(m$R2)], pos=3, cex=0.3)

### Polly 2010 Figure 6 but for A1, instead of R3
	plot(sort(m$A1), pch=21, bg=habitat.colors[m$habitat], type="n")
	rect(xleft=-10, ybottom=mean(m$A1, na.rm=TRUE)-sd(m$A1, na.rm=TRUE), xright=150, ytop=mean(m$A1, na.rm=TRUE)+sd(m$A1, na.rm=TRUE), col=adjustcolor("black", alpha.f=0.25))
	abline(h=mean(m$A1, na.rm=TRUE), lty=2)
	points(sort(m$A1), pch=21, bg=habitat.colors[m$habitat[order(m$A1)]])
	text(seq_len(nrow(m)), sort(m$A1), labels=m$taxon[order(m$A1)], pos=3, cex=0.3)

### Correlation of ratios with body mass
	this.ratio <- m$R2
	plot(m$BM, this.ratio, pch=21, bg=habitat.colors[m$habitat])
	text(m$BM, this.ratio, labels=m$taxon, pos=3, cex=0.3)
	abline(lm(this.ratio ~ m$BM))
	cor.test(m$BM, this.ratio)

	this.ratio <- m$R3
	plot(m$BM, this.ratio, pch=21, bg=habitat.colors[m$habitat])
	text(m$BM, this.ratio, labels=m$taxon, pos=3, cex=0.3)
	abline(lm(this.ratio ~ m$BM))
	cor.test(m$BM, this.ratio)
	
### PCA	################################################################
this.measurements <- data.frame(taxon=m$taxon, 
                                ast.length=rowMeans(m[,c("Length..medial.", "Length..lateral.")], na.rm=TRUE), 
                                ast.width=rowMeans(m[,c("Width..between.prox.condyles.", "Width..between.distal.condyles.navicular.facet.")], na.rm=TRUE), 
                                m$Calcaneoastragular.facet.to.cuboid.facet, 
                                m$Length.to.Calcaneoastragular.facet,
                                m$Length.to.Sustentacular.Process)

#2018/9/27 Need to remove 88929 and 47115 due to missing measurements for calcaneus
this.measurements <- this.measurements[!rownames(this.measurements) == "88929",]
this.measurements <- this.measurements[!rownames(this.measurements) == "47115",]

m.rem <- m
m.rem <- m.rem[!m.rem$Catalog.Number == "88929",]
m.rem <- m.rem[!m.rem$Catalog.Number == "47115",]
#remove rhinos
m.rem <- m.rem[!m.rem$Order == "Perissodactyla",]

this.scaled <- scale(this.measurements[,-1])

pca.log <- prcomp(log(this.measurements[,-1])) 
pca.scaled <- prcomp(this.scaled) 

pca <- pca.log
this.x <- pca$x[,2]
this.y <- pca$x[,3]
plot(this.x, this.y, pch=21, cex=1, bg=habitat.colors[m.rem$habitat]) #don't use m since i am removing things from this.measures (do removal before the parsign of data)
text(this.x, this.y, labels=this.measurements$taxon, pos=3, cex=0.3)

#pca with all measures
this.measurements <- data.frame(taxon=m$taxon, Amed.length=m$Length..medial., Alat.length=m$Length..lateral.,
                                Aprox.width=m$Width..between.prox.condyles.,
                                Adist.width=m$Width..between.distal.condyles.navicular.facet.,
                                m$Calcaneoastragular.facet.to.cuboid.facet, 
                                m$Length.to.Calcaneoastragular.facet)
#2018/9/27 Need to remove 88929 and 47115 due to missing measurements for calcaneus
this.measurements <- this.measurements[complete.cases(this.measurements),]

this.scaled <- scale(this.measurements[,-1])
complete.cases(this.measurements)
pca.log <- prcomp(log(this.measurements[,-1])) 
pca.scaled <- prcomp(this.scaled) 

pca <- pca.log
this.x <- pca$x[,2]
this.y <- pca$x[,3]
plot(this.x, this.y, pch=21, cex=1, bg=habitat.colors[m$habitat])
text(this.x, this.y, labels=m$taxon, pos=3, cex=0.3)
pca.log$rotation

#######################################################################################
###Tyler's Recomendation

#do factor analysis
#kernal local fisher's linear discriminant (lot of overlap in dist. to separate out)

##PCA
library(ggfortify)
dropCol <- c(1:3,12:14)
df <- log(m[,-dropCol])
habitat <- m$habitat

df <- scale(df)

df.pca <- prcomp(df)
plot(df.pca$x[,2], df.pca$x[,3], pch = 16, col = habitat.colors[m$habitat])

#PC1 vs PC2
autoplot(df.pca, data = m, colour = 'habitat', label = FALSE, label.size = 3, 
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3)

#PC2 vs PC3
autoplot(df.pca, data = m, x=2, y=3, choices = c(2:3), colour = 'habitat', 
         label = FALSE, label.size = 3, 
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3)

##Factor Analysis
d.factanal <- factanal(df, factors = 3, scores = 'regression')
#F1 and F2
autoplot(d.factanal, data = m, colour = 'habitat',
         loadings = TRUE, loadings.label = TRUE, loadings.label.size  = 3)

#F2 and F3
autoplot(d.factanal, x=2, y=3, data = m, colour = 'habitat',
         loadings = TRUE, loadings.label = TRUE, loadings.label.size  = 3)

#K-means
set.seed(1)
#set-symbols for pch
ksym <- gsub("OH",16, m$habitat)
ksym <- gsub("MH",17, ksym)
ksym <- gsub("CH",18, ksym)
ksym <- as.numeric(ksym)

df.kmean <- kmeans(df, 3)

kcol <- gsub(1, "red", df.kmean$cluster)
kcol <- gsub(2, "darkgreen", kcol)
kcol<- gsub(3, "blue", kcol)

autoplot(df.kmean, data = m, size = 0.1) +
  geom_point(shape = ksym, colour = kcol, size = 2) 

#Cluster
library(cluster)
df.clara <- clara(df, 3)
kcol <- gsub(1, "red", df.clara$clustering)
kcol <- gsub(2, "darkgreen", kcol)
kcol<- gsub(3, "blue", kcol)

autoplot(df.clara, data = m, size = 0.1) +
  geom_point(shape = ksym, colour = kcol, size = 2) 

df.fanny <- fanny(df, 3)
kcol <- gsub(1, "red", df.fanny$clustering)
kcol <- gsub(2, "darkgreen", kcol)
kcol<- gsub(3, "blue", kcol)

autoplot(df.fanny, data = m, frame = TRUE, size = 0.1)+
  geom_point(shape = ksym, colour = kcol, size = 2) 

df.pam <- pam(df, 3)
kcol <- gsub(1, "red", df.pam$clustering)
kcol <- gsub(2, "darkgreen", kcol)
kcol<- gsub(3, "blue", kcol)
autoplot(df.pam, data = m, frame = TRUE, frame.type = 'norm', size = 0.1)+
  geom_point(shape = ksym, colour = kcol, size = 2) 

#Kernal/Local Fisher Discriminant Analysis
library(lfda)
library(MASS)

lda.tst <- lda(x = as.matrix(df), grouping = m$habitat)
plot(lda.tst, col = habitat.colors[m$habitat])

# Local Fisher Discriminant Analysis (LFDA)
model1 <- lfda(df, m$habitat, 3, metric="plain")
autoplot(model1, data = m, frame = TRUE, frame.colour = 'habitat', size=0.1)+
  geom_point(shape = ksym) 
autoplot(model1, data = m, x=2, y=3,frame = TRUE, frame.colour = 'habitat', size=0.1)+
  geom_point(shape = ksym) 

model2 <- klfda(kmatrixGauss(df), m$habitat, 3, metric="plain")
autoplot(model2, data = m, frame = TRUE, frame.colour = 'habitat', size=0.1)+
  geom_point(shape = ksym) 
autoplot(model2, data = m,x=2,y=3, frame = TRUE, frame.colour = 'habitat', size=0.1)+
  geom_point(shape = ksym) 

# Semi-supervised Local Fisher Discriminant Analysis (SELF)
model3 <- self(df, m$habitat, beta = 0.1, r = 3, metric="plain")
autoplot(model3, data = m, frame = TRUE, frame.colour = 'habitat', size=0.1)+
  geom_point(shape = ksym) 
autoplot(model3, data = m, x=2,y=3, frame = TRUE, frame.colour = 'habitat', size=0.1)+
  geom_point(shape = ksym) 

############################################################################

require(phytools)
#read in tree
testTree <- read.nexus("/Users/emdoughty/Downloads/Appendix 3 alignment.nex")
testTreeTips <- testTree$tip.label

getTipTaxon <- function(treetips){
  tipList <- strsplit(treetips,split = "_")
  tipMat <- matrix(nrow = length(tipList), ncol = 2)
  for(ii in seq(1, length(tipList),1))
  {
    tipMat[ii,1] <- tipList[[ii]][1]
    tipMat[ii,2] <- tipList[[ii]][2]
  }
  
  labelTipsNew <- paste(tipMat[,1], tipMat[,2], sep="_")
  return (labelTipsNew)
}
labelTipsNew <- getTipTaxon(testTree$tip.label)
#how many entries are present per species name
new.table <- table(labelTipsNew)
removeVec <- vector()

for(yy in seq(1, length(new.table),1))
{
  if(new.table[yy] > 1)
  {
    names(new.table[yy])
    tip.options <- which(labelTipsNew == names(new.table[yy]))
    tip.remove <- sample(tip.options,size = length(tip.options)-1)
    removeVec <- append(removeVec, tip.remove)
  }
}

#check removeVec against names
testTree$tip.label[removeVec]

#collapse/remove species with multiple entries
testTree_crop <- drop.tip(testTree, removeVec)

#check if it worked
newNameList <- getTipTaxon(testTree_crop$tip.label)
table(newNameList)

testTree_crop$tip.label <- newNameList

#find which species are not on the tree
unique(m$taxon[!m$taxon %in% newNameList])
unique(testTree_crop$tip.label[!testTree_crop$tip.label %in% m$taxon])

#remove from each dataset after removing incompletes in dataset
m <- m[complete.cases(m),]
m.tree <- m[m$taxon %in% newNameList,]
df.tree <- df[rownames(df) %in% newNameList,]
testTree_crop <- drop.tip(testTree_crop,tip = testTree_crop$tip.label[!testTree_crop$tip.label %in% m.tree$taxon])

artTree <- testTree_crop
#align so that tips and m.tree are in same order
m.tree <- m.tree[artTree$tip.label,]
df.tree <- df.tree[artTree$tip.label,]

plot(artTree, cex = 0.05)
df.phylo.pca <- phyl.pca(artTree, df.tree, method="BM", mode="cov")
summary(df.phylo.pca)
plot(df.phylo.pca)
plot(df.phylo.pca$S[,1],df.phylo.pca$S[,2], col = habitat.colors[m.tree$habitat], 
     pch =16)
text(df.phylo.pca$S[,1],df.phylo.pca$S[,2], labels=artTree$tip.label, cex = 0.5)
#try supervised and unsupervied techniques
#pca, lda, factor anaysis, k-means clustering
#linear modeling

###############################################################################################
###############################################################################################


##phylo pca
#measure need to be indep or do size correction  (phyl.resid(body mass)) correct for all mrohp except max bm, used residuals in PCA
df.resid <- phyl.resid(artTree, df.tree[,"BM"], df.tree[,1:7])

#look at correlation of variables
#use correl of my matrix and correl test each var
#-1 is a difference matrix then plot it out
#want correlated variables for pca


# Visualizing a correlation matrix using Multidimensional Scaling
# MDS can be also used to reveal a hidden pattern in a correlation matrix.
#
# Correlation actually measures similarity, but it is easy to transform it to a 
# measure of dissimilarity. Distance between objects can be calculated as 1 - res.cor.

normal <- m
rownames(normal) <- NULL
res.cor <- cor(column_to_rownames(normal[,c(1,3:11)], var = "taxon")[-1], method = "spearman")
mds.cor <- (1 - res.cor) %>% cmdscale() %>% as_tibble()
colnames(mds.cor) <- c("dimension1", "dimension2")
ggscatter(mds.cor, x = "dimension1", y = "dimension2", size = 1, label = colnames(res.cor), repel = T)

autoplot(object = 1- res.cor)


###############################################################################################

m.R2 <- m[order(m$R2),]
m.R3 <- m[order(m$R3),]
m.A1 <- m[order(m$A1),]


quartz(height = 6, width = 10)
par(mar=c(7,6,3,3))
#par(mfrow=c(4,1))

#### Polly 2010 Figure 4c
# habitat.colors <- rainbow(length(unique(m$habitat)))
plot(m$R2, m$R3, xlim=c(0.35,0.65), ylim=c(1.2, 1.5), pch=21, cex=1, 
     bg=habitat.colors[m$habitat], xlab = "Gear Ratio (R2)", ylab = "Gear Ratio (R3)")
# plot(m$R2, m$R3, pch=21, cex=1, bg=habitat.colors[m$habitat])
#text(m$R2, m$R3, labels=m$taxon, pos=3, cex=0.5)
abline(lm(m$R3 ~ m$R2))

### Polly 2010 Figure 6
plot(sort(m$R3), pch=21, bg=habitat.colors[m$habitat], type="n", xlab=NA,xaxt="n",ylim=c(1,1.5), ylab = "Gear Ratio (R3)")
rect(xleft=-10, ybottom=mean(m$R3, na.rm=TRUE)-sd(m$R3, na.rm=TRUE), xright=150, ytop=mean(m$R3, na.rm=TRUE)+sd(m$R3, na.rm=TRUE), col=adjustcolor("black", alpha.f=0.25))
abline(h=mean(m$R3, na.rm=TRUE), lty=2)
points(sort(m$R3), pch=21, bg=habitat.colors[m$habitat[order(m$R3)]])
#text(seq_len(nrow(m)), sort(m$R3), labels=m$taxon[order(m$R3)], pos=1, cex=0.45, adj = 1,srt = 90, par("usr")[3]+1.6)
text(seq_len(nrow(m)), par("usr")[3]-0.07, labels = m$taxon[order(m$R3)], srt = 90, cex = 0.6, pos = 1, xpd = TRUE)

### Polly 2010 Figure 6 but for R2, instead of R3
plot(sort(m$R2), pch=21, bg=habitat.colors[m$habitat],xlab = NA, xaxt="n", type="n", ylab = "Gear Ratio (R2)")
rect(xleft=-10, ybottom=mean(m$R2, na.rm=TRUE)-sd(m$R2, na.rm=TRUE), xright=150, ytop=mean(m$R2, na.rm=TRUE)+sd(m$R2, na.rm=TRUE), col=adjustcolor("black", alpha.f=0.25))
abline(h=mean(m$R2, na.rm=TRUE), lty=2)
points(sort(m$R2), pch=21, bg=habitat.colors[m$habitat[order(m$R2)]])
text(seq_len(nrow(m)), par("usr")[3]-0.035, labels = m$taxon[order(m$R2)], srt = 90, cex = 0.6, pos = 1, xpd = TRUE)

### Polly 2010 Figure 6 but for A1, instead of R3
plot(sort(m$A1), pch=21, bg=habitat.colors[m$habitat], type="n", xlab = NA, xaxt="n", ylab = "Area of Astragalus (Aa)")
rect(xleft=-10, ybottom=mean(m$A1, na.rm=TRUE)-sd(m$A1, na.rm=TRUE), xright=150, ytop=mean(m$A1, na.rm=TRUE)+sd(m$A1, na.rm=TRUE), col=adjustcolor("black", alpha.f=0.25))
abline(h=mean(m$A1, na.rm=TRUE), lty=2)
points(sort(m$A1), pch=21, bg=habitat.colors[m$habitat[order(m$A1)]])
text(seq_len(nrow(m)), par("usr")[3]-0.1, labels = m$taxon[order(m$R2)], srt = 90, cex = 0.6, pos = 1, xpd = TRUE)



##################
#PCA
quartz(height = 6, width = 10)

autoplot(df.pca, data = m, colour = 'habitat', label = FALSE, label.size = 3, 
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3)+
  scale_color_manual(breaks = c("CH", "MH", "OH"),
                     values=c("dodgerblue", "chartreuse", "goldenrod1"))
autoplot(df.pca, data = m, x=2, y=3, choices = c(2:3), colour = 'habitat', 
               label = FALSE, label.size = 3, 
             loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3) +
  scale_color_manual(breaks = c("CH", "MH", "OH"),
                     values=c("dodgerblue", "chartreuse", "goldenrod1"))

#PHYLO
quartz(height = 6, width = 10)
plot(df.phylo.pca$S[,1],df.phylo.pca$S[,2], col = habitat.colors[m.tree$habitat], 
     pch =16, xlab = "PC1 (98.73%)", ylab = "PC2 (0.67 %)")
points(df.phylo.pca$L[,1:2])
pointLabel(df.phylo.pca$L[,2:3],             # set position of labels
           labels=rownames(df.phylo.pca$L),  # print labels
           cex=0.75 )

plot(df.phylo.pca$S[,2],df.phylo.pca$S[,3], col = habitat.colors[m.tree$habitat], 
     pch =16, xlab = "PC2 (0.67%)", ylab = "PC3 (0.45%)")
points(df.phylo.pca$L[,2:3])
pointLabel(df.phylo.pca$L[,2:3],             # set position of labels
           labels=rownames(df.phylo.pca$L),  # print labels
           cex=0.75 )
pointLabel(df.phylo.pca$S[,2:3],             # set position of labels
           labels=rownames(df.phylo.pca$S),  # print labels
          cex=0.75)

#loadings
quartz(height = 6, width = 10)
par(mfrow=c(2,2))


plot(df.phylo.pca$L[,1],df.phylo.pca$L[,2], xlab="PC1", ylab="PC2", main = "Phylogenetic PCA")
pointLabel(df.phylo.pca$L[,1:2],             # set position of labels
           labels=rownames(df.phylo.pca$L),  # print labels
           cex=0.75 )
plot(df.phylo.pca$L[,2],df.phylo.pca$L[,3],xlab="PC2", ylab="PC3")
pointLabel(df.phylo.pca$L[,2:3],             # set position of labels
           labels=rownames(df.phylo.pca$L),  # print labels
           cex=0.75 )

### LDA
quartz(height = 6, width = 10)
plot(lda.tst_MHCH, col = habitat.colors[m$habitat], xlim = c(-2.5,2.5))
plot(lda.tst_MHOH, col = habitat.colors[m$habitat], xlim = c(-2.5,2.5))

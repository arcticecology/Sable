
#The script to add norway to the conservative model because only RA is available


library(rioja)
library(analogue)
library(vegan)
library(labdsv)
setwd("C:\\Users\\Andrew\\Desktop\\Sable\\R\\surface")
species3<-read.csv(file="dataset_merge.csv", row.names=1) #read the csv file

species1<-species3[, -cbind(118:119) ] #this removes the taxa that are not good for modelling, undiff tanytarsini and such

species2 <- species1 / rowSums(species1) * 100

speciesX <- species2[, which(colSums(species2) != 0)]
species<-speciesX[, -cbind(116:117) ] #this removes the taxa that are not good for modelling, undiff tanytarsini and such

N <- 2
M <- 2
i <- colSums(speciesX >= M) >= N ## L1
species <- speciesX[, i, drop = FALSE]



tspecies <- decostand(species, method="hellinger") #transform data

#dca (b-diversity)
vare.dca <- decorana(tspecies)
vare.dca
summary(vare.dca)
plot(vare.dca)
plot(vare.dca, type = "n", xlim=c(-3,3), ylim=c(3,-3))
points(vare.dca, display = "sites", col = "black", cex = 1, pch = 21, bg = "red")
text(vare.dca, display="species", pos=2)


## load the data
coreinput_PD03<-read.csv(file="PD03m.csv", row.names=1)

coresp<-coreinput_PD03[, -cbind(1:3) ] #use this for cores
coresp1 <- coresp / rowSums(coresp) * 100


coresp2<-coresp1[, -cbind(117:118) ]



cols_to_keep <- intersect(colnames(species),colnames(coresp2))

PD03 <- coresp2[,cols_to_keep, drop=FALSE]

coreinput_PD67<-read.csv(file="PD67m.csv", row.names=1)

cores67<-coreinput_PD67[, -cbind(1:3) ] #use this for cores
coresp67 <- cores67 / rowSums(cores67) * 100

coresp677<-coresp67[, -cbind(117:118) ]

cols_to_keep <- intersect(colnames(species),colnames(coresp677))

PD67 <- coresp677[,cols_to_keep, drop=FALSE]


coreinput_PD37<-read.csv(file="PD37m.csv", row.names=1)
cores37<-coreinput_PD37[, -cbind(1:2) ] #use this for cores
coresp37 <- cores37 / rowSums(cores37) * 100
cols_to_keep <- intersect(colnames(PD67),colnames(coresp37))
PD37a <- coresp37[,cols_to_keep, drop=TRUE]




coreinput_PD38<-read.csv(file="PD38m.csv", row.names=1)
cores38<-coreinput_PD38[, -cbind(1:2) ] #use this for cores
coresp38 <- cores38 / rowSums(cores38) * 100
cols_to_keep <- intersect(colnames(PD67),colnames(coresp38))
PD38 <- coresp38[,cols_to_keep, drop=TRUE]





PD37 <-PD37a[ -cbind(9:15), ]
## Fit the timetrack ordination

core<-PD67
mod <- timetrack(species, core, transform = "hellinger",
                 method = "rda")
mod

## Plot the timetrack
plot(mod, ptype = "b", col = c("black", "red"), lwd = 2)

## Other options (reorder the time track)
ord <- rev(seq_len(nrow(core)))
plot(mod, choices = 2:3, order = ord, ptype = "b",
     col = c("forestgreen", "orange"), lwd = 2)


## scores and fitted methods
## IGNORE_RDIFF_BEGIN
head(fitted(mod, type = "passive"))
head(scores(mod, type = "passive"))
## IGNORE_RDIFF_END

my.pca <- rda(tspecies, scale=FALSE)
biplot(my.pca, scaling = 3, type = c("text", "points"))
text(my.pca, display = "sites", scaling = 1, cex = 0.8, col = "darkcyan")
## predict locations in timetrack for new observations

mod3 <- timetrack(species, PD03, transform = "hellinger",
                 method = "rda")
mod4 <- timetrack(species, PD67, transform = "hellinger",
                  method = "rda", scaling=3)
mod5 <- timetrack(species, PD37, transform = "hellinger",
                  method = "rda", scaling=3)
plot(mod5)

mod6 <- timetrack(species, PD38, transform = "hellinger",
                  method = "rda", scaling=3)
plot(mod6)

mod <- timetrack(species, core, transform = "hellinger",
                 method = "rda", labels=species)
mod
plot(mod)

## Plot the timetrack
plot(mod, ptype = "b", col = c("black", "red"), lwd = 2)
plot(mod3, ptype = "b", col = c("black", "red"), lwd = 2)
plot(mod4, ptype = "b", col = c("black", "red"), lwd = 2)
plot(mod5, ptype = "b", col = c("black", "red"), lwd = 2)
library(survival)
library(grid)
library(gridGraphics)

plot(mod, type = "n", ptype = "b", xlim=c(-0.6,0.5), ylim=c(-0.6,0.5))

# capture the plotted output as a grob of PD67
grid.echo()
grid.grab() -> k

# pull out the data from the grob..
k$children$`graphics-plot-1-points-1`$x -> x
k$children$`graphics-plot-1-points-1`$y -> y


plot(mod5, type = "n", ptype = "b", xlim=c(-0.6,0.5), ylim=c(-0.6,0.5))

# capture the plotted output as a grob of PD37
grid.echo()
grid.grab() -> k2

# pull out the data from the grob..
k2$children$`graphics-plot-1-points-1`$x -> x2
k2$children$`graphics-plot-1-points-1`$y -> y2

plot(mod6, type = "n", ptype = "b", xlim=c(-0.6,0.5), ylim=c(-0.6,0.5))

# capture the plotted output as a grob of PD38
grid.echo()
grid.grab() -> k3

# pull out the data from the grob..
k3$children$`graphics-plot-1-points-1`$x -> x3
k3$children$`graphics-plot-1-points-1`$y -> y3

#plot(mod, type = "n", ptype = "b")
#points(mod, which = "ordination", col = "grey", pch = 19, cex = 0.7)
#lines(x,y, col=1)
#points(x,y)
plot(mod3, xlim=c(-0.6,0.5), ylim=c(-0.6,0.5))
points(mod3, which = "passive", col = "red")
points(mod4, which = "passive", col = "blue")
points(mod, which = "ordination", col = "grey", pch = 19, cex = 0.7)
lines(x,y, col="blue", lty=2)

points(mod5, which = "passive", col = "orange", pch = 21, cex = 1)
lines(x2,y2, col="orange", lty=3)

points(mod6, which = "passive", col = "black", pch = 21, cex = 1)
lines(x3,y3, col="black", lty=4)

text(my.pca, display = "sites", scaling=3, pos=4, cex = 0.8, col = "darkcyan")
text(my.pca, display="species", scaling=3)
legend("bottomleft", inset =0.05,legend=c("PD03", "PD67", "PD37", "PD38"),
       col=c("red", "blue", "orange", "black"), lty=1:4, cex=0.8,
       box.lty=0)



tiff("SableSurface_b.tiff", width = 7, height = 7, units = 'in', res = 300)
plot(mod3, xlim=c(-0.6,0.5), ylim=c(-0.6,0.5), lwd=2)
points(mod3, which = "passive", col = "red")
points(mod4, which = "passive", col = "blue")
points(mod, which = "ordination", col = "grey", pch = 19, cex = 0.7)
lines(x,y, col="blue", lty=2, lwd=2)

points(mod5, which = "passive", col = "brown", pch = 21, cex = 1)
lines(x2,y2, col="brown", lty=3, lwd=3)

points(mod6, which = "passive", col = "black", pch = 21, cex = 1)
lines(x3,y3, col="black", lty=4, lwd=2)

text(my.pca, display = "sites", scaling=3, pos=4, cex = 0.8, col = "darkcyan")
#text(my.pca, display="species", scaling=3)
legend("bottomleft", inset =0.05,legend=c("PD03", "PD67", "PD37", "PD38"),
       col=c("red", "blue", "brown", "black"), lty=1:4, cex=0.8, lwd=2,
       box.lty=0)
dev.off()



#stratigraphy
#fos<-PD03
fos<-PD67
#fos<-PD37
#fos<-PD38

yPD03<-coreinput_PD03[1:2]
colnames(yPD03) <- c("Depth","Year")
yPD67<-coreinput_PD67[1:2]
colnames(yPD67) <- c("Depth","Year")
yPD37<-coreinput_PD37[1:2]
colnames(yPD37) <- c("Depth","Year")
yPD38<-coreinput_PD38[1:2]
colnames(yPD38) <- c("Depth","Year")

yDepth<-yPD67$Depth
yYear<- yPD67$Year


#yDepth<-yPD03$Depth
#yYear<- yPD03$Year

#yDepth<-yPD37$Depth
#yYear<- yPD37$Year

#yDepth<-yPD38$Depth
#yYear<- yPD38$Year

nms <- colnames(fos)
nms3<-gsub("\\.", " ", nms)

diss <- dist(sqrt(fos/100)^2)
clust <- chclust(diss, method="coniss")

diss <- dist(sqrt(fos/100)^2)

clust <- chclust(diss)

bstick(clust, 5)
library(NbClust)

library(ggpalaeo)
sp <- strat.plot(fos, scale.percent=TRUE, yvar=yDepth, 
                 title="PG Lake", ylabel="", srt.xlabel=45,  col.bar=1,
                 plot.line=FALSE, plot.bar = TRUE, lwd.bar=4,  y.rev=TRUE,
                 xLeft=0.18, xRight=0.75, y.axis=0.1, x.names=nms3, yTop=0.68, cex.yaxis=1,
                 cex.ylabel=1, cex.xlabel=1.0, cex.axis=1, clust=clust)


secondary_scale(sp, yvar = yDepth, yvar2 = yYear, n = 15, ylabel2 = "age (CE)")

addClustZone(sp, clust, 3, col="black", lty=2, lwd=2)


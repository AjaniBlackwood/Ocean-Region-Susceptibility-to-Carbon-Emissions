
##Ajani Blackwood
##Ocean Region Clustering R code

##This takes the lattitude and longitude of different locations within the ocean
## and clusters them based on sea surface salinity and temperature changes

#Clustering
library(ggplot2)

rm(list=ls(all=T))
dev.off()

pacman::p_load(pacman, rio) 
setwd('C:/Users/ajani/Documents/Data Analysis/mywork/project/datasets')

####Loading in the matched up climatology data###

climate_means = read.csv("climatology_matchedlabels.csv")

data_avg_in = read.csv('Model_meanSSS_SST_noNA.csv')

lat = data_avg_in[,1]
lon = data_avg_in[,2]

#This data_avg contains everything but the lat and lon information
data_avg = data_avg_in[, -c(1,2)]

#Individual non standardized average differences across models
diff_SST = data_avg[,1]
diff_SSS = data_avg[,2]

# standardize variables SST and SSS since units are different
diff_SST_standard = (diff_SST - mean(diff_SST))/sd(diff_SST) 
diff_SSS_standard = (diff_SSS - mean(diff_SSS))/sd(diff_SSS) 

#Make a new dataframe with standardized variables
data_avg_standard = data.frame(diff_SST_standard,diff_SSS_standard)


#WSSE, average silhouette width function to help select the numbers of clusters
set.seed(123)

library(cluster)
wsse <- function(k) {
  kmeans(data_avg_standard, k, nstart = 1000)$tot.withinss  # k-mean function
}

avg_sil <- function(kk) {
  hm <- agnes(data_avg_standard, diss=F, method="average") # construct the tree with ward's method
  labels=cutree(as.hclust(hm), k=kk)  # kk is the input number of clusters
  ss <- silhouette(labels, dist(data_avg_standard)) 
  return(mean(ss[, 3]))
}


##Make a range of values to test for the two clustering metrics
k.values <- 2:30

#testing out GAP statistics (for mor confidence when choosing number of clusters for WSSE)

gap_2 <- clusGap(data_avg_standard, FUN = kmeans,nstart = 50, K.max = 30, B = 20, d.power = 2,
                    spaceH0 = c("scaledPCA", "original"),
                    verbose = interactive())



plot(gap_2, type = "b", xlab = "k", ylab = expression(Gap[k]),
      main = "Ocean Region clustering Choices with Gap Statistics", do.arrows = TRUE,
     cex.axis = 1.5,
      arrowArgs = list(col="green3", length=1/16, angle=90, code=3))


library('purrr')  # purrr package contains the map function to execute functions with a range of input values (k,values)
#avg_sil_values <- map_dbl(k.values, avg_sil) # map_dbl execute the function avg_sil with inputs of k.values from 2 to n
#wsse_values <- map_dbl(k.values, wsse) # same as the above
#install.packages('mclust')

library(mclust) #mclust package contains the Gaussian Mixture model functions
###Use Mclust function to produce a whole range of solutions for cluster number 1 to 15

ts.mclust=Mclust(data_avg_standard, G = k.values)
###check out which cluster option has the biggest BIC value
maxBIC = apply(ts.mclust$BIC,1,max)

#pdf('Q1.pdf', width=7, height=4)
par(mar=c(4, 3, 4, 1), mgp=c(1.5, 0.3, 0), las=0) 
#par(mfrow=c(1,3))



##### Show a plot for WSSE Vs the number of clusters

plot(k.values, wsse_values, type='b', pch=16, 
      main ="Ocean Region Clustering Choices (WSSE)", col='red',xlab="Number of Clusters",ylab="WSSE", 
      cex.lab=1.2,cex.axis = 1.5)
# #plot(k.values, avg_sil_values,type='b', pch=16, col='blue',xlab="Cluster Number",ylab="Mean Silhouette Width", cex.lab=1.2)
plot(k.values, maxBIC, type='b', pch=16,col='cyan',
      main = "Ocean Region Clustering Choices (BIC)", 
     xlab="Number of Clusters",ylab="Maximum BIC", cex.lab=1.5, cex.axis = 1.75)

###### generate cluster labels based on two different methods

#Kmeans wants around 21 separations 
#try 11 from gap stats
#pm = pam(data_avg_standard, 21, diss = F)
pm = pam(data_avg_standard, 11, diss = F)
pm.labels = pm$cluster

#GMM wants about 23 separations or 11ish
#try 7 or 11 from slowing rate of change for BIC
#ts.mclust_choice=Mclust(data_avg_standard, G = 23)
ts.mclust_choice=Mclust(data_avg_standard, G = 11)
mm.labels = ts.mclust_choice$classification

library(fields)

#mm and pm labels ordered the same
Regions <- as.character(mm.labels)
Region_label_mm <- as.character(sort(unique(mm.labels)))
Region_label_pm <- as.character(sort(unique(pm.labels)))

#col_pm=fields::tim.colors(21)




######Plot the clusters for PAM and GMM

#par(mfrow=c(1,2))
plot(diff_SST_standard, diff_SSS_standard, bg=rainbow(11)[pm.labels], pch=23, 
     main = "PAM Representation for Ocean Region Separation (11 Regions)",
     xlab="Standardized Ocean Temperature Change",ylab="Standardized Ocean Salinity Change",
     xlim = c(-2,7),
     ylim = c(-8,3),
     cex.lab=1.2)
legend("bottom",legend = Region_label_pm, title = "Categorized Regions",
       col= rainbow(11) , lwd=5, cex=0.9, horiz = T)

par(mar=c(4, 5, 3, 1), mgp=c(3, 0.7, 0), las=0)
plot(diff_SST, diff_SSS, bg=rainbow(11)[pm.labels], pch=23, 
     main = "PAM Representation for Ocean Region Separation (11 Regions)",
     xlab="Ocean Temperature Change (Celsius)",ylab="Ocean Salinity Change (PPT)",
     xlim = c(1.4,4.1),
     ylim = c(-2.2,0.5),
     cex.lab=1.7,cex.axis = 1.7)
legend("bottom",legend = Region_label_pm, title = "Categorized Regions",
       col= rainbow(11) , lwd=5, cex=0.9, horiz = T)


plot(diff_SST_standard, diff_SSS_standard, bg=rainbow(11)[mm.labels], pch=23, 
     main = "GMM Representation for Ocean Region Separation (11 Regions)",
     xlab="Standardized Ocean Temperature Change",ylab="Standardized Ocean Salinity Change",
     xlim = c(-2,7),
     ylim = c(-8,3),
     cex.lab=1.2)
legend("bottom",legend = Region_label_mm, title = "Categorized Regions",
       col= rainbow(11), lwd=5, cex=0.9, horiz = T)

plot(diff_SST, diff_SSS, bg=rainbow(11)[mm.labels], pch=23, 
     main = "GMM Representation for Ocean Region Separation (11 Regions)",
     xlab="Ocean Temperature Change (Celsius)",ylab="Ocean Salinity Change (PPT)",
     xlim = c(1.4,4.1),
     ylim = c(-2.2,0.5),
     cex.lab=1.7,cex.axis = 1.7)
legend("bottom",legend = Region_label_pm, title = "Categorized Regions",
       col= rainbow(11) , lwd=5, cex=0.9, horiz = T)




###### Calculations for salinity and temperature
#Will use the original values not standardized
#variables are diff_SST and diff_SSS

### calculate the median for each cluster center
#Units are already in Celsius so no conversions needed
#Look at it for PM
SST_avg_pm =tapply(diff_SST, pm.labels, median)
SSS_avg_pm =tapply(diff_SSS, pm.labels, median)
Cluster_center_pm = data.frame(SST_avg_pm,SSS_avg_pm)

#Look at it for GMM
SST_avg_gmm =tapply(diff_SST, mm.labels, median)
SSS_avg_gmm =tapply(diff_SSS, mm.labels, median)
Cluster_center_gmm = data.frame(SST_avg_gmm,SSS_avg_gmm)

xmark = paste("Cluster", seq(1,11))
#par(mfrow=c(1,2))
#Barplots comparing Temperature and Salinity Change


label_x = "Region (Cluster) #"
label_temp = "Median Temperature Change (Celsius) "
label_salinity = "Median Salinity Change (PPT) "
plot_colors <- c("Blue","Orange")
text <- c("Median SST")
limit_temp = (c(-0.5,5))
limit_SSS = c(-2,0.5)
visibility = 1.4

###### Plot two separate plots for the temperature and salinity per method

#PAM temperature plot
barplot(t(as.matrix(Cluster_center_pm[1])), beside=F, col=c("Blue"),
        xlab = label_x , ylab = label_temp, main = "PAM Median SST for Each Cluster", 
        ylim = limit_temp,
        cex.axis = visibility, cex.main =  visibility, cex.lab = visibility,cex.names = visibility)
# legend("top",legend = text, text.width = max(sapply(text, strwidth)) *1.3,
#        col=plot_colors, lwd=7, cex=1.3, horiz = F)


#PAM salinity plot
barplot(t(as.matrix(Cluster_center_pm[2])), beside=F, col=c("Green"),
        xlab = label_x , ylab = label_salinity, main = "PAM Median SSS for Each Cluster", 
        ylim = limit_SSS,
        cex.axis = visibility, cex.main =  visibility, cex.lab = visibility,cex.names = visibility)

#GMM temperature plot
barplot(t(as.matrix(Cluster_center_gmm[1])), beside=F, col=c("Blue"),
        xlab = label_x , ylab = label_temp, main = "GMM Median SST for Each Cluster", 
        ylim = limit_temp,
        cex.axis = visibility, cex.main =  visibility, cex.lab = visibility,cex.names = visibility)

#GMM salinity plot
barplot(t(as.matrix(Cluster_center_gmm[2])), beside=F, col=c("Green"),
        xlab = label_x , ylab = label_salinity, main = "GMM Median SSS for Each Cluster", 
        ylim = limit_SSS,
        cex.axis = visibility, cex.main =  visibility, cex.lab = visibility,cex.names = visibility)


###### See how each cluster deviates from the grand median

change_rel_gmm = data.frame()
change_rel_pam = data.frame()

###two for loops calculate the deviation of each cluster center
###from the grand median
###one for loop for gmm one for PAM
for(i in 1:length(Cluster_center_gmm)){
  for(j in 1:length(Cluster_center_gmm[,1])){
    dev_gmm = abs(Cluster_center_gmm[j,i] - median(Cluster_center_gmm[,i]))
    change_rel_gmm[i,j] = dev_gmm
    
  }
}

for(i in 1:length(Cluster_center_pm)){
  for(j in 1:length(Cluster_center_pm[,1])){
    dev_pm = abs(Cluster_center_pm[j,i] - median(Cluster_center_pm[,i]))
    change_rel_pam[i,j] = dev_pm
    
  }
}

## relative changes for SST and SSS as calculated above
change_rel_gmm = t(change_rel_gmm)
change_rel_gmm = data.frame(change_rel_gmm)
colnames(change_rel_gmm) = c("Relative delta SST","Relative delta SSS")
rownames(change_rel_gmm) = c(1:11)

change_rel_pam = t(change_rel_pam)
change_rel_pam = data.frame(change_rel_pam)
colnames(change_rel_pam) = c("Relative delta SST","Relative delta SSS")
rownames(change_rel_pam) = c(1:11)

#The metric showing relative change




library("Rainbow")
#graphically show the relative changes and show the reference values
grand_SSTgmm = median(Cluster_center_gmm[,1])
barplot(t(as.matrix(change_rel_gmm[1])), beside=F, col=c("Blue"),
        xlab = label_x, ylab = "Absolute Deviation SST (Celsius)", 
        main = "Absolute Deviation from Grand Regional Median SST by Cluster (GMM)", 
        ylim = c(0,1.5), 
        cex.axis = visibility, cex.main =  visibility, cex.lab = visibility,cex.names = visibility)
mtext(bquote(paste("Grand Median SST (Celsius): "," = ", .(round(grand_SSTgmm,3)))),line=3)

grand_SSSgmm = median(Cluster_center_gmm[,2])
barplot(t(as.matrix(change_rel_gmm[2])), beside=F, col=c("Green"),
        xlab = label_x, ylab = "Absolute Deviation SSS (PPT)", 
        main = "Deviation from Grand Regional Median SSS by Cluster (GMM)", 
        ylim = c(0,1), 
        cex.axis = visibility, cex.main =  visibility, cex.lab = visibility,cex.names = visibility)
mtext(bquote(paste("Grand Median SSS (PPT): "," = ", .(round(grand_SSSgmm,3)))),line=3)

grand_SSTpm = median((Cluster_center_pm[,1]))
barplot(t(as.matrix(change_rel_pam[1])), beside=F, col=c("Blue"),
        xlab = label_x, ylab = "Absolute Deviation SST (Celsius)", 
        main = "Absolute Deviation from Grand Regional Median SST by Cluster (PAM)", 
        ylim = c(0,1.5), 
        cex.axis = visibility, cex.main =  visibility, cex.lab = visibility,cex.names = visibility)
mtext(bquote(paste("Grand Median SST (Celsius): "," = ", .(round(grand_SSTpm,3)))),line=3)

grand_SSSpm = median((Cluster_center_pm[,2]))
barplot(t(as.matrix(change_rel_pam[2])), beside=F, col=c("Green"),
        xlab = label_x, ylab = "Absolute Deviation SSS (PPT)", 
        main = "Absolute Deviation from Grand Regional Median SSS by Cluster (PAM)", 
        ylim = c(0,1), 
        cex.axis = visibility, cex.main =  visibility, cex.lab = visibility,cex.names = visibility)
mtext(bquote(paste("Grand Median SSS (PPT): "," = ", .(round(grand_SSSpm,3)))),line=3)



############the following shows how to plot cluster labels on a map
par(mar=c(4, 3, 4, 1), mgp=c(1.5, 0.3, 0), las=0) 
library('pracma')


########Plotting on world map

library("ggplot2")

nstall.packages(("ggspatial"))
library("ggspatial")

install.packages(("rnaturalearth"))
#require("rnaturalearth")
library("rnaturalearth")

install.packages(("rnaturalearthdata"))
require("rnaturalearthdata")
library("rnaturalearthdata")

install.packages(("rgeos"))
library("rgeos")


##### For PAM , Mapped labels

#limit <- max(abs(w1_loc$sst_avg)) * c(-1, 1)
limit <-  c(1, 11)

world <- ne_countries(scale = "medium", returnclass = "sf")

w1_loc = data.frame(pm.labels, lat, lon)

ggplot() + ggtitle ("Regional Separation by Cluster (PAM)") + 
  
  geom_raster( data = w1_loc , aes(x = lon, y = lat, fill = as.factor(pm.labels))) +
  
  geom_sf(fill="transparent", data = world) +
  
  coord_sf(xlim = c(min(lon), max(lon)), ylim = c(min(lat), max(lat))) +
  scale_fill_manual(name = "Separation", values = rainbow(11)) + 
  
  #scale_colour_gradientn(colors=rainbow(7)) +
  
  #scale_colour_manual(values = rainbow(21)) +
  
  #scale_fill_distiller(name = "Cluster #", palette = "Spectral", type = "div", limit = limit) +
  
  theme(text=element_text(size=18,  family="Arial"))


##### For GMM, Mapped labels

#limit <- max(abs(w1_loc$sst_avg)) * c(-1, 1)
limit <-  c(1, 11)

world <- ne_countries(scale = "medium", returnclass = "sf")

w1_loc = data.frame(mm.labels, lat, lon)

ggplot() + ggtitle ("Regional Separation by Cluster (GMM)") + 
  
  geom_raster( data = w1_loc , aes(x = lon, y = lat, fill = as.factor(mm.labels))) +
  
  geom_sf(fill="transparent", data = world) +
  
  coord_sf(xlim = c(min(lon), max(lon)), ylim = c(min(lat), max(lat))) +
  scale_fill_manual(name = "Separation", values = rainbow(11)) + 
  
  #scale_colour_gradientn(colors=rainbow(7)) +
  
  #scale_colour_manual(values = rainbow(21)) +
  
  #scale_fill_distiller(name = "Cluster #", palette = "Spectral", type = "div", limit = limit) +
  
  theme(text=element_text(size=18,  family="Arial")) 
  



###Before and after modeling for temperature and salinity###

##For PAM labels find the cluster centers for climatology (temperature) 
climatology_temp_pm = data.frame(climate_means[,3], pm.labels)
climatology_groupedtemp_pm = aggregate(x = climatology_temp_pm, 
                                    by = list(pm.labels), 
                                    FUN = "mean", na.rm = TRUE)

##For PAM labels find the cluster centers for climatology (salinity)
climatology_salinity_pm = data.frame(climate_means[,4], pm.labels)
climatology_groupedsalinity_pm = aggregate(x = climatology_salinity_pm, 
                                       by = list(pm.labels), 
                                       FUN = "mean", na.rm = TRUE)


##For GMM labels find the cluster centers for climatology (temperature) 
climatology_temp_gmm = data.frame(climate_means[,3], mm.labels)
climatology_groupedtemp_gmm = aggregate(x = climatology_temp_gmm, 
                                       by = list(mm.labels), 
                                       FUN = "mean", na.rm = TRUE)

##For GMM labels find the cluster centers for climatology (salinity)
climatology_salinity_gmm = data.frame(climate_means[,4], mm.labels)
climatology_groupedsalinity_gmm = aggregate(x = climatology_salinity_gmm, 
                                           by = list(mm.labels), 
                                           FUN = "mean", na.rm = TRUE)




#####  Find the "after" temperatures and salinity values for PAM
after_temp_pm = c()
after_salinity_pm = c()
for(i in 1:length(climatology_groupedtemp_pm[,2])){#length is the same regardless
  
  after_temp_pm[i] = climatology_groupedtemp_pm[,2][i] + Cluster_center_pm[,1][i]
  after_salinity_pm[i] = climatology_groupedsalinity_pm[,2][i] + Cluster_center_pm[,2]
}
after_temp_pm = data.frame(after_temp_pm)
after_salinity_pm = data.frame(after_salinity_pm)


######  Find the "after" temperatures and salinity values for GMM
after_temp_gmm = c()
after_salinity_gmm = c()
for(i in 1:length(climatology_groupedtemp_gmm[,2])){#length is the same regardless
  
  after_temp_gmm[i] = climatology_groupedtemp_gmm[,2][i] + Cluster_center_gmm[,1][i]
  after_salinity_gmm[i] = climatology_groupedsalinity_gmm[,2][i] + Cluster_center_gmm[,2]
}
after_temp_gmm = data.frame(after_temp_gmm)
after_salinity_gmm = data.frame(after_salinity_gmm)


#######put before and after in a dataframe
temp_before_after_pm = data.frame(climatology_groupedtemp_pm[,2],after_temp_pm)
temp_before_after_gmm = data.frame(climatology_groupedtemp_gmm[,2], after_temp_gmm)
salinity_before_after_pm = data.frame(climatology_groupedsalinity_pm[,2], after_salinity_pm)
salinity_before_after_gmm = data.frame(climatology_groupedsalinity_gmm[,2], after_salinity_gmm)

par(mar=c(4, 4, 2, 1))
par(mfrow=c(1,2))
xtick = seq(2, 11*3, by = 3)
xmark = paste("Cluster", seq(1,11))
text <- c("Climatology Temperature", "Predicted Temperature")
barplot(t(as.matrix(temp_before_after_pm)), beside=T, col=c("Blue","Orange"),
        xlab = label_x, ylab = "Sea Surface Temperature (Celsius)", 
        main = "Predicted Future SST by Cluster (PAM)", ylim = c(0,36),
        cex.axis = visibility, cex.main =  visibility, 
        cex.lab = visibility,cex.names = visibility)
abline(h = 27,col = "Red")
axis(side=1, at=xtick, labels=xmark, col='black', cex.axis=1., tck=-0.02)
legend("top",legend = text, text.width = max(sapply(text, strwidth))*1.3,
       col=c("Blue","Orange"), lwd=7, cex=1.3, horiz = F)


barplot(t(as.matrix(temp_before_after_gmm)), beside=T, col=c("Blue","Orange"),
        xlab = label_x, ylab = "Sea Surface Temperature (Celsius)", main = "GMM Method", ylim = c(0,36),
        cex.axis = visibility, cex.main =  visibility, 
        cex.lab = visibility,cex.names = visibility)
abline(h = 27,col = "Red")
axis(side=1, at=xtick, labels=xmark, col='black', cex.axis=1., tck=-0.02)
legend("top",legend = text, text.width = max(sapply(text, strwidth))*1.3,
       col=c("Blue","Orange"), lwd=7, cex=1.3, horiz = F)


###### show plots for salinity before and after

text_salinity <- c("Climatology Salinity", "Predicted Salinity")
barplot(t(as.matrix(salinity_before_after_pm)), beside=T, col=c("Green","Yellow"),
        xlab = label_x, ylab = "Sea Surface Salinity (PPT)", 
        main = "Predicted Future SSS by Cluster (PAM)", ylim = c(0,45),
        cex.axis = visibility, cex.main =  visibility, 
        cex.lab = visibility,cex.names = visibility)
#abline(h = 27,col = "Red")
axis(side=1, at=xtick, labels=xmark, col='black', cex.axis=1., tck=-0.02)
legend("top",legend = text_salinity, text.width = max(sapply(text_salinity, strwidth))*1.3,
       col=c("Green","Yellow"), lwd=7, cex=1.3, horiz = F)


barplot(t(as.matrix(salinity_before_after_gmm)), beside=T, col=c("Green","Yellow"),
        xlab = label_x, ylab = "Salinity (PPT)", main = "GMM Method", ylim = c(0,45),
        cex.axis = visibility, cex.main =  visibility, 
        cex.lab = visibility,cex.names = visibility)
#abline(h = 27,col = "Red")
axis(side=1, at=xtick, labels=xmark, col='black', cex.axis=1., tck=-0.02)
legend("top",legend = text_salinity, text.width = max(sapply(text_salinity, strwidth))*1.3,
       col=c("Green","Yellow"), lwd=7, cex=1.3, horiz = F)







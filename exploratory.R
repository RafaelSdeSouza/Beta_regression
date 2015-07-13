# Exploration of the dataset
require(AMADA)
data.1= read.table(file="FiBY_escape_data_all.dat",header=FALSE)
colnames(data.1)<-c("redshift","fEsc","Mvir","Mstar","Mgas","QHI","sfr_gas",
                    "sfr_stars","ssfr_gas","ssfr_stars","baryon_fraction",
                    "spin","age_star_mean","age_star_max","age_star_min")

data.2<-data.1[,-2]

cor1<-Corr_MIC(data.2,method="pearson")
plotgraph(cor1)

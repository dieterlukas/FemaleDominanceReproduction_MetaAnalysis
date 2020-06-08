#Meta analysis of effect of rank on reproductive success in female social mammals


library(metafor)
library(ape)
library(dplyr)
library(readr)  

#There are five parts: 1) loading the data; 2) analyses with metafor; 3) analyses with rethinking; 4) figures; 5) simulation to show differences in effect sizes when assuming linear hierarchy versus castes


# 1) loading the data
# all data are obtained from online sources


# Simulated data for testing the code for the preregistration
# Load simulated data from GitHub
simulateddata<-read_csv(url("https://raw.githubusercontent.com/dieterlukas/FemaleDominanceReproduction_MetaAnalysis/master/SimulatedData_MetaAnalysis_FemaleDominanceReproduction_March2020.csv"))
simulateddata<-data.frame(simulateddata)
dat<-simulateddata
#In the simulated data, there are 57 effect sizes from 26 studies on 21 species. All effects and their variances have been randomly generated; effects from supposedly the same study were assumed to be similar, but we do not expect a species signal or phylogenetic signal

# For the actual analyses, we will use the full data file
# Load input data from GitHub
inputdata<-read_csv(url("https://raw.githubusercontent.com/dieterlukas/FemaleDominanceReproduction_MetaAnalysis/master/InputData_MetaAnalysis_FemaleDominanceReproduction_March2020.csv"))
inputdata<-data.frame(inputdata)
dat<-inputdata
#In the actual data, there are 444 effect sizes from 187 studies on 86 species. No analyses have yet been performed with these data (01 April 2020)


#These are the variable names
#[1] "EffectRef"            "StudyRef"             "Latin"                "SpeciesRef"          
#[5] "MeasureType_Category" "n"                    "r_value"              "Fishers_Z_r"         
#[9] "variance"             "No_of_years"          "Populationtype"       "GroupSize"           
#[13] "StudyPopulation"      "Rank_Cont_Cat"        "RankDetermined"    

# Phylogenetic tree
# For the actual data, we will construct a consensus tree based on the species in the final dataset from the mammalian phylogeny at vertlife.org/phylosubsets/  ; the nexus file for the tree will be stored at github
# For the fake data, we are using the phylogeny provided by Rolland et al. 2014 Plos Biology
mammaltree<-read.tree(url("https://journals.plos.org/plosbiology/article/file?id=10.1371/journal.pbio.1001775.s001&type=supplementary"))




# 2) analyses building models with the package metafor

#based model: determine the overall effect size, accounting for the measurement error as reflected by the sampling variance
meta_simple <- rma.mv(yi=Fishers_Z_r, V=variance, data=dat)
meta_simple


#Test for publication bias
#We plot the effect sizes over the inverse of their variances (which reflect their precision) to determine whether there are signs of publication bias
#These signs would be evident if at low precision (towards the bottom of the y-axis) the effect size estimates tend to only show extreme values
#For the fake data we do not expect such a signal since all data were randomly generated. In fact, we expect the opposite, in that we will see a large spread of effect sizes even at high precision because the effect sizes were generated randomly. 
funnel(res_all, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), 
       yaxis="seinv",xlab="effect size (Zr)",ylab="precision (1/SE)",digits=c(1,0),xlim=c(-1,1),ylim=c(3,22))


#Determine the overall effect size in a model that assumes that measurements from the same study, the same population, and the same species are not independent and that each measurement might have a unique approach.
#For the fake data, we do not expect any dependence by study (because the sample size is small) or by population/species (because data were randomly generated)
meta_randomfactors <- rma.mv(yi=Fishers_Z_r, V=variance, data=dat, method="REML", random=list(~1|StudyRef,~1|StudyPopulation,~1|SpeciesRef,~1|EffectRef))
meta_randomfactors

#Check potential dependence from each of the random factors separately
#model with species as random term
meta_species <- rma.mv(yi=Fishers_Z_r, V=variance, data=dat, slab=StudyRef, random=~1|SpeciesRef)
meta_species 

#model with study as random term
meta_study <- rma.mv(yi=Fishers_Z_r, V=variance, data=dat, slab=EffectRef, random=~ 1|StudyRef)
meta_study 

#model with population as random term
meta_population <- rma.mv(yi=Fishers_Z_r, V=variance, data=dat, slab=EffectRef, random=~ 1|StudyPopulation)
meta_population 


#Build a model that incorporates potential similarites in effect sizes due to shared phylogenetic history among species

#First we need to match the species in the data to the species included in the phylogenetic tree. 
#The files mtree and mdata will contain the phylogeny and data for the set of species included both in the study and the phylogeny
spp_obs<-unique(dat$Latin)
spp_obs<-matrix(nrow=length(spp_obs))
rownames(spp_obs)<-unique(dat$Latin)

missing<-treedata(mammaltree,data=spp_obs,warnings=FALSE)
mtree<-missing$phy
speciesnames<-mtree$tip.label
mdata<-dat[dat$Latin %in% speciesnames,]

#Next, we convert the phylogeny into a covariance matrix that reflects the shared phylogenetic history among all pairs of species in the sample
#There are two approaches: first assuming a Brownian Motion model of evolution, and next a model that includes the parameter Pagel's lambda to account for potential changes in the rate of evolution
Rbm<-corBrownian( phy = mtree )
Rlambda<-corPagel(0.75,phy=mtree,fixed=FALSE)
VCV_bm <- vcv(Rbm)
VCV_pagel<-vcv(Rlambda)

#We fit the phylogenetic model, including also the potential independence from study reference and the potential for each effect to have a unique measurement.
#First the model using Brownian motion
meta_Brownian <- rma.mv(yi=Fishers_Z_r, V=variance, 
                  random = list(~1 | Latin, ~1 | StudyRef,  ~1 | EffectRef),
                  R = list(Latin = VCV_bm),
                  method = "REML", 
                  data = mdata)
meta_Brownian


#Second the model with Pagel's lambda
meta_Lambda <- rma.mv(yi=Fishers_Z_r, V=variance, 
                  random = list(~1 | Latin, ~1 | StudyRef,  ~1 | EffectRef),
                  R = list(Latin = VCV_pagel),
                  method = "REML", 
                  data = mdata)
meta_Lambda

# We now built the models to test for the influence of potential moderators on the effect sizes
# To test the predictions, we will change the moderator to match the respective variable
# We will test one variable at a time to assess the different predictions


#Models without taking phylogeny into account

#model assessing the influence of difference on group sizes on the measured effect sizes, with study reference as random term
meta_groupsize <- rma.mv(yi=Fishers_Z_r, V=variance, data=dat, mods=~ GroupSize, random=~1| StudyRef)
meta_groupsize

#model assessing whether rank was measured categorical or continuous influences the effect sizes, with the various factors that could lead to dependency included 
meta_rankmeasurement <- rma.mv(yi=Fishers_Z_r, V=variance, data=dat, mods=~ Rank_Cont_Cat, random=list(~1|StudyRef,~1|StudyPopulation,~1|SpeciesRef,~1|EffectRef))
meta_rankmeasurement

#Models taking shared phylogenetic history into account

#model assessing whether rank was measured categorical or continuous influences the effect sizes, with the phylogenetic similarity among species added
meta_Lambda_rankmeasurement <- rma.mv(yi=Fishers_Z_r, V=variance, mods=~ Rank_Cont_Cat,
                      random = list(~1 | Latin, ~1 | StudyRef,  ~1 | EffectRef),
                      R = list(Latin = VCV_pagel),
                      method = "REML", 
                      data = mdata)
meta_Lambda_rankmeasurement


# 3) analyses building models with the package rethinking. For details on how to install the package see: https://github.com/rmcelreath/rethinking

library(rethinking)
library(ape)
library(geiger)

#If some variables end up having missing data, we will filter these cases out before any analysis
dstan <- dat[ complete.cases(dat$Rank_Cont_Cat,dat$GroupSize),]

#A non-phylogenetic model to determine the overall effect size
#We create the structured data file for the model 
dat_list_simple <- list(
  N_spp = nrow(dstan),
  Z_obs = dstan$Fishers_Z_r,
  Z_sd = dstan$variance
)

#We run the model assuming that the sampling variance reflects an error in the measurements of the true effect sizes
rethinking_simple <- ulam(
  alist(
    Z_obs ~dnorm ( Z_true , Z_sd),
    vector[N_spp]:Z_true ~ dnorm( mu , sigma_sq ),
    mu ~ normal(0,1),
    sigma_sq~dexp(1)
  ),data=dat_list_simple)
precis(rethinking_simple)


#A non-phylogenetic model to determine the influence of group size on the variation in effect sizes
#We create the structured data file for the model 
dat_list_groupsize <- list(
  N_spp = nrow(dstan),
  Z_obs = dstan$Fishers_Z_r,
  Z_sd = dstan$variance, 
  G = standardize(standardize(log(dstan$GroupSize)))
)

#We run the model
rethinking_groupsize <- ulam(
  alist(
    Z_obs ~dnorm ( Z_true , Z_sd),
    vector[N_spp]:Z_true ~ dnorm( mu , sigma_sq ),
    mu <- a+bG*G,
    a~normal(0,1),
    bG~normal(0,0.5),
    sigma_sq~dexp(1)
  ),data=dat_list_groupsize)
precis(rethinking_groupsize)


#Next we prepare for the phylogenetic models
# We again match the species in the dataset to the phylogenetic tree
spp_obs<-unique(dstan$Latin)
spp_obs<-matrix(nrow=length(spp_obs))
rownames(spp_obs)<-unique(dstan$Latin)

missing<-treedata(mammaltree,data=spp_obs,warnings=FALSE)
mtree<-missing$phy
speciesnames<-mtree$tip.label
mdata<-dstan[dstan$Latin %in% speciesnames,]

#Next, we convert the phylogeny into a covariance matrix that reflects the shared phylogenetic history among all pairs of species in the sample
#There are two approaches: first assuming a Brownian Motion model of evolution, and next a model that includes the parameter Pagel's lambda to account for potential changes in the rate of evolution
Rbm<-corBrownian( phy = mtree )
Rlambda<-corPagel(0.75,phy=mtree,fixed=FALSE)
VCV_bm <- vcv(Rbm)
VCV_pagel<-vcv(Rlambda)


#We create the structured data file for the model

#example model to test for the influence of group size on effect sizes
dat_list <- list(
  N_spp = nrow(mdata),
  Z_obs = mdata$Fishers_Z_r,
  Z_sd = mdata$variance, 
  G = standardize(standardize(log(mdata$GroupSize)))
)

#The shared phylogenetic history between pairs of species is reflected in a covariance matrix
#The covariance matrix is standardized by the maximum distance from the shared common ancestor to the tips
#The diagonals (comparison of a population with itself) therefore has the value 1
#Two species that are in separate lineages that split off at the very first node after the common ancestor have entries of 0 in the matrix
#In meta-analyses, we usuallyhave repeated observations from the same species. 
#In this model, we assume that all of these to have the same expected similarity of 0.99. That is, we do not further differentiate whether observations come from the same or from different populations.
dat_list$V <- VCV_bm[ mdata$Latin,mdata$Latin ]
dat_list$R <-dat_list$V / max(VCV_bm)
dat_list$R<-dat_list$R*0.99
diag(dat_list$R)<-1

#We run the model.
rethinking_phylogenetic_groupsize <- ulam(
  alist(
    Z_obs ~dnorm ( Z_true , Z_sd),
    vector[N_spp]:Z_true ~ multi_normal( mu , SIGMA ),
    mu <- a+bG*G,
    matrix[N_spp,N_spp]: SIGMA <- R * sigma_sq,
    a~normal(0,1),
    bG~normal(0,0.5),
    sigma_sq~exponential(1)
  ),data=dat_list)
#For the simulated data, there is no phylogenetic effect. Accordingly, model indicates a large error for the species-level estimates: they are supposed to be more similar than they actually are.
precis(rethinking_phylogenetic_groupsize)

#We run another model for comparison on the restricted dataset (species included in phylogeny) without including the shared phylogenetic history.
rethinking_nonphylogenetic_groupsize <- ulam(
  alist(
    Z_obs ~dnorm ( Z_true , Z_sd),
    vector[N_spp]:Z_true ~ dnorm( mu , sigma_sq ),
    mu <- a+bG*G,
    a~normal(0,1),
    bG~normal(0,0.5),
    sigma_sq~dexp(1)
  ),data=dat_list)
precis(rethinking_nonphylogenetic_groupsize)


#Next, we also want to control for the potential similarity whether or not samples come from the same population within a species
populationmatrix<-dist(mdata$StudyPopulation,diag=TRUE,upper=TRUE)
populationmatrix[populationmatrix>0]<-1
populationmatrix<-1-populationmatrix
populationmatrix<-0.99*populationmatrix
populationmatrix<-as.matrix(populationmatrix)
diag(populationmatrix)<-1
dat_list$P<-populationmatrix

rethinking_phylogenetic_population_groupsize <- ulam(
  alist(
    Z_obs ~dnorm ( Z_true , Z_sd),
    vector[N_spp]:Z_true ~ multi_normal( mu , SIGMA ),
    mu <- a+bG*G,
    matrix[N_spp,N_spp]: SIGMA <- R * sigma_sq + P * sigma_pop,
    a~normal(0,1),
    bG~normal(0,0.5),
    sigma_sq~exponential(1),
    sigma_pop~exponential(1)
  ),data=dat_list)

precis(rethinking_phylogenetic_population_groupsize)


#Instead of the Brownian motion model, use a general Gaussian process. This is also helpful to show that there is no phylogenetic signal in the simulated dataset
Dmat<-cophenetic(mtree)
dat_list$Dmat<-Dmat[ mdata$Latin,mdata$Latin ]/max(Dmat)
rethinking_phylogenetic_gaussian_groupsize <- ulam(
  alist(
    Z_obs ~dnorm ( Z_true , Z_sd),
    vector[N_spp]:Z_true ~ multi_normal( mu , SIGMA ),
    mu <- a+bG*G,
    matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
    a~normal(0,1),
    bG~normal(0,0.5),
    etasq~exponential(1),
    rhosq~exponential(1)
  ),data=dat_list,chains=4)
precis(rethinking_phylogenetic_gaussian_groupsize)

#To check the phylogenetic signal, we plot the estimate covariance(similarity) in effect sizes as species are more distantly related
#A phylogenetic signal in the plot would be shown by high covariances among species that have small phylogenetic distances and declining covariances as phylogenetic distance increases
#For the simulated dataset, most of the estimates show a covariance of 0, and no change as phylogenetic distance increases
post<-extract.samples(rethinking_phylogenetic_groupsize)
plot(NULL,xlim=c(0,max(dat_list$Dmat)),ylim=c(0,1.5),xlab="phylogenetic distance",ylab="covariance")
for (i in 1:50) curve(post$etasq[i]*exp(-post$rhosq[i]*x^2),add=TRUE,col=grau())




#Here, we expand the phylogenetic model to include a nested effect
#In the analyses, we might expect that patterns differ between cooperative and plural breeders, or that group size does not have an influence in captive studies where individuals receive food but in wild populations
#Here, we estimate whether patterns differ between studies that measured rank on a linear hierarchy and those that classified individuals into rank categories.
#In particular, we assume that the influence of group size on the effect sizes might be different depending on the measurement method.
dat_list <- list(
  N_spp = nrow(mdata),
  Z_obs = mdata$Fishers_Z_r,
  Z_sd = mdata$variance, 
  G = standardize(standardize(log(mdata$GroupSize))),
  R = as.integer((mdata$Rank_Cont_Cat=="Categorical")+1)
)
Dmat<-cophenetic(mtree)
dat_list$Dmat<-Dmat[ mdata$Latin,mdata$Latin ]/max(Dmat)

rethinking_phylogenetic_gaussian_groupsize <- ulam(
  alist(
    Z_obs ~dnorm ( Z_true , Z_sd),
    vector[N_spp]:Z_true ~ multi_normal( mu , SIGMA ),
    mu <- a + bG*G*r[R],
    r[R] ~ dnorm(0,1),
    matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
    a~normal(0,1),
    bG~normal(0,0.5),
    etasq~exponential(1),
    rhosq~exponential(1)
  ),data=dat_list,chains=4)
precis(rethinking_phylogenetic_gaussian_groupsize)



dat_list <- list(
  N_spp = nrow(mdata),
  Z_obs = mdata$Fishers_Z_r,
  Z_sd = mdata$variance, 
  G = standardize(standardize(log(mdata$GroupSize))),
  M = as.integer(as.factor(dat$MeasureType_Category))
)
Dmat<-cophenetic(mtree)
dat_list$Dmat<-Dmat[ mdata$Latin,mdata$Latin ]/max(Dmat)

rethinking_phylogenetic_gaussian_groupsize <- ulam(
  alist(
    Z_obs ~dnorm ( Z_true , Z_sd),
    vector[N_spp]:Z_true ~ multi_normal( mu , SIGMA ),
    mu <- a + bG*G*r[M],
    r[M] ~ dnorm(0,1),
    matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01),
    a~normal(0,1),
    bG~normal(0,0.5),
    etasq~exponential(1),
    rhosq~exponential(1)
  ),data=dat_list,chains=4)
precis(rethinking_phylogenetic_gaussian_groupsize)







#Next, we expand the phylogenetic Gaussian model to also include the population identiy
#Instead of the Brownian motion model, use a general Gaussian process. We can also use this more flexible model for the other distance matrices
populationmatrix<-dist(mdata$StudyPopulation,diag=TRUE,upper=TRUE)
populationmatrix[populationmatrix>0]<-1
populationmatrix<-as.matrix(populationmatrix)
diag(populationmatrix)<-0
dat_list$Pmat<-populationmatrix

rethinking_phylogenetic_population_gaussian_groupsize <- ulam(
  alist(
    Z_obs ~dnorm ( Z_true , Z_sd),
    vector[N_spp]:Z_true ~ multi_normal( mu , SIGMA ),
    mu <- a+bG*G,
    matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01)+cov_GPL2( Pmat, eta1sq, rho1sq, 0.01),
    a~normal(0,1),
    bG~normal(0,0.5),
    etasq~exponential(1),
    rhosq~exponential(1),
    eta1sq~exponential(1),
    rho1sq~exponential(1)
  ),data=dat_list,chains=4)
precis(rethinking_phylogenetic_population_gaussian_groupsize)

#To check the phylogenetic signal, we plot the estimate covariance(similarity) in effect sizes as species are more distantly related
#A phylogenetic signal in the plot would be shown by high covariances among species that have small phylogenetic distances and declining covariances as phylogenetic distance increases
#For the simulated dataset, most of the estimates show a covariance of 0, and no change as phylogenetic distance increases
post<-extract.samples(rethinking_phylogenetic_population_gaussian_groupsize)
plot(NULL,xlim=c(0,max(dat_list$Dmat)),ylim=c(0,1.5),xlab="phylogenetic distance",ylab="covariance")
for (i in 1:50) curve(post$etasq[i]*exp(-post$rhosq[i]*x^2),add=TRUE,col=grau())

post<-extract.samples(rethinking_phylogenetic_population_gaussian_groupsize)
plot(NULL,xlim=c(0,max(dat_list$Pmat)),ylim=c(0,1.5),xlab="population distance",ylab="covariance")
for (i in 1:50) curve(post$eta1sq[i]*exp(-post$rho1sq[i]*x^2),add=TRUE,col=grau())







#We can also expand this approach for other variables that might be shared, for example how rank was measured
#Instead of the Brownian motion model, use a general Gaussian process. We can also use this more flexible model for the other distance matrices

rankmatrix<-dist(as.numeric(mdata$Rank_Cont_Cat=="Categorical"),diag=TRUE,upper=TRUE)
rankmatrix<-as.matrix(rankmatrix)
dat_list$Cmat<-rankmatrix

rethinking_phylogenetic_rank_gaussian_groupsize <- ulam(
  alist(
    Z_obs ~dnorm ( Z_true , Z_sd),
    vector[N_spp]:Z_true ~ multi_normal( mu , SIGMA ),
    mu <- a+bG*G,
    matrix[N_spp,N_spp]: SIGMA <- cov_GPL2( Dmat, etasq, rhosq, 0.01)+cov_GPL2( Cmat, eta1sq, rho1sq, 0.01),
    a~normal(0,1),
    bG~normal(0,0.5),
    etasq~exponential(1),
    rhosq~exponential(1),
    eta1sq~exponential(1),
    rho1sq~exponential(1)
  ),data=dat_list,chains=4)
precis(rethinking_phylogenetic_rank_gaussian_groupsize)

#To check the signal of similarity according to the use of the same approach to measure rank, we plot the estimate covariance(similarity) in effect sizes as species are either the same or different
#A  signal in the plot would be shown by high covariances among measurements that used the same method (distance of zero) and low covariances when the method differs (distance of 1)
#For the simulated dataset, most of the estimates show a covariance of 0, and no change 

post<-extract.samples(rethinking_phylogenetic_rank_gaussian_groupsize)
plot(NULL,xlim=c(0,max(dat_list$Cmat)),ylim=c(0,1.5),xlab="population distance",ylab="covariance")
for (i in 1:50) curve(post$eta1sq[i]*exp(-post$rho1sq[i]*x^2),add=TRUE,col=grau())








#To show that in the simulated dataset the lack of a phylogenetic signal is not only due to the high variation within species
#but also due to the absence of closely related species being similar, we run another model taking only a single (here the first) entry for each species
mdata_small<-mdata[c(1,6,7,8,12,14,16,17,19,24,36,37,40,42,43,45,46,49,52),]

dat_list <- list(
  N_spp = nrow(mdata_small),
  Z_obs = mdata_small$Fishers_Z_r,
  Z_sd = mdata_small$variance, 
  G = standardize(standardize(log(mdata_small$GroupSize)))
)

dat_list$V <- VCV_bm[ mdata_small$Latin,mdata_small$Latin ]
dat_list$R <-dat_list$V / max(VCV_bm)

rethinking_phylogenetic_groupsize_oneperspp <- ulam(
  alist(
    Z_obs ~dnorm ( Z_true , Z_sd),
    vector[N_spp]:Z_true ~ multi_normal( mu , SIGMA ),
    mu <- a+bG*G,
    matrix[N_spp,N_spp]: SIGMA <- R * sigma_sq,
    a~normal(0,1),
    bG~normal(0,0.5),
    sigma_sq~exponential(1)
  ),data=dat_list)
precis(rethinking_phylogenetic_groupsize_oneperspp)







# 4) Visualisation

# We use the package tidymeta because of it's flexibility
# more info: https://github.com/malcolmbarrett/tidymeta

library(tidymeta)
library(ggplot2)
library(dplyr)
library(broom)


# Forest plot showing distribution of effect sizes across studies
forestplot_simple <- dat[with(dat,order(Fishers_Z_r)),] %>% 
  meta_analysis(yi = Fishers_Z_r, sei = variance, slab = StudyRef)

forest_plot(forestplot_simple)


# Forest plot showing difference in the distribution of effect sizes between groups

forestplot_groupcomparison <- dat[with(dat,order(Fishers_Z_r)),] %>% 
  group_by(Rank_Cont_Cat) %>% 
  meta_analysis(yi = Fishers_Z_r, sei = variance, slab = StudyRef)

fp <- forestplot_groupcomparison %>% 
  forest_plot(group = Rank_Cont_Cat)

fp



# 5) Simulation to show how effect sizes change when either 
# cooperative breeders are analysed with an approach assuming a linear hierarchy
# or in plural breeders females are assigned into two classes of dominants and subordinates


#Simulate 10,000 times the reproductive success of 12 individuals (two groups of six individuals)
#Individuals are in groups of singular cooperative breeders, single dominant has most of the offspring
#Reproductive success is determined by chance, on average top ranking individual has 4.5, all others 2 offspring

cooperativebreeder_effectsizes<-matrix(nrow=10000,ncol=3)
colnames(output)<-c("lm","t_singular","t_groups")

for (i in 1:10000){
  ind1<-round(rnorm(1,mean=2,sd=1),0)
  ind2<-round(rnorm(1,mean=2,sd=1),0)
  ind3<-round(rnorm(1,mean=2,sd=1),0)
  ind4<-round(rnorm(1,mean=2,sd=1),0)
  ind5<-round(rnorm(1,mean=2,sd=1),0)
  ind6<-round(rnorm(1,mean=4.5,sd=1),0)
  ind1a<-round(rnorm(1,mean=2,sd=1),0)
  ind2a<-round(rnorm(1,mean=2,sd=1),0)
  ind3a<-round(rnorm(1,mean=2,sd=1),0)
  ind4a<-round(rnorm(1,mean=2,sd=1),0)
  ind5a<-round(rnorm(1,mean=2,sd=1),0)
  ind6a<-round(rnorm(1,mean=4.5,sd=1),0)
  success<-c(ind1,ind2,ind3,ind4,ind5,ind6,ind1a,ind2a,ind3a,ind4a,ind5a,ind6a)
  rank<-c(1:6,1:6)
  correlation<-lm(success~rank)
  cooperativebreeder_effectsizes[i,1]<-sqrt(summary(correlation)$r.squared)
  dominant<-c(ind6,ind6a)
  subordinates<-c(ind2,ind3,ind4,ind5,ind1,ind2a,ind3a,ind4a,ind5a,ind1a)
  cooperativebreeder_effectsizes[i,2]<-sqrt(t.test(dominant,subordinates)$statistic^2/(t.test(dominant,subordinates)$statistic^2+10))
  dominanthalf<-c(ind6,ind5,ind4,ind6a,ind5a,ind4a)
  subordinatehalf<-c(ind1,ind2,ind3,ind1a,ind2a,ind3a)
  cooperativebreeder_effectsizes[i,3]<-sqrt(t.test(dominanthalf,subordinatehalf)$statistic^2/(t.test(dominanthalf,subordinatehalf)$statistic^2+10))
}

#negative values mean that the t-test has higher effect sizes than the linear regression
hist(cooperativebreeder_effectsizes[,1]-cooperativebreeder_effectsizes[,2])
median(cooperativebreeder_effectsizes[,1]-cooperativebreeder_effectsizes[,2])


#Simulate 10,000 times the reproductive success of 12 individuals (two groups of six individuals)
#Individuals are in a linear hierarchy: going up one position means on average 0.5 more offspring
#Reproductive success is determined by chance, on average top ranking individual has 4.5, bottom 2 offspring

linearhierarchy_effectsizes<-matrix(nrow=10000,ncol=3)
colnames(output)<-c("lm","t_singular","t_groups")

for (i in 1:10000){
ind1<-round(rnorm(1,mean=2,sd=1),0)
ind2<-round(rnorm(1,mean=2.5,sd=1),0)
ind3<-round(rnorm(1,mean=3,sd=1),0)
ind4<-round(rnorm(1,mean=3.5,sd=1),0)
ind5<-round(rnorm(1,mean=4,sd=1),0)
ind6<-round(rnorm(1,mean=4.5,sd=1),0)
ind1a<-round(rnorm(1,mean=2,sd=1),0)
ind2a<-round(rnorm(1,mean=2.5,sd=1),0)
ind3a<-round(rnorm(1,mean=3,sd=1),0)
ind4a<-round(rnorm(1,mean=3.5,sd=1),0)
ind5a<-round(rnorm(1,mean=4,sd=1),0)
ind6a<-round(rnorm(1,mean=4.5,sd=1),0)
success<-c(ind1,ind2,ind3,ind4,ind5,ind6,ind1a,ind2a,ind3a,ind4a,ind5a,ind6a)
rank<-c(1:6,1:6)
correlation<-lm(success~rank)
linearhierarchy_effectsizes[i,1]<-sqrt(summary(correlation)$r.squared)
dominant<-c(ind6,ind6a)
subordinates<-c(ind2,ind3,ind4,ind5,ind1,ind2a,ind3a,ind4a,ind5a,ind1a)
linearhierarchy_effectsizes[i,2]<-sqrt(t.test(dominant,subordinates)$statistic^2/(t.test(dominant,subordinates)$statistic^2+10))
dominanthalf<-c(ind6,ind5,ind4,ind6a,ind5a,ind4a)
subordinatehalf<-c(ind1,ind2,ind3,ind1a,ind2a,ind3a)
linearhierarchy_effectsizes[i,3]<-sqrt(t.test(dominanthalf,subordinatehalf)$statistic^2/(t.test(dominanthalf,subordinatehalf)$statistic^2+10))
}

#positive values mean that the linear regression has higher effect sizes than the t-test
hist(linearhierarchy_effectsizes[,1]-linearhierarchy_effectsizes[,2])
median(linearhierarchy_effectsizes[,1]-linearhierarchy_effectsizes[,2])



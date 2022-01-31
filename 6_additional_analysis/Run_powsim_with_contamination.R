#!/usr/bin/env Rscript

library(optparse)
library(tidyverse)
library(stringr)
library(powsimR)

## contamination functions
source('/data/share/htp/prime-seq_Paper/AmbientRNA/scripts/helper_functions.R')

packageVersion("powsimR")

option_list <- list(
  make_option(c("--input_mat"), type="character", default=NULL, 
              help="Input count matrix", metavar="character"),
  make_option(c("--cores"), type="integer", default=1, 
              help="cores to use", metavar="integer"),
  make_option(c("--nsims"), type="integer", default=10, 
              help="number of simulations", metavar="integer"),
  make_option(c("--p_DE"), type="numeric", default=0.1, 
              help="fraction of genes to be simulted as DE", metavar="numeric"),
  make_option(c("--n_sample"), type="character", default="6,8,10,12,14", 
              help="vector specifying the number of samples to simulate", metavar="character"),
  make_option(c("--ngenes"), type="integer", default="20000", 
              help="vector specifying the number of samples to simulate", metavar="integer"),
  make_option(c("--delta_val"), type="numeric", default=1, 
              help="LFC cutoff for genes to be deemed biologically important", metavar="numeric"),
  make_option(c("--outdir"), type="character", default=NULL, 
              help="directory to which the output should be written", metavar="character"),
  make_option(c("--ambprop"), type="numeric", default=NULL, 
              help="Contamination fraction", metavar="numeric")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
opt$n_sample<-as.integer(str_split(opt$n_sample,pattern = ",",simplify = T))

print("These are the input parameters:")
opt

#opt<-list(input_mat="/data/share/htp/Foxp2/FoxP2_Atlas/Novogene/Analysis/Foxp2_Atlas/powsim/SMI/SMI.KO.Rds",cores=10,nsims=20,p_DE=0.1,n_sample=c(4,6,8,10,12),delta_val=1,outbase="SMI_KI",outdir="/data/share/htp/Foxp2/FoxP2_Atlas/Novogene/Analysis/Foxp2_Atlas/powsim/SMI/")

print(paste("Running powsim with a contamination fraction of:",opt$ambprop))

input_mat<-readRDS(opt$input_mat)

## add contamination
if(opt$ambprop>0){
ObsAmbCounts <-get_ambient_sample(countData = input_mat,
                                  AmbProp = opt$ambprop)

#Add contaminating counts to endogenous counts --> mixture

input_mat <- input_mat + ObsAmbCounts

}

print("Estimating parameters")

estparam<- estimateParam(countData = input_mat,
                         Normalisation = "MR",
                         RNAseq = "bulk", 
                         Protocol = 'UMI',
                         Distribution = 'NB',
                         GeneFilter = 0.1, 
                         SampleFilter = 5,
                         sigma = 1.96, 
                         NCores = opt$cores, 
                         verbose = TRUE)


print(paste0("Finished Parameter estimation at ",Sys.time()))

p.lfc <- function(x) sample(c(-1.5,1.5), size=x,replace=T)*rgamma(x, shape = 1, rate = 2)


setupres <- Setup(nsims = opt$nsims,
                  ngenes = opt$ngenes,
                  p.DE = opt$p_DE,
                  pLFC = p.lfc,
                  n1 = opt$n_sample, 
                  n2 = opt$n_sample,
                  Thinning = NULL, 
                  LibSize = 'equal',
                  estParamRes = estparam,
                  estSpikeRes = NULL,
                  DropGenes = FALSE,
                  setup.seed = 5299,
                  verbose = FALSE)


start<-Sys.time()
print(paste0("Starting Simulations at ",start))

simres <- simulateDE(SetupRes = setupres,
                     Prefilter = NULL,
                     Imputation = NULL,
                     Normalisation = 'MR',
                     DEmethod = "limma-voom", 
                     DEFilter = FALSE,
                     NCores = opt$cores,
                     verbose = FALSE)

warnings()

print(paste0("Finished Simulations at ",Sys.time()))

print(Sys.time()-start)

print(paste0("Starting Evaluation"))

evalderes<- evaluateDE(simRes = simres,
                       alpha.type = 'adjusted',
                       MTC = 'BH',
                       alpha.nominal = 0.05,
                       stratify.by = 'mean',
                       filter.by = 'mean',
                       strata.filtered = 1,
                       target.by = "lfc",
                       delta =opt$delta_val,
                       Table=F)

print(paste0("Finished Evaluation at ",Sys.time()))

saveRDS(estparam,file = paste0(opt$outdir,opt$ambprop,".powsim.estparam.Rds"))
saveRDS(setupres,file = paste0(opt$outdir,opt$ambprop,".powsim.setupres.Rds"))
saveRDS(simres,file = paste0(opt$outdir,opt$ambprop,".powsim.simres.Rds"))
saveRDS(evalderes,file = paste0(opt$outdir,opt$ambprop,".powsim.evalderes.Rds"))

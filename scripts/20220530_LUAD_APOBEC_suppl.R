source('Z:/projects/codes/mg_base.R')

load('origin_datas/TCGA/luad.tcga.exp.RData')
load('origin_datas/GEO/gse31210.exp.RData')
load('origin_datas/GEO/gse72094.exp.RData')
load('origin_datas/GEO/gse50081.exp.RData')

library(GSVA)
library(GSVAdata)
geneSets <- getGmt("./origin_datas/h.all.v7.5.1.symbols.gmt")

tcga.exp.h.all <- gsva(as.matrix(luad.tcga.exp), geneSets,min.sz=10,method='ssgsea')
tcga.exp.h.all <- data.frame(tcga.exp.h.all,check.names = F)

gse31210.exp.h.all <- gsva(as.matrix(gse31210.exp), geneSets,min.sz=10,method='ssgsea')
gse31210.exp.h.all <- data.frame(gse31210.exp.h.all)

gse72094.exp.h.all <- gsva(as.matrix(gse72094.exp), geneSets,min.sz=10,method='ssgsea')
gse72094.exp.h.all <- data.frame(gse72094.exp.h.all)

gse50081.exp.h.all <- gsva(as.matrix(gse50081.exp), geneSets,min.sz=10,method='ssgsea')
gse50081.exp.h.all <- data.frame(gse50081.exp.h.all)

dir.create('origin_datas/GSVA/')
save(tcga.exp.h.all,file = 'origin_datas/GSVA/tcga.exp.h.all.RData')
save(gse31210.exp.h.all,file = 'origin_datas/GSVA/gse31210.exp.h.all.RData')
save(gse72094.exp.h.all,file = 'origin_datas/GSVA/gse72094.exp.h.all.RData')
save(gse50081.exp.h.all,file = 'origin_datas/GSVA/gse50081.exp.h.all.RData')

############
library('GSVA')
library(GSEABase)
gmtFile='origin_datas/c2.cp.kegg.v7.5.1.symbols.gmt'
c2KEGG <- getGmt(gmtFile,
                 collectionType=BroadCollection(category="c2"),
                 geneIdType=SymbolIdentifier())

tcga.kegg.ssgsea <- gsva(as.matrix(luad.tcga.exp), 
                         c2KEGG,
                         method = 'ssgsea',
                         min.sz = 10,
                         max.sz = 500,
                         verbose = TRUE)
save(tcga.kegg.ssgsea,file = 'origin_datas/GSVA/tcga.kegg.ssgsea.RData')

##################
########
immu_ssgsea=function(exp,isTCGA=T){
  
  imm.genes=read.csv(paste0(baseFolder,'/26_immu_signature_pmid_27855702.txt'),sep = '\t',stringsAsFactors = F)
  all.list=list()
  for(s in unique(imm.genes$SetName)){
    inds=which(imm.genes$SetName==s)
    gs=GSEABase::GeneSet(setName=s, setIdentifier=paste0(s,"101")
                         ,geneIds=imm.genes$Gene[inds],GSEABase::SymbolIdentifier()) 
    all.list <- c(all.list, list(gs))
  }
  gsc <- GSEABase::GeneSetCollection(all.list)
  fl <- tempfile()
  GSEABase::toGmt(gsc, fl)
  c2immue=GSEABase::getGmt(fl)
  ssGSEA.immue <- GSVA::gsva(as.matrix(exp), c2immue,method='ssgsea',
                             min.sz=1, max.sz=Inf, verbose=TRUE)
  ssGSEA.immue=t(ssGSEA.immue)
  if(isTCGA){
    rnames=gsub('\\.','-',row.names(ssGSEA.immue))
    row.names(ssGSEA.immue)=rnames
  }
  return(ssGSEA.immue)
}

load('origin_datas/TCGA/luad.tcga.exp.RData')
load('origin_datas/GEO/gse31210.exp.RData')
load('origin_datas/GEO/gse72094.exp.RData')
load('origin_datas/GEO/gse50081.exp.RData')

#####
tcga.exp.immu.ssgsea=immu_ssgsea(luad.tcga.exp)
gse31210.exp.immu.ssgsea=immu_ssgsea(gse31210.exp)
gse72094.exp.immu.ssgsea=immu_ssgsea(gse72094.exp)
gse50081.exp.immu.ssgsea=immu_ssgsea(gse50081.exp)

save(tcga.exp.immu.ssgsea,file='origin_datas/immune/tcga.exp.immu.ssgsea.RData')
save(gse31210.exp.immu.ssgsea,file='origin_datas/immune/gse31210.exp.immu.ssgsea.RData')
save(gse72094.exp.immu.ssgsea,file='origin_datas/immune/gse72094.exp.immu.ssgsea.RData')
save(gse50081.exp.immu.ssgsea,file='origin_datas/immune/gse50081.exp.immu.ssgsea.RData')

library(IOBR)
library(EPIC)
library(estimate)

### CIBERSORT
tcga.exp.cibersort<-deconvo_cibersort(eset=luad.tcga.exp,arrays=T)
gse31210.exp.cibersort<-deconvo_cibersort(eset=gse31210.exp,arrays=T)
gse72094.exp.cibersort<-deconvo_cibersort(eset=gse72094.exp,arrays=T)
gse50081.exp.cibersort<-deconvo_cibersort(eset=gse50081.exp,arrays=T)

save(tcga.exp.cibersort,file='origin_datas/immune/tcga.exp.cibersort.RData')
save(gse31210.exp.cibersort,file='origin_datas/immune/gse31210.exp.cibersort.RData')
save(gse72094.exp.cibersort,file='origin_datas/immune/gse72094.exp.cibersort.RData')
save(gse50081.exp.cibersort,file='origin_datas/immune/gse50081.exp.cibersort.RData')

#### ESTIMATE
tcga.exp.estimate<-deconvo_estimate(eset=luad.tcga.exp)
gse31210.exp.estimate<-deconvo_estimate(eset=gse31210.exp)
gse72094.exp.estimate<-deconvo_estimate(eset=gse72094.exp)
gse50081.exp.estimate<-deconvo_estimate(eset=gse50081.exp)

save(tcga.exp.estimate,file='origin_datas/immune/tcga.exp.estimate.RData')
save(gse31210.exp.estimate,file='origin_datas/immune/gse31210.exp.estimate.RData')
save(gse72094.exp.estimate,file='origin_datas/immune/gse72094.exp.estimate.RData')
save(gse50081.exp.estimate,file='origin_datas/immune/gse50081.exp.estimate.RData')

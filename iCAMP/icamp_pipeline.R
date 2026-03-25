# the script is mainly based on the script in github website (https://github.com/DaliangNing/iCAMP1).

rm(list = ls())
### step1:input and packages
setwd("path/to/your/project") # Run this script from the project root directory
library(iCAMP)
library(ape)

t0=Sys.time() 

### step2:input tables
MAG = read.table("all_MAG_covem.csv", header = TRUE, sep = ",", row.names = 1, check.names = FALSE)
MAG <- MAG[, colSums(MAG != 0) > 0]
treat=read.table("mapping.csv",header=T,sep=",",row.names=1) #treatment
tre.file="rooted_tree.tree" #tree file
clas=read.table("anonotation.csv", header = TRUE, sep=",", row.names = 1,
                as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                check.names = FALSE) #classification

i_treatment<-as.numeric(1)#specify the column used to select group
prefix<-"mag_tpm_"#this is a label used for each job in HPC
rand.time<-as.numeric(200)# set random time for iCAMP
nworker<-as.numeric(8) #number of CPUs in your server nworker<-as.numeric(32)
memory.G <- as.numeric(128)#memory (G) used in server memory.G <- as.numeric(900)
bin.size.limit<-as.numeric(12) #bin size for iCAMP, usually 12-48
sig.index<-"Confidence"#the index for null model significance test
taxo.metric<-"bray"#taxonomic beta diversity index
detail.null<-as.logical("TRUE")  #detail.null # this is very important          
is.logical(detail.null)
save.wd<-"path/to/your/project/icamp"
if(!dir.exists(save.wd)){dir.create(save.wd)}
MAG <- round(MAG, digits = 0)
phylo.metric<-"bMPD"#character, the metric for phylogenetic null model analysis"

tree=read.tree(file = tre.file)
comm<-t(MAG)
sampid.check=match.name(rn.list=list(comm=comm,treat=treat))
treat=sampid.check$treat
comm=sampid.check$comm
comm=comm[,colSums(comm)>0,drop=FALSE]  

###step3:match OTU IDs in OTU table and tree file
spid.check=match.name(cn.list=list(comm=comm),rn.list=list(clas=clas),tree.list=list(tree=tree))

# for the example data, the output should be "All match very well".
# for your data files, if you have not matched the IDs before, the unmatched OTUs will be removed
comm=spid.check$comm
clas=spid.check$clas
tree=spid.check$tree

tree$root.edge <- 0
is.rooted(tree)#true 

library(iCAMP)
print("start calculate pd.big")
setwd(save.wd)

###step4:running iCAMP
if(!file.exists("pd.desc")) {
  pd.big=iCAMP::pdist.big(tree = tree, wd=save.wd, nworker = nworker, memory.G = memory.G,treepath.file=paste0(prefix,"path.rda"),pd.spname.file=paste0(prefix,"pd.taxon.name.csv"),pd.backingfile=paste0(prefix,"pd.bin"),pd.desc.file=paste0(prefix,"pd.desc"))
  
}else{
  
  pd.big=list()
  pd.big$tip.label=read.csv(paste0(save.wd,"/pd.taxon.name",prefix,".csv"),row.names = 1,stringsAsFactors = FALSE)[,1]
  pd.big$pd.wd=save.wd
  pd.big$pd.file=paste0(prefix,".pd.desc")
  pd.big$pd.name.file=paste0(prefix,".pd.taxon.name.csv")
}

print("start calculate iCAMP.big")

icres=iCAMP::icamp.big(comm=comm, pd.desc = pd.big$pd.file, pd.spname=pd.big$tip.label,
                       pd.wd = pd.big$pd.wd, rand = rand.time, tree=tree,
                       prefix = prefix, ds = 0.2, pd.cut = NA, sp.check = TRUE,
                       phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                       phylo.metric = "bMPD", sig.index=sig.index, bin.size.limit = bin.size.limit, 
                       nworker = nworker, memory.G = memory.G, rtree.save = FALSE, detail.save = TRUE, 
                       qp.save = FALSE, detail.null = FALSE, ignore.zero = TRUE, output.wd = save.wd, 
                       correct.special = TRUE, unit.sum = rowSums(comm), special.method = "depend",
                       ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = "no",meta.ab = NULL)
save(icres,file = paste0(prefix,".iCAMP.save_important.rda")) 

# load("D:/2025Workshop/05icamp/output/mag_tpm_.iCAMP.save_important.rda")

i=as.numeric(1)
treat.use=treat[,i,drop=FALSE]

###step5:iCAMP bin level statistics
icbin=icamp.bins(icamp.detail = icres$detail,treat = treat,
                 clas=clas,silent=FALSE, boot = TRUE,
                 rand.time = rand.time,between.group = TRUE)
save(icbin,file = paste0(prefix,".iCAMP.Summary.rda")) # just to archive the result. rda file is automatically compressed, and easy to load into R.
write.csv(icbin$Pt,file = paste0(prefix,".ProcessImportance_EachGroup.csv"),row.names = FALSE)
write.csv(icbin$Ptk,file = paste0(prefix,".ProcessImportance_EachBin_EachGroup.csv"),row.names = FALSE)
write.csv(icbin$Ptuv,file = paste0(prefix,".ProcessImportance_EachTurnover.csv"),row.names = FALSE)
write.csv(icbin$BPtk,file = paste0(prefix,".BinContributeToProcess_EachGroup.csv"),row.names = FALSE)

write.csv(data.frame(ID=rownames(icbin$Class.Bin),icbin$Class.Bin,stringsAsFactors = FALSE),
          file = paste0(prefix,".Taxon_Bin.csv"),row.names = FALSE)
write.csv(icbin$Bin.TopClass,file = paste0(prefix,".Bin_TopTaxon.csv"),row.names = FALSE)

###step6  Bootstrapping test
# please specify which column in the treatment information table.
i=as.numeric(1)
treat.use=treat[,i,drop=FALSE]
icamp.result=icres$CbMPDiCBraya
icboot=iCAMP::icamp.boot(icamp.result = icamp.result,treat = treat.use,rand.time = rand.time,
                         compare = TRUE,silent = FALSE,between.group = TRUE,ST.estimation = TRUE)
save(icboot,file=paste0(prefix,".iCAMP.Boot.",colnames(treat)[i],".rda"))
write.csv(icboot$summary,file = paste0(prefix,".iCAMP.BootSummary.",colnames(treat)[i],".csv"),row.names = FALSE)
write.csv(icboot$compare,file = paste0(prefix,".iCAMP.Compare.",colnames(treat)[i],".csv"),row.names = FALSE)
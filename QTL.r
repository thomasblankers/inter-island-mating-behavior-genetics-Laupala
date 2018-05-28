### merge linkage maps from HC1 and HC2 ####
## hc1_map_geno and hc2_map_geno are the genotypes for each marker in the linkage maps for the HC1 and HC2 cross families (see supplemetary info)
CerEuk_LG1_LPmerged<-LPmerge(list(hc1_map_geno[hc1_map_geno["LG"]==1,c("locus","cM")],hc2_map_geno[hc2_map_geno["LG"]==1,c("locus","cM")]), max.interval=1:10, weights=c(1.15,1))
CerEuk_LG1_LPmerged_final<-CerEuk_LG1_LPmerged[3]
CerEuk_LG2_LPmerged<-LPmerge(list(hc1_map_geno[hc1_map_geno["LG"]==2,c("locus","cM")],hc2_map_geno[hc2_map_geno["LG"]==2,c("locus","cM")]), max.interval=1:10, weights=c(1,1))
CerEuk_LG2_LPmerged_final<-CerEuk_LG2_LPmerged[5]
CerEuk_LG3_LPmerged<-LPmerge(list(hc1_map_geno[hc1_map_geno["LG"]==3,c("locus","cM")],hc2_map_geno[hc2_map_geno["LG"]==3,c("locus","cM")]), max.interval=1:10, weights=c(1,5))
CerEuk_LG3_LPmerged_final<-CerEuk_LG3_LPmerged[10]
CerEuk_LG4_LPmerged<-LPmerge(list(hc1_map_geno[hc1_map_geno["LG"]==4,c("locus","cM")],hc2_map_geno[hc2_map_geno["LG"]==4,c("locus","cM")]), max.interval=1:10, weights=c(1,1))
CerEuk_LG4_LPmerged_final<-CerEuk_LG4_LPmerged[10]
CerEuk_LG5_LPmerged<-LPmerge(list(hc1_map_geno[hc1_map_geno["LG"]==5,c("locus","cM")],hc2_map_geno[hc2_map_geno["LG"]==5,c("locus","cM")]), max.interval=1:10, weights=c(1.1,1))
CerEuk_LG5_LPmerged_final<-CerEuk_LG5_LPmerged[6]
CerEuk_LG6_LPmerged<-LPmerge(list(hc1_map_geno[hc1_map_geno["LG"]==6,c("locus","cM")],hc2_map_geno[hc2_map_geno["LG"]==6,c("locus","cM")]), max.interval=1:10, weights=c(1,1))
CerEuk_LG6_LPmerged_final<-CerEuk_LG6_LPmerged[1]
CerEuk_LG7_LPmerged<-LPmerge(list(hc1_map_geno[hc1_map_geno["LG"]==7,c("locus","cM")],hc2_map_geno[hc2_map_geno["LG"]==7,c("locus","cM")]), max.interval=1:10, weights=c(1,1.2))
CerEuk_LG7_LPmerged_final<-CerEuk_LG7_LPmerged[10]
CerEuk_LGX_LPmerged<-LPmerge(list(hc1_map_geno[hc1_map_geno["LG"]=="X",c("locus","cM")],hc2_map_geno[hc2_map_geno["LG"]=="X",c("locus","cM")]), max.interval=1:10, weights=c(1,1))
CerEuk_LGX_LPmerged_final<-CerEuk_LGX_LPmerged[2]

CerEuk_LPmerged_map<-rbind(CerEuk_LG1_LPmerged_final[[1]],
CerEuk_LG2_LPmerged_final[[1]],
CerEuk_LG3_LPmerged_final[[1]],
CerEuk_LG4_LPmerged_final[[1]],
CerEuk_LG5_LPmerged_final[[1]],
CerEuk_LG6_LPmerged_final[[1]],
CerEuk_LG7_LPmerged_final[[1]],
CerEuk_LGX_LPmerged_final[[1]])

########################################################### r/QTL ######################################################################


########################################################### HC1 ######################################################################

CerEukHC1<-read.cross(format="csv", file="CerEuk_HC1_CrossFile_NEW.csv", estimate.map=FALSE)
CerEukHC1_jitter<-jittermap(CerEukHC1)
CerEukHC1_jitter<-sim.geno(CerEukHC1_jitter, step=1, n.draws=1000, error.prob=0.001)
CerEukHC1_jitter<-calc.genoprob(CerEukHC1_jitter, step=1, error.prob=0.001)

## scanone

qtl_CerEukHC1_imp<-scanone(CerEukHC1_jitter,method="imp", pheno.col="pr")
qtl_CerEukHC1_hk<-scanone(CerEukHC1_jitter,method="hk", pheno.col="pr")

qtl_CerEukHC1_imp_1000perm<-scanone(CerEukHC1_jitter,method="imp", pheno.col="pr", n.perm=1000, perm.Xsp=TRUE) 
qtl_CerEukHC1_hk_1000perm<-scanone(CerEukHC1_jitter,method="hk", pheno.col="pr", n.perm=1000, perm.Xsp=TRUE)



## cim

cim_CerEukHC1_imp<-cim(CerEukHC1_jitter, pheno.col="pr", method="imp", n.marcovar=3, window=20)
cim_CerEukHC1_hk<-cim(CerEukHC1_jitter, pheno.col="pr", method="hk", n.marcovar=3, window=20)

cim_CerEukHC1_hk_2<-cim(CerEukHC1_jitter, pheno.col="pr", method="hk", n.marcovar=2, window=20)

cim_CerEukHC1_hk_2_1000perm<-cim(CerEukHC1_jitter, pheno.col="pr", method="hk", n.marcovar=2, window=20, n.perm=1000)

# plot CIM
cairo_pdf("CerEukHC1_CIM_2_20.pdf", width=8, height=3.5)
plot(cim_CerEukHC1_hk_2);abline(h=summary(cim_CerEukHC1_hk_2_1000perm)[1,1])
dev.off()


qtl_CerEukHC1<-c(1,2,3,4,5,6,7,"X")

## scantwo
twoqtl_CerEukHC1_imp<-scantwo(CerEukHC1_jitter, pheno.col="pr", method="imp", chr=qtl_CerEukHC1, clean.output=TRUE)
twoqtl_CerEukHC1_hk<-scantwo(CerEukHC1_jitter, pheno.col="pr", method="hk", chr=qtl_CerEukHC1, clean.output=TRUE)

# check for potential interactions:

plot(twoqtl_CerEukHC1_imp, lower="cond-int")
plot(twoqtl_CerEukHC1_hk, lower="cond-int")
# check for secondary peaks at specific chromosomes:

plot(twoqtl_CerEukHC1_imp, chr=1, lower="cond-int", upper="cond-add")
plot(twoqtl_CerEukHC1_hk, chr=2, lower="cond-int", upper="cond-add")

## mutiple QTL mapping

# permutate scantwo to obtain penalized LOD thresholds
twoqtl_CerEukHC1_hk_1000perm<-scantwo(CerEukHC1_jitter, pheno.col="pr", method="hk", chr=qtl_CerEukHC1, clean.output=TRUE, n.perm=1000, perm.Xsp=FALSE) # takes about 40 minutes
penalties=calc.penalties(twoqtl_CerEukHC1_hk_1000perm, alpha=0.05)

    main    heavy    light 
3.446001 5.603086 3.434136

# build MQM
CerEukHC1_init_qtl<-makeqtl(CerEukHC1_jitter,chr= 2, pos= 115.5)
summary(fitqtl(CerEukHC1_jitter, qtl=CerEukHC1_init_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula= y ~ Q1))
CerEukHC1_ref_qtl<-refineqtl(CerEukHC1_jitter, pheno.col="pr", method="imp", qtl = CerEukHC1_init_qtl)
summary(fitqtl(CerEukHC1_jitter, qtl=CerEukHC1_ref_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula= y ~ Q1))
CerEukHC1_add_init_qtl <- addqtl(CerEukHC1_jitter, qtl=CerEukHC1_init_qtl, pheno.col="pr", method="imp", formula= y ~ Q1)

CerEukHC1_expand_qtl<-addtoqtl(CerEukHC1_jitter,CerEukHC1_init_qtl, 1, 72.0)
summary(fitqtl(CerEukHC1_jitter, qtl=CerEukHC1_expand_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula= y ~ Q1 + Q2))
CerEukHC1_expand_qtl<-refineqtl(CerEukHC1_jitter, pheno.col="pr", method="imp", qtl = CerEukHC1_expand_qtl)
summary(fitqtl(CerEukHC1_jitter, qtl=CerEukHC1_expand_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula= y ~ Q1 + Q2))
CerEukHC1_add_expand_qtl <- addqtl(CerEukHC1_jitter, qtl=CerEukHC1_expand_qtl, pheno.col="pr", method="imp", formula= y ~ Q1 + Q2)
CerEukHC1_expand_qtl_addint<-addint(CerEukHC1_jitter, pheno.col="pr",  qtl=CerEukHC1_expand_qtl, method="imp")

# plot MQM
cairo_pdf("CerEukHC1_fitqtl.pdf", width=8, height=3.5)
plotLodProfile(CerEukHC1_expand2_qtl,showallchr=TRUE); abline(h=3.446001)
dev.off()

# obtain 1-LOD and 1.5-LOD intervals
rbind(
lodint(CerEukHC1_expand2_qtl, drop=1, qtl.index=1),
lodint(CerEukHC1_expand2_qtl, drop=1.5 qtl.index=1),
lodint(CerEukHC1_expand2_qtl, drop=1, qtl.index=2),
lodint(CerEukHC1_expand2_qtl, drop=1.5, qtl.index=2))

# calculate true number of loci using Otto & Jones (2000) formula;  see QTL_power_detect.r for details
true_loci_HC1<-otto_jones_ci(D=0.8381,M=0.168,nd=2,alpha=0.05,amin=0.116,res=4,res2=2,max.loci=100)

########################################################### HC2 ######################################################################

CerEukHC2<-read.cross(format="csv", file="CerEuk_HC2_CrossFile_NEW.csv", estimate.map=FALSE)
CerEukHC2_jitter<-jittermap(CerEukHC2)
CerEukHC2_jitter<-sim.geno(CerEukHC2_jitter, step=1, n.draws=1000, error.prob=0.001)
CerEukHC2_jitter<-calc.genoprob(CerEukHC2_jitter, step=1, error.prob=0.001)

## scanone

qtl_CerEukHC2_imp<-scanone(CerEukHC2_jitter,method="imp", pheno.col="pr")
qtl_CerEukHC2_hk<-scanone(CerEukHC2_jitter,method="hk", pheno.col="pr")

qtl_CerEukHC2_imp_1000perm<-scanone(CerEukHC2_jitter,method="imp", pheno.col="pr", n.perm=1000, perm.Xsp=TRUE)
qtl_CerEukHC2_hk_1000perm<-scanone(CerEukHC2_jitter,method="hk", pheno.col="pr", n.perm=1000, perm.Xsp=TRUE)


## cim


cim_CerEukHC2_hk_3<-cim(CerEukHC2_jitter, pheno.col="pr", method="hk", n.marcovar=3, window=20)

cim_CerEukHC2_hk_3_1000perm<-cim(CerEukHC2_jitter, pheno.col="pr", method="hk", n.marcovar=3, window=20, n.perm=1000)

# create marker covariates from cim QTL
qtl1 <- pull.geno(fill.geno(CerEukHC2_jitter))[,find.marker(CerEukHC2_jitter,chr=1,pos=63)] 
qtl2 <- pull.geno(fill.geno(CerEukHC2_jitter))[,find.marker(CerEukHC2_jitter,chr=2,pos=59)]
qtl3 <- pull.geno(fill.geno(CerEukHC2_jitter))[,find.marker(CerEukHC2_jitter,chr=5,pos=38.4)]


# do X specific permutations for CIM
cim_CerEukHC2_hk_1000perm_Xsp<-scanone(CerEukHC2_jitter,method="hk", pheno.col="pr", addcovar=cbind(qtl1,qtl2,qtl3), n.perm=1000, perm.Xsp=TRUE)

cairo_pdf("CerEukHC2_CIM_3_20.pdf", width=8, height=3.5)
plot(cim_CerEukHC2_hk_3);abline(h=summary(cim_CerEukHC2_hk_3_1000perm)[1,1]); abline(h=summary(cim_CerEukHC2_hk_1000perm_Xsp)[[2]][1,1], lty="dotted")
dev.off()


qtl_CerEukHC2<-c(1,2,3,4,5,6,7,"X")

## scantwo
twoqtl_CerEukHC2_imp<-scantwo(CerEukHC2_jitter, pheno.col="pr", method="imp", chr=qtl_CerEukHC2, clean.output=TRUE)
twoqtl_CerEukHC2_hk<-scantwo(CerEukHC2_jitter, pheno.col="pr", method="hk", chr=qtl_CerEukHC2, clean.output=TRUE)

# check for potential interactions:

plot(twoqtl_CerEukHC2_imp, lower="cond-int")
plot(twoqtl_CerEukHC2_hk, lower="cond-int")

# check for secondary peaks at specific chromosomes:
	
plot(twoqtl_CerEukHC2_imp, chr=1, lower="cond-int", upper="cond-add")
plot(twoqtl_CerEukHC2_imp, chr=2, lower="cond-int", upper="cond-add")
plot(twoqtl_CerEukHC2_imp, chr=5, lower="cond-int", upper="cond-add")
plot(twoqtl_CerEukHC2_imp, chr="X", lower="cond-int", upper="cond-add")

## MQM

# permuate scantwo
twoqtl_CerEukHC2_hk_1000perm<-scantwo(CerEukHC2_jitter, pheno.col="pr", method="hk", chr=qtl_CerEukHC2, clean.output=TRUE, n.perm=1000, perm.Xsp=FALSE) # takes about 40 minutes
penalties=calc.penalties(twoqtl_CerEukHC2_hk_1000perm, alpha=0.05)

    main    heavy    light 
3.350607 5.530137 3.396078

# build MQM
CerEukHC2_init_qtl<-makeqtl(CerEukHC2_jitter,chr= 2, pos= 59)
summary(fitqtl(CerEukHC2_jitter, qtl=CerEukHC2_init_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1))
CerEukHC2_init_qtl<-refineqtl(CerEukHC2_jitter, pheno.col="pr", method="imp", qtl = CerEukHC2_init_qtl)
summary(fitqtl(CerEukHC2_jitter, qtl=CerEukHC2_init_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula = y ~ Q1))
CerEukHC2_add_init_qtl <- addqtl(CerEukHC2_jitter, qtl=CerEukHC2_init_qtl, pheno.col="pr", method="imp", formula= y ~ Q1)

CerEukHC2_expand_qtl<-addtoqtl(CerEukHC2_jitter,CerEukHC2_init_qtl, 1, 65)
summary(fitqtl(CerEukHC2_jitter, qtl=CerEukHC2_expand_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula= y ~ Q1 + Q2)) 
CerEukHC2_expand_qtl<-refineqtl(CerEukHC2_jitter, pheno.col="pr", method="imp", qtl = CerEukHC2_expand_qtl)
summary(fitqtl(CerEukHC2_jitter, qtl=CerEukHC2_expand_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula= y ~ Q1 + Q2))
CerEukHC2_add_expand_qtl <- addqtl(CerEukHC2_jitter, qtl=CerEukHC2_expand_qtl, pheno.col="pr", method="imp", formula= y ~ Q1 + Q2)
CerEukHC2_expand_qtl_addint<-addint(CerEukHC2_jitter, pheno.col="pr",  qtl=CerEukHC2_expand_qtl, method="imp")

CerEukHC2_expand2_qtl<-addtoqtl(CerEukHC2_jitter,CerEukHC2_expand_qtl, 5, 38.4)
summary(fitqtl(CerEukHC2_jitter, qtl=CerEukHC2_expand2_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula= y ~ Q1 + Q2 + Q3)) 
CerEukHC2_expand2_qtl<-refineqtl(CerEukHC2_jitter, pheno.col="pr", method="imp", qtl = CerEukHC2_expand2_qtl)
summary(fitqtl(CerEukHC2_jitter, qtl=CerEukHC2_expand2_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula= y ~ Q1 + Q2 + Q3))
CerEukHC2_add_expand2_qtl <- addqtl(CerEukHC2_jitter, qtl=CerEukHC2_expand2_qtl, pheno.col="pr", method="imp", formula= y ~ Q1 + Q2 + Q3)
CerEukHC2_expand2_qtl_addint<-addint(CerEukHC2_jitter, pheno.col="pr",  qtl=CerEukHC2_expand2_qtl, method="imp")

CerEukHC2_expand3_qtl<-addtoqtl(CerEukHC2_jitter,CerEukHC2_expand2_qtl, "X", 22)
summary(fitqtl(CerEukHC2_jitter, qtl=CerEukHC2_expand3_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula= y ~ Q1 + Q2 + Q3 + Q4)) 
CerEukHC2_expand3_qtl<-refineqtl(CerEukHC2_jitter, pheno.col="pr", method="imp", qtl = CerEukHC2_expand3_qtl)
summary(fitqtl(CerEukHC2_jitter, qtl=CerEukHC2_expand3_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula= y ~ Q1 + Q2 + Q3 + Q4))
CerEukHC2_add_expand3_qtl <- addqtl(CerEukHC2_jitter, qtl=CerEukHC2_expand3_qtl, pheno.col="pr", method="imp", formula= y ~ Q1 + Q2 + Q3 + Q4)
CerEukHC2_expand3_qtl_addint<-addint(CerEukHC2_jitter, pheno.col="pr",  qtl=CerEukHC2_expand3_qtl, method="imp")

cairo_pdf("CerEukHC2_fitqtl.pdf", width=8, height=3.5)
plotLodProfile(CerEukHC2_expand3_qtl,showallchr=TRUE); abline(h=3.350607)
dev.off()


rbind(
lodint(CerEukHC2_expand3_qtl, drop=1,  qtl.index=1),
lodint(CerEukHC2_expand3_qtl, drop=1.5,  qtl.index=1),
lodint(CerEukHC2_expand3_qtl, drop=1, qtl.index=2),
lodint(CerEukHC2_expand3_qtl, drop=1.5, qtl.index=2),
lodint(CerEukHC2_expand3_qtl, drop=1, qtl.index=3),
lodint(CerEukHC2_expand3_qtl, drop=1.5, qtl.index=3),
lodint(CerEukHC2_expand3_qtl, drop=1, qtl.index=4),
lodint(CerEukHC2_expand3_qtl, drop=1.5, qtl.index=4))

# true number of loci
true_loci_HC2<-otto_jones_ci(D=0.8381,M=0.1586,nd=4,alpha=0.05,amin=0.1026,res=4,res2=2,max.loci=100)


########################################################### LPMerge ######################################################################


CerEuk_LPmerge<-read.cross(format="csv", file="CerEuk_LPmerged_map_reorder_CrossFile_NEW.csv", estimate.map=FALSE)
CerEuk_LPmerge_jitter<-jittermap(CerEuk_LPmerge)
CerEuk_LPmerge_jitter<-sim.geno(CerEuk_LPmerge_jitter, step=1, n.draws=1000, error.prob=0.001)
CerEuk_LPmerge_jitter<-calc.genoprob(CerEuk_LPmerge_jitter, step=1, error.prob=0.001)

## scanone

qtl_CerEuk_LPmerge_imp<-scanone(CerEuk_LPmerge_jitter,method="imp", addcovar=as.numeric(CerEuk_LPmerge_jitter$pheno$cross),pheno.col="pr")
qtl_CerEuk_LPmerge_hk<-scanone(CerEuk_LPmerge_jitter,method="hk", addcovar=as.numeric(CerEuk_LPmerge_jitter$pheno$cross), pheno.col="pr")

qtl_CerEuk_LPmerge_imp_1000perm<-scanone(CerEuk_LPmerge_jitter,method="imp", addcovar=as.numeric(CerEuk_LPmerge_jitter$pheno$cross), pheno.col="pr", n.perm=1000, perm.Xsp=TRUE)
qtl_CerEuk_LPmerge_hk_1000perm<-scanone(CerEuk_LPmerge_jitter,method="hk", addcovar=as.numeric(CerEuk_LPmerge_jitter$pheno$cross), pheno.col="pr", n.perm=1000, perm.Xsp=TRUE)

## cim

cim_CerEuk_LPmerge_imp<-cim(CerEuk_LPmerge_jitter, pheno.col="pr", method="imp", n.marcovar=3, window=10)
cim_CerEuk_LPmerge_hk<-cim(CerEuk_LPmerge_jitter, pheno.col="pr", method="hk", n.marcovar=3, window=30)


cim_CerEuk_LPmerge_hk_3<-cim(CerEuk_LPmerge_jitter, pheno.col="pr", method="hk", n.marcovar=3, window=20)


cim_CerEuk_LPmerge_hk_3_1000perm<-cim(CerEuk_LPmerge_jitter, pheno.col="pr", method="hk", n.marcovar=3, window=20, n.perm=1000)

qtl1 <- pull.geno(fill.geno(CerEuk_LPmerge_jitter))[,find.marker(CerEuk_LPmerge_jitter,chr=2,pos=57.6)] 
qtl2 <- pull.geno(fill.geno(CerEuk_LPmerge_jitter))[,find.marker(CerEuk_LPmerge_jitter,chr=1,pos=68.4)] 
qtl3 <- pull.geno(fill.geno(CerEuk_LPmerge_jitter))[,find.marker(CerEuk_LPmerge_jitter,chr=5,pos=27.8)] 

# do X specific permutations for CIM
cim_CerEuk_LPmerge_hk_1000perm_Xsp<-scanone(CerEuk_LPmerge_jitter,method="hk", pheno.col="pr", addcovar=cbind(qtl1,qtl2,qtl3), n.perm=1000, perm.Xsp=TRUE)

cairo_pdf("CerEukLPmerge_CIM_3_20.pdf", width=8, height=3.5)
plot(cim_CerEuk_LPmerge_hk_3, ylim=c(0,35));abline(h=4.15); abline(h= 2.54, lty="dotted")
dev.off()

qtl_CerEuk_LPmerge<-c(1,2,3,4,5,6,"X")

## scantwo
twoqtl_CerEuk_LPmerge_imp<-scantwo(CerEuk_LPmerge_jitter, pheno.col="pr", method="imp", chr=qtl_CerEuk_LPmerge, clean.output=TRUE)
twoqtl_CerEuk_LPmerge_hk<-scantwo(CerEuk_LPmerge_jitter, pheno.col="pr", method="hk", chr=qtl_CerEuk_LPmerge, clean.output=TRUE)

# check for potential interactions:

plot(twoqtl_CerEuk_LPmerge_imp, lower="cond-int")
plot(twoqtl_CerEuk_LPmerge_hk, lower="cond-int")

# check for secondary peaks at specific chromosomes:

plot(twoqtl_CerEuk_LPmerge_hk, chr=1, lower="cond-int", upper="cond-add")
plot(twoqtl_CerEuk_LPmerge_hk, chr=2, lower="cond-int", upper="cond-add")
plot(twoqtl_CerEuk_LPmerge_hk, chr=3, lower="cond-int", upper="cond-add")
plot(twoqtl_CerEuk_LPmerge_hk, chr=4, lower="cond-int", upper="cond-add")
plot(twoqtl_CerEuk_LPmerge_hk, chr=5, lower="cond-int", upper="cond-add")
plot(twoqtl_CerEuk_LPmerge_hk, chr=6, lower="cond-int", upper="cond-add")
plot(twoqtl_CerEuk_LPmerge_hk, chr="X", lower="cond-int", upper="cond-add")

twoqtl_CerEuk_LPmerge_hk_1000perm<-scantwo(CerEuk_LPmerge_jitter, pheno.col="pr", method="hk", chr=qtl_CerEuk_LPmerge, clean.output=TRUE, n.perm=1000, perm.Xsp=FALSE) # takes about 40 minutes
penalties=calc.penalties(twoqtl_CerEuk_LPmerge_hk_1000perm, alpha=0.05)

    main    heavy    light 
3.154491 5.393225 3.362881


CerEuk_LPmerge_init_qtl<-makeqtl(CerEuk_LPmerge_jitter,chr= 2 , pos= 57.6)
summary(fitqtl(CerEuk_LPmerge_jitter, qtl=CerEuk_LPmerge_init_qtl, dropone=TRUE, get.ests=TRUE, covar=data.frame(cross=CerEuk_LPmerge_jitter$pheno$cross), pheno.col="pr", formula= y ~ Q1 + cross))
CerEuk_LPmerge_init_qtl<-refineqtl(CerEuk_LPmerge_jitter, pheno.col="pr", method="imp", covar=data.frame(cross=CerEuk_LPmerge_jitter$pheno$cross), qtl = CerEuk_LPmerge_init_qtl)
summary(fitqtl(CerEuk_LPmerge_jitter, qtl=CerEuk_LPmerge_init_qtl, dropone=TRUE, get.ests=TRUE, covar=data.frame(cross=CerEuk_LPmerge_jitter$pheno$cross), pheno.col="pr", formula= y ~ Q1 ))
CerEuk_LPmerge_add_init_qtl <- addqtl(CerEuk_LPmerge_jitter, qtl=CerEuk_LPmerge_init_qtl, covar=data.frame(cross=CerEuk_LPmerge_jitter$pheno$cross), pheno.col="pr", method="imp", formula= y ~ Q1 + cross)

CerEuk_LPmerge_expand_qtl<-addtoqtl(CerEuk_LPmerge_jitter, CerEuk_LPmerge_init_qtl, 1, 68.4)
summary(fitqtl(CerEuk_LPmerge_jitter, qtl=CerEuk_LPmerge_expand_qtl, dropone=TRUE, get.ests=TRUE, covar=data.frame(cross=CerEuk_LPmerge_jitter$pheno$cross), pheno.col="pr", formula= y ~ Q1 + Q2 + cross))
CerEuk_LPmerge_expand_qtl<-refineqtl(CerEuk_LPmerge_jitter, pheno.col="pr", method="imp", covar=data.frame(cross=CerEuk_LPmerge_jitter$pheno$cross), qtl = CerEuk_LPmerge_expand_qtl)
summary(fitqtl(CerEuk_LPmerge_jitter, qtl=CerEuk_LPmerge_expand_qtl, dropone=TRUE, get.ests=TRUE, covar=data.frame(cross=CerEuk_LPmerge_jitter$pheno$cross), pheno.col="pr", formula= y ~ Q1 + Q2 + cross))
CerEuk_LPmerge_add_expand_qtl <- addqtl(CerEuk_LPmerge_jitter, qtl=CerEuk_LPmerge_expand_qtl, covar=data.frame(cross=CerEuk_LPmerge_jitter$pheno$cross), pheno.col="pr", method="imp", formula= y ~ Q1 + Q2 + cross)
CerEuk_LPmerge_expand_qtl_addint<-addint(CerEuk_LPmerge_jitter, pheno.col="pr",  qtl=CerEuk_LPmerge_expand_qtl, method="imp")


CerEuk_LPmerge_expand2_qtl<-addtoqtl(CerEuk_LPmerge_jitter, CerEuk_LPmerge_expand_qtl, 5, 29)
summary(fitqtl(CerEuk_LPmerge_jitter, qtl=CerEuk_LPmerge_expand2_qtl, dropone=TRUE, get.ests=TRUE, covar=data.frame(cross=CerEuk_LPmerge_jitter$pheno$cross), pheno.col="pr", formula= y ~ Q1 + Q2 + Q3 + cross))
CerEuk_LPmerge_expand2_qtl<-refineqtl(CerEuk_LPmerge_jitter, pheno.col="pr", method="imp", covar=data.frame(cross=CerEuk_LPmerge_jitter$pheno$cross), qtl = CerEuk_LPmerge_expand2_qtl)
summary(fitqtl(CerEuk_LPmerge_jitter, qtl=CerEuk_LPmerge_expand2_qtl, dropone=TRUE, get.ests=TRUE, covar=data.frame(cross=CerEuk_LPmerge_jitter$pheno$cross), pheno.col="pr", formula= y ~ Q1 + Q2 + Q3 + cross))
CerEuk_LPmerge_add_expand2_qtl <- addqtl(CerEuk_LPmerge_jitter, qtl=CerEuk_LPmerge_expand2_qtl, covar=data.frame(cross=CerEuk_LPmerge_jitter$pheno$cross), pheno.col="pr", method="imp", formula= y ~ Q1 + Q2 + Q3 + cross)
CerEuk_LPmerge_expand2_qtl_addint<-addint(CerEuk_LPmerge_jitter, pheno.col="pr",  qtl=CerEuk_LPmerge_expand2_qtl, method="imp")

CerEuk_LPmerge_expand3_qtl<-addtoqtl(CerEuk_LPmerge_jitter, CerEuk_LPmerge_expand2_qtl, "X", 62)
summary(fitqtl(CerEuk_LPmerge_jitter, qtl=CerEuk_LPmerge_expand3_qtl, dropone=TRUE, get.ests=TRUE, covar=data.frame(cross=CerEuk_LPmerge_jitter$pheno$cross), pheno.col="pr", formula= y ~ Q1 + Q2 + Q3 + Q4 + cross))
CerEuk_LPmerge_expand3_qtl<-refineqtl(CerEuk_LPmerge_jitter, pheno.col="pr", method="imp", covar=data.frame(cross=CerEuk_LPmerge_jitter$pheno$cross), qtl = CerEuk_LPmerge_expand3_qtl)
summary(fitqtl(CerEuk_LPmerge_jitter, qtl=CerEuk_LPmerge_expand3_qtl, dropone=TRUE, get.ests=TRUE, covar=data.frame(cross=CerEuk_LPmerge_jitter$pheno$cross), pheno.col="pr", formula= y ~ Q1 + Q2 + Q3 + Q4 + cross))
CerEuk_LPmerge_add_expand3_qtl <- addqtl(CerEuk_LPmerge_jitter, qtl=CerEuk_LPmerge_expand3_qtl, covar=data.frame(cross=CerEuk_LPmerge_jitter$pheno$cross), pheno.col="pr", method="imp", formula= y ~ Q1 + Q2 + Q3 + Q4 + cross)
CerEuk_LPmerge_expand3_qtl_addint<-addint(CerEuk_LPmerge_jitter, pheno.col="pr",  qtl=CerEuk_LPmerge_expand3_qtl, method="imp")

CerEuk_LPmerge_expand4_qtl<-addtoqtl(CerEuk_LPmerge_jitter, CerEuk_LPmerge_expand3_qtl, 3, 71)
summary(fitqtl(CerEuk_LPmerge_jitter, qtl=CerEuk_LPmerge_expand4_qtl, dropone=TRUE, get.ests=TRUE, covar=data.frame(cross=CerEuk_LPmerge_jitter$pheno$cross), pheno.col="pr", formula= y ~ Q1 + Q2 + Q3 + Q4 + Q5 + cross))
CerEuk_LPmerge_expand4_qtl<-refineqtl(CerEuk_LPmerge_jitter, pheno.col="pr", method="imp", covar=data.frame(cross=CerEuk_LPmerge_jitter$pheno$cross), qtl = CerEuk_LPmerge_expand4_qtl)
summary(fitqtl(CerEuk_LPmerge_jitter, qtl=CerEuk_LPmerge_expand4_qtl, dropone=TRUE, get.ests=TRUE, covar=data.frame(cross=CerEuk_LPmerge_jitter$pheno$cross), pheno.col="pr", formula= y ~ Q1 + Q2 + Q3 + Q4 + Q5 + cross))
CerEuk_LPmerge_add_expand4_qtl <- addqtl(CerEuk_LPmerge_jitter, qtl=CerEuk_LPmerge_expand4_qtl, covar=data.frame(cross=CerEuk_LPmerge_jitter$pheno$cross), pheno.col="pr", method="imp", formula= y ~ Q1 + Q2 + Q3 + Q4 + Q5 + cross) 
CerEuk_LPmerge_expand4_qtl_addint<-addint(CerEuk_LPmerge_jitter, pheno.col="pr",  qtl=CerEuk_LPmerge_expand4_qtl, method="imp")

CerEuk_LPmerge_expand5_qtl<-addtoqtl(CerEuk_LPmerge_jitter, CerEuk_LPmerge_expand4_qtl, 4, 40)
summary(fitqtl(CerEuk_LPmerge_jitter, qtl=CerEuk_LPmerge_expand5_qtl, dropone=TRUE, get.ests=TRUE, covar=data.frame(cross=CerEuk_LPmerge_jitter$pheno$cross), pheno.col="pr", formula= y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + cross))
CerEuk_LPmerge_expand5_qtl<-refineqtl(CerEuk_LPmerge_jitter, pheno.col="pr", method="imp", covar=data.frame(cross=CerEuk_LPmerge_jitter$pheno$cross), qtl = CerEuk_LPmerge_expand5_qtl)
summary(fitqtl(CerEuk_LPmerge_jitter, qtl=CerEuk_LPmerge_expand5_qtl, dropone=TRUE, get.ests=TRUE, covar=data.frame(cross=CerEuk_LPmerge_jitter$pheno$cross), pheno.col="pr", formula= y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + cross))
CerEuk_LPmerge_add_expand5_qtl <- addqtl(CerEuk_LPmerge_jitter, qtl=CerEuk_LPmerge_expand5_qtl, covar=data.frame(cross=CerEuk_LPmerge_jitter$pheno$cross), pheno.col="pr", method="imp", formula= y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + cross)
CerEuk_LPmerge_expand5_qtl_addint<-addint(CerEuk_LPmerge_jitter, pheno.col="pr",  qtl=CerEuk_LPmerge_expand5_qtl, method="imp")

CerEuk_LPmerge_expand6_qtl<-addtoqtl(CerEuk_LPmerge_jitter, CerEuk_LPmerge_expand5_qtl, 6, 5.0)
summary(fitqtl(CerEuk_LPmerge_jitter, qtl=CerEuk_LPmerge_expand6_qtl, dropone=TRUE, get.ests=TRUE, covar=data.frame(cross=CerEuk_LPmerge_jitter$pheno$cross), pheno.col="pr", formula= y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + cross))
CerEuk_LPmerge_expand6_qtl<-refineqtl(CerEuk_LPmerge_jitter, pheno.col="pr", method="imp", qtl = CerEuk_LPmerge_expand6_qtl)
summary(fitqtl(CerEuk_LPmerge_jitter, qtl=CerEuk_LPmerge_expand6_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", covar=data.frame(cross=as.numeric(CerEuk_LPmerge$pheno$cross)), formula= y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + cross ))
CerEuk_LPmerge_add_expand6_qtl <- addqtl(CerEuk_LPmerge_jitter, qtl=CerEuk_LPmerge_expand6_qtl, pheno.col="pr", covar=data.frame(cross=CerEuk_LPmerge_jitter$pheno$cross), method="imp", formula= y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + cross)
CerEuk_LPmerge_expand6_qtl_addint<-addint(CerEuk_LPmerge_jitter, pheno.col="pr",  qtl=CerEuk_LPmerge_expand6_qtl, method="imp")

cairo_pdf("CerEukLPmerge_fitqtl.pdf", width=8, height=3.5)
plotLodProfile(CerEuk_LPmerge_expand6_qtl,showallchr=TRUE); abline(h=3.154491)
dev.off()


rbind(
lodint(CerEuk_LPmerge_expand6_qtl, qtl.index=1, drop=1, expandtomarker=FALSE),
lodint(CerEuk_LPmerge_expand6_qtl, qtl.index=1, drop=1.5, expandtomarker=FALSE),
lodint(CerEuk_LPmerge_expand6_qtl, qtl.index=2, drop=1, expandtomarker=FALSE),
lodint(CerEuk_LPmerge_expand6_qtl, qtl.index=2, drop=1.5, expandtomarker=FALSE),
lodint(CerEuk_LPmerge_expand6_qtl, qtl.index=3, drop=1, expandtomarker=FALSE),
lodint(CerEuk_LPmerge_expand6_qtl, qtl.index=3, drop=1.5, expandtomarker=FALSE),
lodint(CerEuk_LPmerge_expand6_qtl, qtl.index=4, drop=1, expandtomarker=FALSE),
lodint(CerEuk_LPmerge_expand6_qtl, qtl.index=4, drop=1.5, expandtomarker=FALSE),
lodint(CerEuk_LPmerge_expand6_qtl, qtl.index=5, drop=1, expandtomarker=FALSE),
lodint(CerEuk_LPmerge_expand6_qtl, qtl.index=5, drop=1.5, expandtomarker=FALSE),
lodint(CerEuk_LPmerge_expand6_qtl, qtl.index=6, drop=1, expandtomarker=FALSE),
lodint(CerEuk_LPmerge_expand6_qtl, qtl.index=6, drop=1.5, expandtomarker=FALSE),
lodint(CerEuk_LPmerge_expand6_qtl, qtl.index=7, drop=1, expandtomarker=FALSE),
lodint(CerEuk_LPmerge_expand6_qtl, qtl.index=7, drop=1.5, expandtomarker=FALSE))

cairo_pdf("CerEuk_LPmerge_effectplots.pdf", height=20, width=3)
par(mfrow = c(7,1))
effectplot(CerEuk_LPmerge_jitter,mname1="1@79.9", pheno.col="pr", ylim=c(2.5,3.5))
effectplot(CerEuk_LPmerge_jitter,mname1="2@59", pheno.col="pr", ylim=c(2.5,3.5))
effectplot(CerEuk_LPmerge_jitter,mname1="3@71", pheno.col="pr", ylim=c(2.5,3.5))
effectplot(CerEuk_LPmerge_jitter,mname1="4@59.8", pheno.col="pr", ylim=c(2.5,3.5))
effectplot(CerEuk_LPmerge_jitter,mname1="5@36.9", pheno.col="pr", ylim=c(2.5,3.5))
effectplot(CerEuk_LPmerge_jitter,mname1="6@5", pheno.col="pr", ylim=c(2.5,3.5))
effectplot(CerEuk_LPmerge_jitter,mname1="X@38", pheno.col="pr", ylim=c(2.5,3.5))
dev.off()

true_loci_LPmerge<-otto_jones_ci(D=0.8381,M=0.0872,nd=7,alpha=0.05,amin=0.0335, res=4, res2=2, max.loci=100)


## export QTL scaffolds for enrichement (all scaffolds within 1-LOD interval from the peak)
cereuk_QTL_1LODinterval<-data.frame(locus=character(), position=numeric(), chr=factor())
for(i in 1:7) {
	temp_df<-data.frame(locus=names(pull.map(CerEuk_LPmerge_jitter)[[as.character(lodint(CerEuk_LPmerge_expand6_qtl,qtl.index=i,drop=1,expandtomarkers=FALSE)$chr[1])]]),position=matrix(pull.map(CerEuk_LPmerge_jitter)[[as.character(lodint(CerEuk_LPmerge_expand6_qtl,qtl.index=i,drop=1,expandtomarkers=FALSE)$chr[1])]]), chr=lodint(CerEuk_LPmerge_expand6_qtl,qtl.index=i,drop=1,expandtomarkers=FALSE)$chr[1])

	from=min(lodint(CerEuk_LPmerge_expand6_qtl,qtl.index=i,expandtomarkers=FALSE)$pos)
	to=max(lodint(CerEuk_LPmerge_expand6_qtl,qtl.index=i,expandtomarkers=FALSE)$pos)

	cereuk_QTL_1LODinterval<-rbind(cereuk_QTL_1LODinterval,temp_df[max(which(temp_df$position<=from)):min(which(temp_df$position>=to)),])
	}
	
cereuk_QTL_1LODinterval$scaffold<-gsub("_.*","",cereuk_QTL_1LODinterval$locus)

write.table(cereuk_QTL_1LODinterval,"cereuk_QTL_1LODinterval.txt", quote=FALSE, row.names=FALSE, sep="\t")

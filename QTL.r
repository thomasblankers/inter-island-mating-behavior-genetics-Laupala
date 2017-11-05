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
cim_CerEukHC1_hk_3<-cim(CerEukHC1_jitter, pheno.col="pr", method="hk", n.marcovar=3, window=20)
cim_CerEukHC1_hk_4<-cim(CerEukHC1_jitter, pheno.col="pr", method="hk", n.marcovar=4, window=20)
cim_CerEukHC1_hk_5<-cim(CerEukHC1_jitter, pheno.col="pr", method="hk", n.marcovar=5, window=20)

cim_CerEukHC1_hk_3_1000perm<-cim(CerEukHC1_jitter, pheno.col="pr", method="hk", n.marcovar=3, window=20, n.perm=1000)

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
plot(twoqtl_CerEukHC1_imp, chr=5, lower="cond-int", upper="cond-add")
plot(twoqtl_CerEukHC1_imp, chr=3, lower="cond-int", upper="cond-add")
plot(twoqtl_CerEukHC1_imp, chr=4, lower="cond-int", upper="cond-add")r
plot(twoqtl_CerEukHC1_imp, chr="X", lower="cond-int", upper="cond-add")

twoqtl_CerEukHC1_hk_1000perm<-scantwo(CerEukHC1_jitter, pheno.col="pr", method="hk", chr=qtl_CerEukHC1, clean.output=TRUE, n.perm=1000, perm.Xsp=FALSE) # takes about 40 minutes
penalties=calc.penalties(twoqtl_CerEukHC1_hk_1000perm, alpha=0.05)

    main    heavy    light 
3.446001 5.603086 3.434136

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

CerEukHC1_expand2_qtl<-addtoqtl(CerEukHC1_jitter,CerEukHC1_expand_qtl, 4, 58)
summary(fitqtl(CerEukHC1_jitter, qtl=CerEukHC1_expand2_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula= y ~ Q1 + Q2 + Q3))# not significant
CerEukHC1_expand2_qtl<-refineqtl(CerEukHC1_jitter, pheno.col="pr", method="imp", qtl = CerEukHC1_expand2_qtl)
summary(fitqtl(CerEukHC1_jitter, qtl=CerEukHC1_expand2_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula= y ~ Q1 + Q2 + Q3))# not significant
CerEukHC1_add_expand2_qtl <- addqtl(CerEukHC1_jitter, qtl=CerEukHC1_expand2_qtl, pheno.col="pr", method="imp", formula= y ~ Q1 + Q2 + Q3) # not significant
CerEukHC1_expand2_qtl_addint<-addint(CerEukHC1_jitter, pheno.col="pr",  qtl=CerEukHC1_expand2_qtl, method="imp")


rbind(
lodint(CerEukHC1_expand2_qtl, drop=1, qtl.index=1),
lodint(CerEukHC1_expand2_qtl, drop=2, qtl.index=1),
lodint(CerEukHC1_expand2_qtl, drop=1, qtl.index=2),
lodint(CerEukHC1_expand2_qtl, drop=2, qtl.index=2),
lodint(CerEukHC1_expand2_qtl, drop=1, qtl.index=3),
lodint(CerEukHC1_expand2_qtl, drop=2, qtl.index=3))


true_loci_HC1<-otto_jones_ci(D=0.8381,M=0.1124,nd=5,alpha=0.05,amin=0.02852,res=4,res2=2,max.loci=100)

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

cim_CerEukHC2_imp<-cim(CerEukHC2_jitter, pheno.col="pr", method="imp", n.marcovar=3, window=20)
cim_CerEukHC2_hk<-cim(CerEukHC2_jitter, pheno.col="pr", method="hk", n.marcovar=3, window=20)

cim_CerEukHC2_hk_2<-cim(CerEukHC2_jitter, pheno.col="pr", method="hk", n.marcovar=2, window=20)
cim_CerEukHC2_hk_3<-cim(CerEukHC2_jitter, pheno.col="pr", method="hk", n.marcovar=3, window=20)
cim_CerEukHC2_hk_4<-cim(CerEukHC2_jitter, pheno.col="pr", method="hk", n.marcovar=4, window=20)
cim_CerEukHC2_hk_5<-cim(CerEukHC2_jitter, pheno.col="pr", method="hk", n.marcovar=5, window=20)

cim_CerEukHC2_hk_4_1000perm<-cim(CerEukHC2_jitter, pheno.col="pr", method="hk", n.marcovar=4, window=20, n.perm=1000)


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
plot(twoqtl_CerEukHC2_imp, chr=3, lower="cond-int", upper="cond-add")
plot(twoqtl_CerEukHC2_imp, chr=4, lower="cond-int", upper="cond-add")
plot(twoqtl_CerEukHC2_imp, chr="X", lower="cond-int", upper="cond-add")

twoqtl_CerEukHC2_hk_1000perm<-scantwo(CerEukHC2_jitter, pheno.col="pr", method="hk", chr=qtl_CerEukHC2, clean.output=TRUE, n.perm=1000, perm.Xsp=FALSE) # takes about 40 minutes
penalties=calc.penalties(twoqtl_CerEukHC2_hk_1000perm, alpha=0.05)

    main    heavy    light 
3.350607 5.530137 3.396078

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


CerEukHC2_expand4_qtl<-addtoqtl(CerEukHC2_jitter,CerEukHC2_expand3_qtl, 6, 55) # not significant
summary(fitqtl(CerEukHC2_jitter, qtl=CerEukHC2_expand4_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula= y ~ Q1 + Q2 + Q3 + Q4 + Q5)) 
CerEukHC2_expand4_qtl<-refineqtl(CerEukHC2_jitter, pheno.col="pr", method="imp", qtl = CerEukHC2_expand4_qtl)
summary(fitqtl(CerEukHC2_jitter, qtl=CerEukHC2_expand4_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula= y ~ Q1 + Q2 + Q3 + Q4 + Q5))
CerEukHC2_add_expand4_qtl <- addqtl(CerEukHC2_jitter, qtl=CerEukHC2_expand4_qtl, pheno.col="pr", method="imp", formula= y ~ Q1 + Q2 + Q3 + Q4 + Q5)
CerEukHC2_expand4_qtl_addint<-addint(CerEukHC2_jitter, pheno.col="pr",  qtl=CerEukHC2_expand4_qtl, method="imp")

CerEukHC2_expand5_qtl<-addtoqtl(CerEukHC2_jitter,CerEukHC2_expand4_qtl, 4, 61) # not significant in this model, but significant in scantwo and scanone
summary(fitqtl(CerEukHC2_jitter, qtl=CerEukHC2_expand5_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula= y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6)) 
CerEukHC2_expand5_qtl<-refineqtl(CerEukHC2_jitter, pheno.col="pr", method="imp", qtl = CerEukHC2_expand5_qtl)
summary(fitqtl(CerEukHC2_jitter, qtl=CerEukHC2_expand5_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", formula= y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6))
#CerEukHC2_add_expand5_qtl <- addqtl(CerEukHC2_jitter, qtl=CerEukHC2_expand5_qtl, pheno.col="pr", method="imp", formula= y ~ Q1 + Q2 + Q3 + Q4 + Q5)
#CerEukHC2_expand5_qtl_addint<-addint(CerEukHC2_jitter, pheno.col="pr",  qtl=CerEukHC2_expand5_qtl, method="imp")


rbind(
lodint(CerEukHC2_expand4_qtl, drop=1,  qtl.index=1),
lodint(CerEukHC2_expand4_qtl, drop=2,  qtl.index=1),
lodint(CerEukHC2_expand4_qtl, drop=1, qtl.index=2),
lodint(CerEukHC2_expand4_qtl, drop=2, qtl.index=2),
lodint(CerEukHC2_expand4_qtl, drop=1, qtl.index=3),
lodint(CerEukHC2_expand4_qtl, drop=2, qtl.index=3),
lodint(CerEukHC2_expand4_qtl, drop=1, qtl.index=4),
lodint(CerEukHC2_expand4_qtl, drop=2, qtl.index=4),
lodint(CerEukHC2_expand4_qtl, drop=1, qtl.index=5),
lodint(CerEukHC2_expand4_qtl, drop=2, qtl.index=5))


rbind(effectplot(CerEukHC2_jitter,mname1="1@67.1", pheno.col="pr", draw=FALSE)$Means,
effectplot(CerEukHC2_jitter,mname1="2@57", pheno.col="pr", draw=FALSE)$Means,
effectplot(CerEukHC2_jitter,mname1="5@48.9", pheno.col="pr", draw=FALSE)$Means,
effectplot(CerEukHC2_jitter,mname1="6@56", pheno.col="pr", draw=FALSE)$Means,
effectplot(CerEukHC2_jitter,mname1="X@32", pheno.col="pr", draw=FALSE)$Means,
effectplot(CerEukHC2_jitter,mname1="4@61.8", pheno.col="pr", draw=FALSE)$Means)


true_loci_HC2<-otto_jones_ci(D=0.8381,M=0.1586,n=2,alpha=0.05,amin=0.1026,res=4,res2=2,max.loci=100)

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

cim_CerEuk_LPmerge_hk_2<-cim(CerEuk_LPmerge_jitter, pheno.col="pr", method="hk", n.marcovar=2, window=20)
cim_CerEuk_LPmerge_hk_3<-cim(CerEuk_LPmerge_jitter, pheno.col="pr", method="hk", n.marcovar=3, window=20)
cim_CerEuk_LPmerge_hk_4<-cim(CerEuk_LPmerge_jitter, pheno.col="pr", method="hk", n.marcovar=4, window=20)
cim_CerEuk_LPmerge_hk_5<-cim(CerEuk_LPmerge_jitter, pheno.col="pr", method="hk", n.marcovar=5, window=20)
cim_CerEuk_LPmerge_hk_6<-cim(CerEuk_LPmerge_jitter, pheno.col="pr", method="hk", n.marcovar=6, window=20)

cim_CerEuk_LPmerge_hk_4_1000perm<-cim(CerEuk_LPmerge_jitter, pheno.col="pr", method="hk", n.marcovar=4, window=20, n.perm=1000)


qtl_CerEuk_LPmerge<-c(1,2,3,4,5,6,"X")

## scantwo
twoqtl_CerEuk_LPmerge_imp<-scantwo(CerEuk_LPmerge_jitter, pheno.col="pr", method="imp", chr=qtl_CerEuk_LPmerge, clean.output=TRUE)
twoqtl_CerEuk_LPmerge_hk<-scantwo(CerEuk_LPmerge_jitter, pheno.col="pr", method="hk", chr=qtl_CerEuk_LPmerge, clean.output=TRUE)

# check for potential interactions:

plot(twoqtl_CerEuk_LPmerge_imp, lower="cond-int")
plot(twoqtl_CerEuk_LPmerge_hk, lower="cond-int")

# potential interaction at c(1,5)

# check for secondary peaks at specific chromosomes:

plot(twoqtl_CerEuk_LPmerge_hk, chr=1, lower="cond-int", upper="cond-add")
plot(twoqtl_CerEuk_LPmerge_hk, chr=2, lower="cond-int", upper="cond-add")
plot(twoqtl_CerEuk_LPmerge_hk, chr=5, lower="cond-int", upper="cond-add")
plot(twoqtl_CerEuk_LPmerge_hk, chr=3, lower="cond-int", upper="cond-add")
plot(twoqtl_CerEuk_LPmerge_hk, chr=4, lower="cond-int", upper="cond-add")
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

CerEuk_LPmerge_expand6_qtl<-addtoqtl(CerEuk_LPmerge_jitter, CerEuk_LPmerge_expand5_qtl, 6, 50)
summary(fitqtl(CerEuk_LPmerge_jitter, qtl=CerEuk_LPmerge_expand6_qtl, dropone=TRUE, get.ests=TRUE, covar=data.frame(cross=CerEuk_LPmerge_jitter$pheno$cross), pheno.col="pr", formula= y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + cross))
CerEuk_LPmerge_expand6_qtl<-refineqtl(CerEuk_LPmerge_jitter, pheno.col="pr", method="imp", qtl = CerEuk_LPmerge_expand6_qtl)
summary(fitqtl(CerEuk_LPmerge_jitter, qtl=CerEuk_LPmerge_expand6_qtl, dropone=TRUE, get.ests=TRUE, pheno.col="pr", covar=data.frame(cross=as.numeric(CerEuk_LPmerge$pheno$cross)), formula= y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + cross ))
CerEuk_LPmerge_add_expand6_qtl <- addqtl(CerEuk_LPmerge_jitter, qtl=CerEuk_LPmerge_expand6_qtl, pheno.col="pr", covar=data.frame(cross=CerEuk_LPmerge_jitter$pheno$cross), method="imp", formula= y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + cross)
CerEuk_LPmerge_expand6_qtl_addint<-addint(CerEuk_LPmerge_jitter, pheno.col="pr",  qtl=CerEuk_LPmerge_expand6_qtl, method="imp")



rbind(
lodint(CerEuk_LPmerge_expand6_qtl, qtl.index=1, drop=1, expandtomarker=TRUE),
lodint(CerEuk_LPmerge_expand6_qtl, qtl.index=1, drop=2, expandtomarker=TRUE),
lodint(CerEuk_LPmerge_expand6_qtl, qtl.index=2, drop=1, expandtomarker=TRUE),
lodint(CerEuk_LPmerge_expand6_qtl, qtl.index=2, drop=2, expandtomarker=TRUE),
lodint(CerEuk_LPmerge_expand6_qtl, qtl.index=3, drop=1, expandtomarker=TRUE),
lodint(CerEuk_LPmerge_expand6_qtl, qtl.index=3, drop=2, expandtomarker=TRUE),
lodint(CerEuk_LPmerge_expand6_qtl, qtl.index=4, drop=1, expandtomarker=TRUE),
lodint(CerEuk_LPmerge_expand6_qtl, qtl.index=4, drop=2, expandtomarker=TRUE),
lodint(CerEuk_LPmerge_expand6_qtl, qtl.index=5, drop=1, expandtomarker=TRUE),
lodint(CerEuk_LPmerge_expand6_qtl, qtl.index=5, drop=2, expandtomarker=TRUE),
lodint(CerEuk_LPmerge_expand6_qtl, qtl.index=6, drop=1, expandtomarker=TRUE),
lodint(CerEuk_LPmerge_expand6_qtl, qtl.index=6, drop=2, expandtomarker=TRUE),
lodint(CerEuk_LPmerge_expand6_qtl, qtl.index=7, drop=1, expandtomarker=TRUE),
lodint(CerEuk_LPmerge_expand6_qtl, qtl.index=7, drop=2, expandtomarker=TRUE))


rbind(effectplot(CerEuk_LPmerge_jitter,mname1="1@79.9", pheno.col="pr", draw=FALSE)$Means,
effectplot(CerEuk_LPmerge_jitter,mname1="2@59", pheno.col="pr", draw=FALSE)$Means,
effectplot(CerEuk_LPmerge_jitter,mname1="3@71", pheno.col="pr", draw=FALSE)$Means,
effectplot(CerEuk_LPmerge_jitter,mname1="4@59.8", pheno.col="pr", draw=FALSE)$Means,
effectplot(CerEuk_LPmerge_jitter,mname1="5@36.9", pheno.col="pr", draw=FALSE)$Means,
effectplot(CerEuk_LPmerge_jitter,mname1="6@5", pheno.col="pr", draw=FALSE)$Means,
effectplot(CerEuk_LPmerge_jitter,mname1="X@38", pheno.col="pr", draw=FALSE)$Means)


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


####### annotate the linkage map #####

### see transcriptome notes for transcriptome
setwd("E:/Data/Laupala/RNAseq_Transcriptome/v057")
scaffold_transcripts<-read.delim("LaCer057_Scaffolds_Transcripts_NoDuplicates.txt")
annotation<-read.delim("LaCer057_transcriptome_annotation.txt")
scaffold_transcripts_annotation<-merge(scaffold_transcripts,annotation,"transcript", all.x=TRUE)
setwd("E:/Data/Laupala/Rwork")
CerEuk_LPmerged_map<-read.delim("CerEuk_LPmerge_map.txt")# this is the merged linkage map usign the LPmerge algorithm (See above)
CerEuk_LPmerged_map$scaffold<-paste0("Lko051",CerEuk_LPmerged_map$scaffold)
Genome_scaffold_translation<-read.delim("E:/Data/Laupala/ReferenceAssembly/Lko_v051tov057_translations.txt")
CerEuk_LPmerged_map_translate<-merge(CerEuk_LPmerged_map,Genome_scaffold_translation,"scaffold", all.x=TRUE, sort=FALSE)
# only masked scaffolds, so no need to worry about names or positions
CerEuk_LPmerged_map<-read.delim("CerEuk_LPmerge_map.txt")# this is the merged linkage map usign the LPmerge algorithm (See above)
CerEuk_LPmerged_map$scaffold<-factor(paste0("Lko057",CerEuk_LPmerged_map$scaffold))

CerEuk_LPmerged_map_annotated<-merge(CerEuk_LPmerged_map, scaffold_transcripts_annotation, "scaffold", all.x=TRUE, sort=FALSE) 
write.table(CerEuk_LPmerged_map_annotated,"CerEuk_LPmerged_map_annotated.txt",sep="\t",quote=FALSE, row.names=FALSE)

#merge the map with the QTL scaffold list and limit proximity info to LOD1 interval
CerEuk_LPmerged_map$order<-seq(1:nrow(CerEuk_LPmerged_map))
QTL_genes<-read.delim("CerEuk_QTL_scaffolds.txt")
CerEuk_LPmerged_map<-merge(CerEuk_LPmerged_map,QTL_genes,"scaffold",all.x=TRUE)
CerEuk_LPmerged_map<-CerEuk_LPmerged_map[order(CerEuk_LPmerged_map$order),]
CerEuk_LPmerged_map[which(CerEuk_LPmerged_map$proximity=="LOD2"),"proximity"]<-NA
CerEuk_LPmerged_map_annotated<-merge(CerEuk_LPmerged_map, scaffold_transcripts_annotation, "scaffold", all.x=TRUE, sort=FALSE) 
CerEuk_LPmerged_map_annotated<-CerEuk_LPmerged_map_annotated[order(CerEuk_LPmerged_map_annotated$order),]
CerEuk_LPmerged_map_annotated$LG<-CerEuk_LPmerged_map_annotated$LG.x
CerEuk_LPmerged_map_annotated$LG.y<-NULL


#### Are QTL regions enriched for a particular GO? ####

source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("topGO")
library("topGO")


geneID2GO <- readMappings(file = "E:/Data/Laupala/RNAseq_Transcriptome/v057/LaCer057_transcriptome_annotation_GOs.txt")
geneUniverse <- names(geneID2GO) # all genes with GO term
QTL_genes<-read.delim("E:/Data/Laupala/CerEukCross/CerEuk_QTL_scaffolds_transcripts057_annotated.txt") # target genes/transcripts of interest
QTL_genes<-QTL_genes[-which(QTL_genes$proximity == "LOD2"),] # target genes/transcripts of interest ,only keep genes within 1 LOD of the peak

geneList <- factor(as.integer(geneUniverse %in% as.vector(QTL_genes$transcript))) # specify for each gene in the gene universe whether it is present in the list of target genes

names(geneList) <- geneUniverse

myGOdata_BP <- new("topGOdata", description="CerEuk_QTL", ontology="BP",allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO) # create GOdata object for BP, MF, CC
resultFisher_BP_pc <- runTest(myGOdata_BP, algorithm="parentchild", statistic="fisher") # parentchild (Grossman et al 2007) assess the annotation of the parental nodes of the significant term to correct for the inheritance problem

allRes_BP_pc <- GenTable(myGOdata_BP, pcFisher = resultFisher_BP_pc, topNodes=6426) # get all nodes
allRes_BP_pc <- allRes_BP_pc[which(allRes_BP_pc$Significant > 0),] # only keep the GOs present in the outlier set
allRes_BP_pc$fdr_pvalue<-p.adjust(allRes_BP_pc$pcFisher,"fdr",nrow(allRes_BP_pc)) # do FDR correction, not really neccessary, but conservative

cairo_pdf(width=20, height=30,"E:/Data/Laupala/CerEukCross/CerEuk_GOenrichment_SigNodes.pdf")
showSigOfNodes(myGOdata_BP, score(resultFisher_BP_pc), firstSigNodes = 10, useInfo = "all") # plot a scheme connecting the enriched GO terms with their parent terms - other nodes in the network
dev.off()
showSigOfNodes(myGOdata_BP, score(resultFisher_BP_pc), firstSigNodes = 10, wantedNodes=c(, putWN=TRUE useInfo = "all") # plot a scheme connecting the enriched GO terms with their parent terms - other nodes in the network

## what are some of the candidate loci?
sg <- sigGenes(myGOdata_BP)

locomotion = "GO:0031987"
locomotion_genes <- genesInTerm(myGOdata_BP, locomotion)[[1]]
locomotion_genes_sg <-locomotion_genes[which(locomotion_genes %in% sg)]
locomotion_genes_sg <- merge(data.frame(transcript=locomotion_genes_sg),annotation, all.x=TRUE)

brain = "GO:0048854"
brain_genes <- genesInTerm(myGOdata_BP, brain)[[1]]
brain_genes_sg <-brain_genes[which(brain_genes %in% sg)]
brain_genes_sg <- merge(data.frame(transcript=brain_genes_sg),annotation, all.x=TRUE)

photoreceptor = "GO:0048056"
photoreceptor_genes <- genesInTerm(myGOdata_BP, photoreceptor)[[1]]
photoreceptor_genes_sg <-photoreceptor_genes[which(photoreceptor_genes %in% sg)]
photoreceptor_genes_sg <- merge(data.frame(transcript=photoreceptor_genes_sg),annotation, all.x=TRUE)

neuromuscular = "GO:0050905"
neuromuscular_genes <- genesInTerm(myGOdata_BP, neuromuscular)[[1]]
neuromuscular_genes_sg <-neuromuscular_genes[which(neuromuscular_genes %in% sg)]
neuromuscular_genes_sg <- merge(data.frame(transcript=neuromuscular_genes_sg),annotation, all.x=TRUE)


neurotransmitter = "GO:0007269"
neurotransmitter_genes <- genesInTerm(myGOdata_BP, neurotransmitter)[[1]]
neurotransmitter_genes_sg <-neurotransmitter_genes[which(neurotransmitter_genes %in% sg)]
neurotransmitter_genes_sg <- merge(data.frame(transcript=neurotransmitter_genes_sg),annotation, all.x=TRUE)

reproduction = "GO:0044702"
reproduction_genes <- genesInTerm(myGOdata_BP, reproduction)[[1]]
reproduction_genes_sg <-reproduction_genes[which(reproduction_genes %in% sg)]
reproduction_genes_sg <- merge(data.frame(transcript=reproduction_genes_sg),annotation, all.x=TRUE)

steroid = "GO:0006694"
steroid_genes <- genesInTerm(myGOdata_BP, steroid)[[1]]
steroid_genes_sg <-steroid_genes[which(steroid_genes %in% sg)]
steroid_genes_sg <- merge(data.frame(transcript=steroid_genes_sg),annotation, all.x=TRUE)


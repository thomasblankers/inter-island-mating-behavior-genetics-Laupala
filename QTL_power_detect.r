## these are function to analyze the power to detect QTL and to estimate the true number of loci in QTL mapping studies using intercross data (all formulae) and backcross (only Otto & Jones formulae)
## lynch_walsh_sample calculates the minimum F2 sample size needed to detect a QTL of a given effect size under different type 1 and type 2 errors (equation 15.37, Lynch & Wlash, 1998, Genetics and Analysis of Quantitative Traits)
## lynch_walsh_detect approximates the smallest effect size one can detect in a QTL experiment given the F2 sample size based on Lynch & Walsh, 1998
## otto_jones_detect approximates the smalles effect size one can detect (threshold theta) based on equation 11 in Otto & Jones, 2000, Detecting the undetected: estimating the total number of loci underlying a quantitative trait, Genetics 156.
## otto_jones_trueQTLnumber calculates the expected true number of loci underlying a trait given the approximated theta, the minimum and average effect size, the parental difference delta z, and the number of detected QTL. Note that these values can be both standardized (delta Z = 0.5, and minimum effect = a/parental_difference) and on the observed scale (latter is preferred)
## otto_jones_ci estimates the confidence limits around otto_jones_trueQTLnumber
## I recommend trying the examples in Example 11 in Chapeter 15 of Lynch and Walsh to get a feel for the first two functions and the example in Table 5 in the Otto & Jones paper to get a feel for how the latter two functions work.

## lynch & walsh sample and effect size base on standardized effects (fractions of F2 segregation variance r2f2)
lynch_walsh_sample<-function(effect.size=0.5, alpha=0.05, beta=0.1, k=0) {
#effect.size = expected effect size of QTL
# alpha = type I error, false positive rate
# beta = typ II error, 1-expected power to detect a QTL
# k = dominant/recessive effect; 0 means fully additive, positive implies dominance effects, negative recessive effects
	 prob.beta=qnorm(1-beta, mean = 0, sd = 1, log = FALSE)
	 prob.alpha=qnorm(1-(0.5*alpha), mean = 0, sd = 1, log = FALSE)
	 sample.size=((1-effect.size)/effect.size)*(((prob.alpha/(sqrt(1-effect.size)))+prob.beta)^2)*(1+((k^2)/2))
	 sample.size
	 }
	 
lynch_walsh_detect<-function(sample.size=100, alpha=0.05, beta=0.1, k=0,res=0.001) {
#sample.size=sample size of QTL experiment
#res is resolution at which to estimate the minimum effect size of detectable QTL

	prob.beta=qnorm(1-beta, mean = 0, sd = 1, log = FALSE)
	prob.alpha=qnorm(1-(0.5*alpha), mean = 0, sd = 1, log = FALSE)
	vector.sample.size<-data.frame(effectsize=1, samplesize=lynch_walsh_sample(1,alpha,beta,k))
	for(i in seq(from=1-res, to=0+res, by=-res)) {
		effect.size=i
		vector.sample.size=rbind(vector.sample.size,data.frame(effectsize=i, samplesize=lynch_walsh_sample(i,alpha,beta,k)))
		}
	detect=min(vector.sample.size[which(vector.sample.size$samplesize<sample.size),"effectsize"])
	detect
	}

lynch_walsh_detect<-function(sample.size=100, alpha=0.05, beta=0.1, k=0,res=0.001) {
#sample.size=sample size of QTL experiment
#res is resolution at which to estimate the minimum effect size of detectable QTL

	prob.beta=qnorm(1-beta, mean = 0, sd = 1, log = FALSE)
	prob.alpha=qnorm(1-(0.5*alpha), mean = 0, sd = 1, log = FALSE)
	vector.sample.size<-data.frame(effectsize=1, samplesize=lynch_walsh_sample(1,alpha,beta,k))
	for(i in seq(from=1-res, to=0+res, by=-res)) {
		effect.size=i
		vector.sample.size=rbind(vector.sample.size,data.frame(effectsize=i, samplesize=lynch_walsh_sample(i,alpha,beta,k)))
		}
	detect=min(vector.sample.size[which(vector.sample.size$samplesize<sample.size),"effectsize"])
	detect
	}

## alternatives for lynch & walsh sample size using the measurement scale of the phenotype (a, where a is the additive effect, ie the difference between the average phenotypic value in one parent relative to the hybrid average):
lynch_walsh_sample<-function(a=1, f2var= 6.15, alpha=0.05, beta=0.1, k=0, crosstype="backcross") {
#a = expected effect size of QTL
# alpha = type I error, false positive rate
# beta = typ II error, 1-expected power to detect a QTL
# k = dominant/recessive effect; 0 means fully additive, positive implies dominance effects, negative recessive effects
	 prob.beta=qnorm(1-beta, mean = 0, sd = 1, log = FALSE)
	 prob.alpha=qnorm(1-(0.5*alpha), mean = 0, sd = 1, log = FALSE)
	 if(crosstype=="intercross") {
		r2f2=(((2*a)/f2var)^2)/8
		sample.size=((1-r2f2)/r2f2)*(((prob.alpha/(sqrt(1-r2f2)))+prob.beta)^2)*(1+((k^2)/2))
		}
	 if(crosstype=="backcross") {
		r2f2=(((a*(1-k))/f2var)^2)/4
		sample.size=((1-r2f2)/r2f2)*(((prob.alpha/(sqrt(1-r2f2)))+prob.beta)^2)
		}
	sample.size
	 }
	 


otto_jones_detect<-function(amin=0.03,M=0.09,nd=7, res=4) {
#theta = detectability threshold	
#D = sum of additive effects across all possible QTL, ie half the parental difference on the original scale (in which case amin should also be on the scale of the measured variable) or standardized ( in which case, always 0.5)
#amin = minimum effect size that was detected in a QTL experiment
#nd = number of detected QTL
#M = average effect size of detected QTL
#res = the resolution used to approximate zero, the larger the more precise, but also the more computation time involved. 

equation11=function(theta){amin-theta-((M-theta)/nd)+((amin*exp((-1*amin*nd)/(M-theta)))/(1-exp((-1*amin*nd)/(M-theta))))}
	repeat{
		theta=sample(seq(from=0, to=amin, by=1*10^(-1*res)),1)
		solution=round(equation11(theta),res)
		if(solution == 0) break
		}
	return(theta)
	}
	
otto_jones_trueQTLnumber<-function(D=0.83,amin=0.03,M=0.09,nd=7, res=4) {
	nqtl=D/(M-otto_jones_detect(amin,M,nd,res))
	return(nqtl)
	}

otto_jones_ci<-function(D= 0.8381, M= 0.0872,nd= 7, alpha= 0.05, amin= 0.0335, res= 4, res2=2, max.loci=100) {
#alpha = significance level, type 1 error
#res = resolution of the detectability threshold and true loic estimation (3 < res < 6 should be relatively quickly and accurate)
#res2 = resolution of the approximation of the solution for equation 9 in Otto & Jones; high values make this particularly slow, try res2=2 or res2=3)
	theta=otto_jones_detect(amin=amin,M=M,nd=nd,res=res)
	n=otto_jones_trueQTLnumber(D=D,amin=amin,M=M,nd=nd, res=res)
	X=qchisq(1-alpha, df=1)
	equation9<-function(n){((-1)*nd) + (((M-theta)*n*nd)/D) - (nd*log(((M-theta)*n)/D)) - (X/2)}
	n_vec<-data.frame(number=1, solution=equation9(1), sign=sign(equation9(1)))
	for( i in seq(from=1+1*10^(-1*res2), to=max.loci, by=1*10^(-1*res2))) {
		n_vec<-rbind(n_vec,data.frame(number=i, solution=equation9(i), sign=sign(equation9(i))))
		}
	plustominus <- which(diff(n_vec$sign)<0)
	minustoplus <- which(diff(n_vec$sign)>0)
	solutions=c(mean(n_vec[plustominus:plustominus+1,"number"]),mean(n_vec[minustoplus:minustoplus+1,"number"]))
	tau=(M-theta)/M
	Mund=M*(1-((1-tau)/(1-exp((-1*(1-tau))/tau)))) # is the effect size of undetected QTL
	result<-list()
	result$theta<-theta; result$Mundetected=Mund; result$trueQTLnumber=n; result$lowerCI=solutions[1]; result$upperCI=solutions[2]
	return(result)
	}


	
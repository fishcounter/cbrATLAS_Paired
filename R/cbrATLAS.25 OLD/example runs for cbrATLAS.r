
### set up for cbrATLAS session
# set up

# save the .First fn, then the files will be automatically cleared & reloaded each time .Rdata is opened.
.First=function(){
	# files have to be loaded individually because R4.5.1 broke it
	invisible({
	rm(list=ls(all=TRUE))
	utils::setWindowTitle(getwd())
	options(contrasts = c("contr.treatment", "contr.poly"), length = 4000)
	library(MASS)
#	istall.packages(c("survival", "flexsurv", "vitality"))
#	install.packages("failCompare_1.2.0.zip")
	library("flexsurv")
	library("vitality")
	source("latest code/AdjSurv multiRel.fn.R")
	source("latest code/atlas2flat.fn.R")
	source("latest code/boot.L.R")
	source("latest code/cjs.fn.R")
	source("latest code/cjs.lik.R")
	source("latest code/cjs.paired.fn.R")
	source("latest code/cjs.taglife.corr.R")
	source("latest code/correct.fn.R")
	source("latest code/dayhr.fn.R")
	source("latest code/harmonic.R")
	source("latest code/mean_tt2site.fn.R")
	source("latest code/paired.cjs.lik.R")
	source("latest code/thist0.R")
	source("latest code/vs2.fn.R")
	})	
}

Below are examples of runs.  
======================================================================================================================
# list of taglife curves that can be selected from failCompare
model=c("weibull","weibull3", "gompertz", "gamma", "lognormal", "llogis", "gengamma","vitality.ku","vitality.4p")

## Sample runs for cbrATLAS
# single release
Note: The current default output is to print the output the screen and return an list that I would like to save in an object.  
	   Save to an object if you only want the output printed at the end.

# allow cbrATLAS to pick best taglife fit
AdjSurv.fn(taglife.file="../microtagtest for Atlas.csv",taghist.file="data/Single_Release_ATLAS_Practice_File_no_censors.csv",taglife.model=NULL, 
	num.release=1,study.type="single",num.bootstrap=1000) 
	
# save output and use weibull to fit taglife	
tmp.out=AdjSurv.fn(taglife.file="data/Taglife_ATLAS_Practice_File.csv",taghist.file="data/Single_Release_ATLAS_Practice_File_no_censors.csv",taglife.model="weibull", 
	    num.release=1,study.type="single",num.bootstrap=1000)

# different taglife curve
AdjSurv.fn(taglife.file="../microtagtest for Atlas.csv",taghist.file="data/Single_Release_ATLAS_Practice_File_no_censors.csv",taglife.model=NULL, 
	num.release=1,study.type="single",num.bootstrap=1000)  # this differs slightly--I think it's due to tag detections beyond the taglife curve--I'll have to ask Jim what ATLAS does with those.

# paired release
AdjSurv.fn(taglife.file="data/Spring2010TagLife.csv",taghist.file="data/Paired_Release_ATLAS_Practice_File.csv",taglife.model="weibull",num.release=2,study.type="paired")

AdjSurv.fn(taglife.file="data/JDay2011_TagLifepooled.csv",taghist.file="data/chinR1_R3_JDay2011.csv",taglife.model="weibull3",num.release=2,study.type="paired") 



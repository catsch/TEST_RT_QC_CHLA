##############################################################################
# Test of Dark
# Catherine Schmechtig 
# September 2020
##############################################################################
library(ncdf4)
library(stringr)
 
source("./read_VSS.R")
source("./RunningFilter.R")
source("./READ_CTD.R")
source("./MLD.R")

uf=commandArgs()

mission  <- uf[2]

liste_to_do=read.table("./liste_all",sep=" ",header=FALSE)

# List of the file to process
LIST_nc=liste_to_do$V1

print(LIST_nc)

# We are working on CHLA 
PARAM_STRING=str_pad("CHLA",64,"right")

# text_file for the whole float
path_out_txt=paste(mission,"_DARK.txt",sep="")

for (IDnc in LIST_nc) {

	# Open the B file 
	filenc=nc_open(IDnc,readunlim=FALSE,write=FALSE)

	# Get the corresponding C file name

	file_in_C=str_replace(IDnc,"/B","/")

	# if B and C are not in the same mode
	if (!file.exists(file_in_C)) file_in_C=str_replace(file_in_C,"profiles/R","profiles/D")
	if (!file.exists(file_in_C)) file_in_C=str_replace(file_in_C,"profiles/D","profiles/R")


	# open the C file 
	filenc_C=nc_open(file_in_C,readunlim=FALSE,write=FALSE)

###################################################################################
####    Read the B file PARAMETER to check the availability of CHLA 
###################################################################################

	PARAMETER=ncvar_get(filenc,"PARAMETER")

	index_param=which(PARAMETER == PARAM_STRING , arr.ind=TRUE)

###	Very IMPORTANT
###	Next iteration if the parameter is not in the file 

	if ( length(index_param)==0 ) {
		
	next

	}

###################################################################################
####    Read the C file and estimate the MLD (for quenching test) 
###################################################################################

#### READ Core file 

	CTD=read_CTD(filenc_C)

# we get        : CTD$PRES
#               : CTD$PSAL
#               : CTD$TEMP

#### Estimation of the MLD 
	
	MLD=CALC_MLD(CTD$PRES, CTD$PSAL , CTD$TEMP)

	if ( is.na(MLD) ) {
		
	next

	}

###################################################################################
####    Read the B file 
###################################################################################

###	studied Profile

	i_param_param =index_param[1]

	i_prof_param = index_param[3]

###	Read the BFILE

	PRES=ncvar_get(filenc,"PRES")

	CHLA=ncvar_get(filenc,"CHLA")

	CYCLE_NUMBER=unique(ncvar_get(filenc,"CYCLE_NUMBER"))

###	working only on the studied profile

	PRES_CHLA=PRES[!is.na(CHLA)]

	CHLA_CHLA=CHLA[!is.na(CHLA)]

	MED_CHLA=rep(NA,length(CHLA_CHLA))

	SPIKE_CHLA_C=rep(FALSE,length(CHLA_CHLA))

	SPIKE_CHLA_A=rep(FALSE,length(CHLA_CHLA))

	RESOLUTION=read_VSS(filenc,i_prof_param)

###     Let s calculate all the different median filters 5,7,11 

	MED_CHLA_5=RunningFilter(2,CHLA_CHLA,na.fill=T, ends.fill=T, Method="Median")

	MED_CHLA_7=RunningFilter(3,CHLA_CHLA,na.fill=T, ends.fill=T, Method="Median")

	MED_CHLA_11=RunningFilter(5,CHLA_CHLA,na.fill=T, ends.fill=T, Method="Median")


###	Calculate the profile of MED_CHLA

	for (i in seq(1,length(CHLA_CHLA))) {

		if ( RESOLUTION[i] < 1 ) MED_CHLA[i]=MED_CHLA_11[i]

		if ( (RESOLUTION[i] >= 1) & (RESOLUTION[i] < 3) ) MED_CHLA[i]=MED_CHLA_7[i]

		if ( RESOLUTION[i] >= 3 ) MED_CHLA[i]=MED_CHLA_5[i]

	}

###	CALCULATE the RESIDUALS

###	Christina's proposal 
	RESID_C=abs(CHLA_CHLA-MED_CHLA)

###	Previous version sliding median of 5 
	RESID_A=abs(CHLA_CHLA-MED_CHLA_5)

###	Calculate the percentile of both methods 

	Q10_C=rep(quantile(RESID_C,0.90),length(CHLA_CHLA))

	Q10_A=rep(2*quantile(RESID_A,0.90),length(CHLA_CHLA))

###	Spike 

	SPIKE_CHLA_C[which(RESID_C>Q10_C)]=TRUE

	SPIKE_CHLA_A[which(RESID_A>Q10_A)]=TRUE

###	Nb spikes 

	NB_SPIKE_C=length(CHLA_CHLA[SPIKE_CHLA_C])

	NB_SPIKE_A=length(CHLA_CHLA[SPIKE_CHLA_A])

###	DARK estimation 

	DARK_C_1=min(CHLA_CHLA[!SPIKE_CHLA_C])

	DEPTH_DARK_C_1=PRES_CHLA[which.min(CHLA_CHLA[!SPIKE_CHLA_C])]

	max_depth=max(PRES_CHLA)

###  	Writing a txt file 
	summary=data.frame(CYCLE_NUMBER,DARK_C_1,DEPTH_DARK_C_1,MLD,max_depth)

	write.table(file=path_out_txt,summary,col.names=F,row.names=F,append=TRUE)
###########################################################################
###	CLOSING the NCFILE
###########################################################################

	nc_close(filenc)

	nc_close(filenc_C)


}




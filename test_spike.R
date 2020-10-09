##############################################################################
# Test of the Spike test 
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
path_out_txt=paste(mission,".txt",sep="")

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

### 	Quenching correction 

###	max de la CHLORO despike dans 0.9*MLD
	NPQ_C_1=max(CHLA_CHLA[!SPIKE_CHLA_C & (PRES_CHLA<0.9*MLD)])

	DEPTH_NPQ_C_1=max(PRES_CHLA[CHLA_CHLA==NPQ_C_1 & (PRES_CHLA<0.9*MLD) & !SPIKE_CHLA_C ])

###	max de la CHLORO filtree 
	NPQ_C_2=max(MED_CHLA[(PRES_CHLA<0.9*MLD)])

	DEPTH_NPQ_C_2=max(PRES_CHLA[MED_CHLA==NPQ_C_2 & (PRES_CHLA<0.9*MLD)])

###	max de la CHLORO despike avec la version actuelle de detection des spikes 
	NPQ_A=max(CHLA_CHLA[!SPIKE_CHLA_A & (PRES_CHLA<0.9*MLD)])

	DEPTH_NPQ_A=max(PRES_CHLA[CHLA_CHLA==NPQ_A & (PRES_CHLA<0.9*MLD) & !SPIKE_CHLA_A])

###	What would be the median value of the CHLA in the quenching Area without quenching 

	MEDIAN_RAW_C_1=median(CHLA_CHLA[PRES_CHLA<DEPTH_NPQ_C_1])

	MEDIAN_RAW_C_2=median(CHLA_CHLA[PRES_CHLA<DEPTH_NPQ_C_2])

	MEDIAN_RAW_A=median(CHLA_CHLA[PRES_CHLA<DEPTH_NPQ_A])

###  	Writing a txt file 
	summary=data.frame(CYCLE_NUMBER,NB_SPIKE_C,NB_SPIKE_A,NPQ_C_1,MEDIAN_RAW_C_1,DEPTH_NPQ_C_1,NPQ_C_2,MEDIAN_RAW_C_2,DEPTH_NPQ_C_2,NPQ_A,MEDIAN_RAW_A,DEPTH_NPQ_A)

###	Adding some plots with the quenching correction 

	CHLA_CHLA_NPQ_C_1=CHLA_CHLA
	
	CHLA_CHLA_NPQ_C_1[PRES_CHLA<DEPTH_NPQ_C_1]=NPQ_C_1

	CHLA_CHLA_NPQ_C_2=CHLA_CHLA
	
	CHLA_CHLA_NPQ_C_2[PRES_CHLA<DEPTH_NPQ_C_2]=NPQ_C_2

	CHLA_CHLA_NPQ_A=CHLA_CHLA
	
	CHLA_CHLA_NPQ_A[PRES_CHLA<DEPTH_NPQ_A]=NPQ_A

	write.table(file=path_out_txt,summary,col.names=F,row.names=F,append=TRUE)
###########################################################################
###	CLOSING the NCFILE
###########################################################################

	nc_close(filenc)

	nc_close(filenc_C)

###########################################################################
##	Some plots : localisation of the spikes
###########################################################################

	path_out_jpeg=paste(substr(IDnc,start=36,stop=49),"jpeg",sep="")

	path_out_zoomjpeg=paste(substr(IDnc,start=36,stop=48),"_zoom.jpeg",sep="")

	path_out_quenchingjpeg=paste(substr(IDnc,start=36,stop=48),"_quench.jpeg",sep="")

	jpeg(file=path_out_zoomjpeg)

#	matplot(CHLA_CHLA,PRES_CHLA,col=8,type="l",ylab="Depth [m]",xlab=expression("Chlorophyll a [mg."*m ^ -3 * "]"),xlim=c(-0.2,max(CHLA_CHLA)+0.5),ylim=rev(c(0, max(PRES_CHLA))))

	matplot(CHLA_CHLA,PRES_CHLA,col=8,type="l",ylab="Depth [m]",xlab=expression("Chlorophyll a [mg."*m ^ -3 * "]"),xlim=c(-0.2,max(CHLA_CHLA)+0.5),ylim=rev(c(0, MLD)))

	matplot(CHLA_CHLA[SPIKE_CHLA_C],PRES_CHLA[SPIKE_CHLA_C],type="p",pch=1, col=1,cex=2,add=TRUE)

	matplot(CHLA_CHLA[SPIKE_CHLA_A],PRES_CHLA[SPIKE_CHLA_A],type="p",pch=1, col=2,cex=3,add=TRUE)
	
	legend("bottomright",c("Chl-a","Spike_C","Spike_A"),pch=c(20,20,20),col=c(8,1,2))

	dev.off()

	jpeg(file=path_out_jpeg)

	matplot(CHLA_CHLA,PRES_CHLA,col=8,type="l",ylab="Depth [m]",xlab=expression("Chlorophyll a [mg."*m ^ -3 * "]"),xlim=c(-0.2,max(CHLA_CHLA)+0.5),ylim=rev(c(0, max(PRES_CHLA))))

	matplot(CHLA_CHLA[SPIKE_CHLA_C],PRES_CHLA[SPIKE_CHLA_C],type="p",pch=1, col=1,cex=2,add=TRUE)

	matplot(CHLA_CHLA[SPIKE_CHLA_A],PRES_CHLA[SPIKE_CHLA_A],type="p",pch=1, col=2,cex=3,add=TRUE)
	
	legend("bottomright",c("Chl-a","Spike_C","Spike_A"),pch=c(20,20,20),col=c(8,1,2))

	dev.off()


	jpeg(file=path_out_quenchingjpeg)

	matplot(CHLA_CHLA,PRES_CHLA,col=8,type="l",ylab="Depth [m]",xlab=expression("Chlorophyll a [mg."*m ^ -3 * "]"),xlim=c(-0.2,max(CHLA_CHLA)+0.5),ylim=rev(c(0, MLD)))

	matplot(CHLA_CHLA_NPQ_C_1,PRES_CHLA,type="l",pch=1, col=1,cex=3,add=TRUE)

	matplot(CHLA_CHLA_NPQ_C_2,PRES_CHLA,type="l",pch=1, col=5,cex=3,add=TRUE)

	matplot(CHLA_CHLA_NPQ_A,PRES_CHLA,type="l",pch=1, col=2,cex=3,add=TRUE)

	
	legend("bottomright",c("Chl-a","NPQ_C_1","NPQ_C_2","NPQ_A"),pch=c(20,20,20,20),col=c(8,1,5,2))

	dev.off()


}




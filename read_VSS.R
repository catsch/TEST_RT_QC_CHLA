read_VSS <- function(filenc,iprof_param){

### Read the pressure in the Bfile 
	PRES=ncvar_get(filenc,"PRES")

### Only work on the studied profile
	PRES_IPROF=PRES[,i_prof_param]

### Init resolution profile
	RESOLUTION_IPROF=rep(NA,length(PRES_IPROF))

### Read the VSS in the B file
	VSS=ncvar_get(filenc,"VERTICAL_SAMPLING_SCHEME")

### Only work on the VSS of the studied profile
	VSS_IPROF=VSS[i_prof_param]

### Separator between every sampling ";"
	VSS_PART=str_split(VSS_IPROF,";")

### Only keep the numbers in the character string and remove blank with str_squish
	sampling=str_split(str_squish(gsub("\\[|\\]|\\.|\\,","",gsub("[a-zA-Z()_= :-]"," ",VSS_PART[[1]][])))," ")
	
### How many sampling zones are in the VSS
	nb_zone=length(VSS_PART[[1]])

### Initialize the sampling information time is the period, max and min are relative to the pressure
	sampling_time=rep(NA,nb_zone)
	sampling_max=rep(NA,nb_zone)
	sampling_min=rep(NA,nb_zone)

### Assigned the sampling information in the correct sampling zone
	for ( i in seq(1,nb_zone) ) {

		sampling_time[i]=as.numeric(sampling[[i]][1])
		sampling_max[i]=as.numeric(sampling[[i]][2])
		sampling_min[i]=as.numeric(sampling[[i]][3])

	} 

### Assigned the sampling period*float ascent speed (10cm/s=0.1dbar/s) in a resolution profile
	for ( i in seq(1,nb_zone)) {

		RESOLUTION_IPROF[which( (PRES_IPROF <= sampling_max[i]) & (PRES_IPROF > sampling_min[i]))]=sampling_time[i]*0.1

	}


### 	Fill the NA with the min of the resolution => pres <=1dbar 

	RESOLUTION_IPROF[is.na(RESOLUTION_IPROF)]=min(sampling_time)*0.1

#	DIFF_PRES_CHLA=PRES_CHLA
#	for ( i in seq(2,length(PRES_CHLA)-1 )) {
#		DIFF_PRES_CHLA[i]=PRES_CHLA[i+1]-PRES_CHLA[i]
#	} 

	return(RESOLUTION_IPROF)

}

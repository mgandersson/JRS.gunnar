###Ideas

##### There are four scripts,user can chose to run all four or a particular script.  Note that each of the four
##### scripts produces an .rda output in the associated /Output file.

####
All1    <- TRUE
Script1 <- FALSE
Script2 <- FALSE
Script3 <- FALSE
Script4 <- FALSE

 ##### Set working directory
  #setwd("/home/jpp/Desktop/Joint Risk Model_2/Model syndr surv/")
  setwd("/home/jpp/Desktop/NowherelandJRS/")  #User Input


################################################################################
################################################################################
################################################################################

if(All1 == TRUE || Script1 == TRUE) {
    if(All1 == TRUE)  print("Running all four scripts") else print("Running ONLY script one")

    
########################################
# Script 1
########################################
  
  #All necessary???
  require(plyr)
  require(stringr)
  require(ggmap)
  require(mapproj)
  require(sp)
  require(rgdal)
  require(rgeos)
  require(maptools)
  require(taRifx.geo)
  require(fields)
  require(zoo) 
  require(MASS)
  library(coda)
  library(gam)
  library(vcd)
  library(VGAM)
  library(pscl)
  library(mvtnorm)
    
  #############################################################################
  # DATASETS 
  #############################################################################      
  #User enters name of projects, names must be consistent with other files  
  projects        = list()
    projects[[1]]   = "Syndrome1"        #User input
    #projects[[2]]   = "Respiratory"  #User input
     
  #For each project make sure file exists and is properly named
  projectsCountData    =  list()  
  for(i in 1:length(projects)){
    a1 <- paste(projects[[i]],"CountData", sep="")
    if(file.exists(paste("Data/ProjectsSyndroms/", a1, ".csv", sep="")) == TRUE){
        projectsCountData[[i]]   = paste("Data/ProjectsSyndroms/", a1, ".csv", sep="")
      } else print(paste("File ", a1, " not found in file Data/ProjectsSyndroms/", sep=""))  
  }     
     
  shapeFileName   ="Data/ShapeFile"
  shapeFileLayer  ="region"    #User input
    if(file.exists(paste(shapeFileName,"/",shapeFileLayer,".shp", sep=""))==FALSE){
        print(paste(shapeFileName,"/",shapeFileLayer,".shp", " doesn't exits", sep=""))}
    
  gridFileName    ="Data/GridFile"
  gridFileLayer   ="grid"    #User input
      if(file.exists(paste(gridFileName,"/",gridFileLayer,".shp", sep=""))==FALSE){
        print(paste(gridFileName,"/",gridFileLayer,".shp", " doesn't exits", sep=""))}
    
  #Output files used in other scripts temporarily saved to file Output
  FilenameGridInfo_Script1       = paste("Output/", shapeFileLayer , "GridInfo_Script1.rda", sep="")  
  FilenameGridInfoCount_Script1  = paste("Output/", shapeFileLayer , "GridInfoCount_Script1.rda", sep="")
  
  #############################################################################
  # SETTINGS and PARAMETERS
  #############################################################################
  
  gridtype           = "VICEgrid"  #User input???
  my.radius.lambert  = 50000       #User input (before 1000)
  
  CRSnum             = 4326        #User input 
  CRSnumSourceGrid   = 4326        #User input 4326 for IA  #3035
  CRSnumPointBuf     = 2154        #User input 
  
  AnimalofInterest     = "horse"   #User input
  CommunityIdentifier  = "area"    #Generic
  
  run.big.loop       = TRUE    #If TRUE entire model runs (takes time)
  save.part1.data    = TRUE    #If TRUE, overwrites saved data
  load.part1.data    = TRUE    #Load existing file?
  save.gridinfoCount = TRUE    #gridinfo + count data per grid
  
  population.source <- "popfile" # "popfile" or "shapefile" 
  
  #############################################################################
  # Read GIS files and adjustment coordination system
  #############################################################################
     #These are the communities:
     x.try1 <- try(readOGR(dsn=shapeFileName,layer=shapeFileLayer))
         if('try-error' %in% class(x.try1))
             {print(paste("Can't open", shapeFileName,"/", shapeFileLayer, sep=""))} else {
               shapeFile <- readOGR(dsn=shapeFileName,layer=shapeFileLayer)}

       ###Added following:
      # Add line if using unprojected shapefile
      # proj4string(shapeFile)=CRS("+init=epsg:3035")
       ### End Add
       shapeFile_WGS84 = spTransform(shapeFile, CRS(paste("+init=epsg:",CRSnum,sep="")))
       
       shapeFile_ID <- sapply(slot(shapeFile_WGS84, "polygons"), function(x) slot(x, "ID")) 
  
     #This produces a grid of county associated animal and geographic data: 
      x.try2 <- try(readOGR(dsn=gridFileName,layer=gridFileLayer))
         if('try-error' %in% class(x.try2))
             {print(paste("Can't open", gridFileName,"/", gridFileLayer, sep=""))} else {
               gridShapeFile <- readOGR(dsn=gridFileName,layer=gridFileLayer)}

       proj4string(gridShapeFile)=CRS(paste("+init=epsg:",CRSnumSourceGrid,sep=""))
       gridShapeFile_WGS84 = spTransform(gridShapeFile, CRS(paste("+init=epsg:",CRSnum,sep="")))
       coordinates_country <- coordinates(gridShapeFile_WGS84)
       grid_country_ID <- sapply(slot(gridShapeFile_WGS84, "polygons"), function(x) slot(x, "ID")) 
  
     SelectedRegion_shape <- shapeFile_WGS84
     names(SelectedRegion_shape)[1] <- "DISTRICT"                                 
  
     #Sets row num to 0 and Slot ID = 0
      empty.data2 <- as.data.frame(cbind(as.character(rep(1,length(SelectedRegion_shape))),
                                       sapply(slot(SelectedRegion_shape, "polygons"), function(x) slot(x, "ID"))))    
      row.names(empty.data2)<-row.names(SelectedRegion_shape)
      SelectedRegionSpP.frame <- SpatialPolygonsDataFrame(SelectedRegion_shape,empty.data2,match.ID = TRUE)
      proj4string(SelectedRegionSpP.frame)=CRS(paste("+init=epsg:",CRSnum,sep=""))
  
      #obtain host data from file
      if( population.source == "popfile"){
       rene.grid.definition.file <-"DATA/hosts/hosts.csv"
       rene.grid <- read.table(rene.grid.definition.file,sep=",",header=TRUE)
        hosts <- rene.grid$Value
      }

      #The plot is required to add the circles created in the next part
      plot(SelectedRegionSpP.frame, pbg="white")
  
  #############################################################################
  # FUNCTIONS
  #############################################################################
  my.pointinfo <- function(my.circle.x,my.circle.y,my.radius.lambert){
  
    #my.radius.lambert<- radius.lambert
    #my.circle.x <-my.grid[k,2]  
    #my.circle.y <-my.grid[k,3]
    
    #These points of interest are gridded data containing e.g. number of horses  
    point.of.interest <- SpatialPoints(data.frame(x = my.circle.x, y = my.circle.y))
      proj4string(point.of.interest)=CRS(paste("+init=epsg:",CRSnum,sep=""))
  
    #A point is the center of a buffer of overlapping communues
    my.point.lambert<- spTransform( point.of.interest, CRS(paste("+init=epsg:",CRSnumPointBuf,sep="")) )
      my.buf.lambert <- gBuffer(my.point.lambert, width = my.radius.lambert, byid=TRUE,quadsegs = 50)
      my.buf.reprojected <- spTransform( my.buf.lambert , CRS(paste("+init=epsg:",CRSnum,sep="")) )
  
    #check which communes are covered by circle (SLOOOWWWWWWWW)
    communes <- over( SelectedRegionSpP.frame,my.buf.reprojected )
     communes[which(!is.na(communes))]
  
    #select the data associated with each cirle
    #get population from shapefile or separate file
    if(population.source== "shapefile"){
    animals.in.neighbors.communes<-SelectedRegion_shape[,AnimalofInterest][[1]][(which(communes==1))] 
       totalanimals <- sum(animals.in.neighbors.communes,na.rm=TRUE)
    } else{
      totalanimals<- sum(hosts[which(communes==1)])
    }
    
       plot(gBoundary(my.buf.reprojected),add=TRUE, col=2)
      
    #numbers of communes?  
    neighbor.communes.CommunityIdentifier <-SelectedRegion_shape[,"DISTRICT"][[1]][(which(communes==1))]
                                                                            
    #get area from shapefile in lambert etendu (THIS IS ONE COMMUNE i.e., NOT neighbors)
       #area of each communes neighbor  
    area.of.neighbors.communes<-sapply(slot(shapeFile, "polygons"),
                                     function(x) slot(x, "area"))[(which(communes==1))]
       #commune polygon in point of interest
       in.commune <-over(SelectedRegionSpP.frame,point.of.interest)
       in.commune.CommunityIdentifier <- SelectedRegion_shape[,"DISTRICT"][[1]][which(in.commune==1)]
    
    #store information
    mylist<-list()
    mylist[["within_Community"]]<- in.commune.CommunityIdentifier
    mylist[["animals.in.range"]]<-totalanimals
    mylist[["neighbors.Communities"]]<-neighbor.communes.CommunityIdentifier
    mylist[["neighbors.area"]]<-area.of.neighbors.communes
    mylist[["total.neighbors.area.km3"]]<-sum(area.of.neighbors.communes,na.rm =TRUE)/1000000
    mylist[["animals.per.km3"]]<-1000000*totalanimals/sum(area.of.neighbors.communes,na.rm=TRUE)
    
    return(mylist)
  }#end function
    
  #Could be another grid type, now just VICEgrid
  if (gridtype =="VICEgrid"){
    my.gridindex    <- as.numeric(grid_country_ID)
    my.xval         <- as.numeric(coordinates_country[,1])
    my.yval         <- as.numeric(coordinates_country[,2])
    my.grid         <- as.data.frame(cbind(my.gridindex,my.xval,my.yval))
    my.gridinfo     <- as.list(rep(NA,length(my.gridindex)))
    radius.lambert  <- my.radius.lambert
  }
                                                                  
  #############################################################################
  # This calls the main loop
  #############################################################################  
  if(run.big.loop == TRUE){  
    for(k in 1:length(my.gridindex)){
      print(100*k/length(my.gridindex))
      my.gridinfo[[k]]  <- my.pointinfo(my.grid[k,2],my.grid[k,3],radius.lambert)
    }
  } else {print("Run big loop is set to FALSE, old data will be loaded")}
  
  
  #############################################################################
  # EXPORT/LOAD my.gridinfo
  #############################################################################
                                                                  
  if(save.part1.data == TRUE){
    part1.data.grid <- list(my.grid,my.gridinfo)
    save(part1.data.grid, file = FilenameGridInfo_Script1)
  }
                                                                                                                                 
  if(load.part1.data == TRUE){
    load(FilenameGridInfo_Script1)
    my.grid <- part1.data.grid[[1]]
    my.gridinfo <- part1.data.grid[[2]]
  }
  
  
  ############################################################################
  #############################################################################
  # Read COUNT data
  #############################################################################
  ############################################################################
  
  for(j in 1:length(projects)){
    
    project <- projects[[j]]
    raw.data <- projectsCountData[[j]]

    x.try3  <- read.table(raw.data,sep=";",header=TRUE,allowEscapes = TRUE, fill = TRUE,)
        if('try-error' %in% class(x.try3))
            {print(paste("Can't open ", raw.data))} else {           
              count.data.points <- read.table(raw.data,sep=";",header=TRUE,allowEscapes = TRUE, fill = TRUE,)}
    
    commune.level.counts <- as.data.frame(cbind(as.character(count.data.points$date_decl),
                                                as.numeric(as.character(count.data.points$N_Community))))
      
    names(commune.level.counts)<-c("date","area1")
  
    ##################################################
    ############ prepare count- vecors for grid points.
    ##################################################
    totcounts.test <- 0; maxtest <- 0
    for(i in 1:length(my.gridinfo)){
      my.neighbors.in.area <- na.omit(my.gridinfo[[i]]$neighbors.Communities)
      my.commune.area   <- my.gridinfo[[i]]$"within_Community"
      
      if(!length(my.neighbors.in.area)==0){
          varname1 <- paste("base", project, "Counts", sep="")  
          temp1 <- as.character(commune.level.counts$date[which(!is.na(match(commune.level.counts$area1, my.neighbors.in.area)))])
          my.gridinfo[[i]][[varname1]] <- temp1
          totcounts.test <- totcounts.test + length(my.gridinfo[[i]][[varname1]])
          if(length(my.gridinfo[[i]][[varname1]])>maxtest){
            maxtest <- length(my.gridinfo[[i]][[varname1]])
          }
      }
    }#end of  for i in length my gridinfo...
    

  }# end for each project
  
  if(save.gridinfoCount == TRUE){
     save(my.gridinfo,file = paste("FilenameGridInfoCount", shapeFileLayer,  "_Script1", sep=""))
  }
  
} #END script 1 loop


########################################
# END Script 1
########################################    

########################################
# Script 2
########################################

if(All1 == TRUE || Script2 == TRUE) {
    if(All1 == TRUE)  print("Running all four scripts, in Script 2") else print("Running ONLY script two")        
  
  dataOutbreaks <-  "Data/Outbreaks/"
  outbreakfilename <- "Syndrome1_outbreakArea_A"
  
  #Names of outbreak file(s):
  numberOutbreakGeographicRegions <- 1  #User input
  namesArea <- list()  
  for(i in 1:numberOutbreakGeographicRegions){
      name1 <- paste(outbreakfilename, i, sep="")
      if(file.exists(paste(dataOutbreaks, name1, ".csv", sep="")) == TRUE){
        print(paste("File", paste(dataOutbreaks, name1, ".csv", sep="")," found", sep=""))
        namesArea[[i]] <- name1
        }else print(paste("File ", name1, " not found in file Data/Outbreaks/", sep=""))
   }    
        
  #Check if syndroms/projects from above in both files
  numberofSyndroms      <- projects
  dataOutbreaksFiles    <- list(); projNames <- list()  
  for(i in 1:numberOutbreakGeographicRegions){
      name1 <- paste(outbreakfilename, i, sep="")
      x.try4 <- read.table(paste(paste(dataOutbreaks, name1, ".csv", sep="")), sep=";", header=T)
      if('try-error' %in% class(x.try4))
             {print(paste("Can't open", paste(dataOutbreaks, name1, ".csv", sep="")), sep="")
        } else {
                dataOutbreaksFiles[[i]]  <- read.table(paste(paste(dataOutbreaks, name1, ".csv", sep="")), sep=";", header=T)
                checknames <- names( dataOutbreaksFiles[[i]])
                projNames[[i]]  <- paste("decl_", projects, sep="")
                if(all(as.character(unlist(projNames[[i]])) %in% checknames) == FALSE) {print(paste("Not all syndroms in file ", paste(dataOutbreaks, name1, ".csv", sep="")))}
          }
   }
  
  #Input from previous script is input for this script 
  baseData  <- paste("FilenameGridInfoCount", shapeFileLayer, "_Script1", sep="")
    
  #Output file, where data from this script is eventually stored   
  outFile   <- paste("Output/",shapeFileLayer, "GridDataCountOutbreaks_Script2.rda", sep="")
  
  #Time series data
  if(file.exists("Data/RefTimelineVice/ref_timeline_VICE.csv") == FALSE) {
                  print(paste("File: Data/RefTimelineVice/ref_timeline_VICE.csv doesn't exist"))
              } else {timeData <- "Data/RefTimelineVice/ref_timeline_VICE.csv"}
  
  ##################################################
  # Code
  ##################################################
  
  #Read in data; load my.gridinfo
  load(baseData)
  
  #Loop over "i" outbreak areas
  for(i in 1:numberOutbreakGeographicRegions){
      
    outbreak.raw  <- dataOutbreaksFiles[[i]]
    outbreak.year <- 2012
  
    #Time series data
    timeLine <- read.table(timeData,sep=";",header=TRUE,allowEscapes = TRUE, fill = TRUE,)  
    outbreak.raw$absweek <- outbreak.raw$week_year_num + (timeLine$week_abs[which(timeLine$Date==paste(outbreak.year,"-01-01",sep=""))]-1)
   #> outbreak.raw[which(outbreak.raw$decl_Syndrome1>0),] # to check

    #####ADDED jpp Aug 13th ########## Need to add DATE in good form to outbreak.raw
       #Get right year
       timeln1  <- subset(timeLine, timeLine$year == outbreak.year)
       #Get right weeks in year
       timeln2 <-  subset(timeln1, timeln1$week_year_num %in% outbreak.raw$week_year_num)[,2:10]
          ###Get first occurrence of event in a given week
          firstDayWeek <- do.call(rbind,(lapply(split(timeln2, timeln2$week_year_num), function(x) x <- x[1,])))  
       outbreak.raw <- merge(outbreak.raw, firstDayWeek[,c("week_year_num", "Date")], all.x=T, by="week_year_num")
    #####END: ADDED jpp Aug 13th #####
    
      #Loop over "j" grids
      for(j in 1:length(my.gridinfo)){          
         withinCommunity     <- my.gridinfo[[j]]$within_Community  
         neighborCommunities <- na.omit(my.gridinfo[[j]]$neighbors.Communities)
              ### TEST neighborCommunities  <- c(neighborCommunities, 31550)
         
         #Loop over "g" number of projects per geographic area
         for(g in 1:length(projects)){
           ####Set = 0      
           my.gridinfo[[j]][[paste(namesArea[[i]],projNames[[i]][g],sep="_")]]  <-rep(0,timeLine$week_abs[length(timeLine$week_abs)])
          # my.gridinfo[[j]][[paste(namesArea[[i]],projNames[[i]][g],"Date",sep="_")]] <-rep(0,timeLine$week_abs[length(timeLine$week_abs)])
           #########   
              temp <- outbreak.raw[which(!is.na(match(outbreak.raw$Community,neighborCommunities)) &outbreak.raw[[projNames[[i]][g]]]>0),]
           if(length(temp[,1])>0){
             #print(paste(i=i, j=j,g=g))
             #Loop over q areas with outbreaks 
             for(q1 in 1:length(temp[,1])){
                my.gridinfo[[j]][[paste(namesArea[[i]],projNames[[i]][g],sep="_")]][temp$absweek[q1]] <- temp[[projNames[[i]][g]]][[q1]]+ my.gridinfo[[j]][[paste(namesArea[[i]],projNames[[i]][g],sep="_")]][temp$absweek[q1]]                
               # my.gridinfo[[j]][[paste(namesArea[[i]],projNames[[i]][g],"_Date",sep="")]][temp$absweek[q1]]  <- as.Date(temp$Date[[q1]])
                
             }
           }#end if 
         }#end syndroms  
       }#end grids
    }#end outbreak areas
  
  
  save(my.gridinfo,file=outFile)
   
} #End Script 2     

########################################
# End Script 2
########################################

########################################
# Script 3
########################################

if(All1 == TRUE || Script3 == TRUE) {
    if(All1 == TRUE)  print("Running all four scripts, in Script 3") else print("Running ONLY script three")  
  
  #Load functions
  x.try5 <- try(source("Functions/surveillanceFunctions.r"))
      if('try-error' %in% class(x.try5))
             {print(paste("Can't open: Functions/surveillanceFunctions.r, sep="))} else {
               source("Functions/surveillanceFunctions.r")}
    
  #These are pointers to data files
  timeData <-"Data/RefTimelineVice/ref_timeline_VICE.csv"
    
  if(file.exists(paste("Data/SurveillanceData/", AnimalofInterest, ".csv", sep="")) == TRUE) {
     animalsPerDistrict   <- paste("Data/SurveillanceData/", AnimalofInterest, ".csv", sep="")
     } else print(paste("Can't find: Data/SurveillanceData/", AnimalofInterest, ".csv", sep=""))

  for(i in 1:length(projects)){
    name1 <- paste("SurveillanceData", projects[[i]], sep="_")  
    if(file.exists(paste("Data/SurveillanceData/", name1, ".csv", sep="")) == TRUE) {
     name2 = paste("raw.data.", projects[[i]], sep="")
     name3 <- paste("Data/SurveillanceData/", name1, ".csv", sep="")     
     delayedAssign(name2, name3)
     } else print(paste("Can't find Data/SurveillanceData/", name1, ".csv", sep=""))
   }

  #If vStuff TRUE, save it for next script  
  vStuff        <-  TRUE  
  
  #These are user choices LARGE (vs small) geographic areas
  outbreak.IDs  <- namesArea
  
  #Choice for calculating the value of evidence
  outbreaktype  <- "Geom.from.pop"   

  #Read in Negbinom Parameters
  if(file.exists("Data/SurveillanceData/NegbinomParameters.csv") == TRUE){
       par1 <- read.csv("Data/SurveillanceData/NegbinomParameters.csv", header=T, sep=";")
     } else {print("File: Data/SurveillanceData/NegbinomParameters.csv not found")}
    
  if(!(any(projects %in% par1$Syndrom)) == TRUE) print(paste("Project", projects, "data not found."))

  for(i in 1:length(projects)){      
       binomOutbreakData <- subset(par1, par1$Syndrom == projects[[i]])
         name1  <- paste("outbreak.mean.", projects[[i]], sep="")
         name2  <- subset(binomOutbreakData, binomOutbreakData$Variable == "Mean")$Data
       delayedAssign(name1, name2)
         name3  <- paste("outbreak.theta.", projects[[i]], sep="")
         name4  <- subset(binomOutbreakData, binomOutbreakData$Variable == "Theta")$Data
       delayedAssign(name3, name4)
   }
    
  #Input from previous script and  the eventual output file
    if(file.exists(paste("Output/GridDataCountOutbreaks", shapeFileLayer, "_Script2.rda", sep=""))==TRUE){
        baseData2  <- paste("Output/GridDataCountOutbreaks", shapeFileLayer, "_Script2.rda", sep="")
    }else {paste("File: Output/GridDataCountOutbreaks", shapeFileLayer, "_Script2.rda not found", sep="")}
    
  outFile    <- paste("Output/", shapeFileLayer, "GridDataPlusOutbreaksSurv_Script3.rda", sep="")
  
  #load base data from previous script
  x.try6 <- try(load(baseData2))
    if('try-error' %in% class(x.try6))
            {print("Can't open baseData2")} else {load(baseData2)} 
  
  #############################################################
  ######### start program #####################################
  #############################################################
  
  #Set up for next section
  for(i in 1:length(projects)){
     name1 <- paste("all.counts.", projects[[i]],sep="") 
     name2 <- paste("all.week.counts.", projects[[i]],sep="")
     name3 <- paste("all.week.time.series.",  projects[[i]],sep="")
     name4 <- paste("all.v.",  projects[[i]],sep="")
     delayedAssign(name1, list())
     delayedAssign(name2, list())
     delayedAssign(name3, list())
     delayedAssign(name4, list())}
    
  for(i in 1:length(projects)){
     name1 <- paste("maxcounts.", projects[[i]],sep="") 
     name2 <- paste("maxweekcounts.", projects[[i]],sep="")
     name3 <- paste("meanweekcounts.",  projects[[i]],sep="")
     name4 <- paste("allcounts.",  projects[[i]],sep="")
     delayedAssign(name1, c())
     delayedAssign(name2, c())
     delayedAssign(name3, c())
     delayedAssign(name4, c())}
    
   gridpop <- rep(NA,length(my.gridinfo))
    
  ##### Multiple loops
  #Loop over grids  
  for(i in 1:length(my.gridinfo)){
    #Short loop over projects  
    for(j in 1:length(projects)){
        name1 <- paste("all.counts.", projects, sep="")
      delayedAssign(name1[[j]], NA)
        name2 <- paste("maxcounts.", projects, sep="")
      delayedAssign(name2[[j]], NA)
    }  

    my.gridinfo[[i]][["out.geom.p"]] <- my.out.params.from.pop(my.gridinfo[[i]][["animals.in.range"]])

    #Short loop over projects
    for(q in 1:length(projects)){
        varName1 <- paste("base", projects[[q]], "Counts", sep="")    
      if(length(my.gridinfo[[i]][[varName1]])!=0){
          outDates1 <- subset(my.gridinfo[[i]][[varName1]], my.gridinfo[[i]][[varName1]] != "0")
          my.return <- make.time.series(outDates1)
          } else{my.return <- make.time.series(c())}

      if(length(my.return[[1]]) != 0){ 
      eval(parse(text=(paste("all.week.time.series.", projects[[q]], "[[", i, "]] <- ", my.return[[1]], sep=""))))}
      eval(parse(text=(paste("all.counts.", projects[[q]], "[[", i, "]] <- ", my.return[[2]], sep=""))))
              temp1 <- eval(parse(text=(paste("all.counts.", projects[[q]], "[[", i, "]] <- ", my.return[[2]], sep=""))))
      eval(parse(text=(paste("all.week.counts.", projects[[q]], "[[", i, "]] <- ", my.return[[3]], sep=""))))
              temp2 <- eval(parse(text=(paste("all.week.counts.", projects[[q]], "[[", i, "]] <- ", my.return[[3]], sep=""))))
              temp3 <- eval(parse(text=(paste("allcounts.", projects[[q]], sep=""))))
      eval(parse(text=(paste("allcounts.", projects[[q]], "[[", i, "]] <- ", c(temp3 , my.return[[3]]), sep=""))))
      eval(parse(text=(paste("maxcounts.", projects[[q]], "[[", i, "]] <- ", max(temp1), sep=""))))
      eval(parse(text=(paste("maxweekcounts.", projects[[q]], "[[", i, "]] <- ", max(temp2), sep=""))))
      eval(parse(text=(paste("meanweekcounts.", projects[[q]], "[[", i, "]] <- ", max(temp2), sep=""))))          

      gridpop[[i]] <- my.gridinfo[[i]]$animals.in.range    
    }
  } #End gridinfo Loop p
    print("Check:  Many things calculated, but only a few used?")

  #Load surveillance data, we know from above that files exists
   firsttwoYears <- list(); histCounts <- list(); all.baseline.counts.X <- list() 
   for(i in 1:length(projects)){
      projx <- read.csv(paste("Data/SurveillanceData/SurveillanceData_", projects[[i]],".csv", sep=""), sep=";")
      all.baseline.counts.X[[i]] <- make.time.series2(projx$date_decl, t.start="1.01.2006",t.end="31.12.2013",level="week")[[3]]
      tx <- seq(1:length(all.baseline.counts.X[[i]]))
      week.of.year <- timeLine$week_year_num[match(tx,timeLine$week_abs)]
      t.year <- timeLine$year[match(tx,timeLine$week_abs)]
      firsttwoYears[[i]] <- all.baseline.counts.X[[i]][1:106]  #Needs to be generalized
      histCounts[[i]]    <- c(firsttwoYears[[i]],all.baseline.counts.X[[i]])
    }
    
  ### prepare model inputs
  histMean.X <- list()  
  for(i in 1:length(histCounts)){  
  #histmean is based on the total data and total population
     q1 <- list() 
     for(r in 1:length(tx)){
      q1[[r]] <- mean(histCounts[[i]][(r+43):(r+96)])  #Needs to be generalized 
     }
   histMean.X[[i]] <- q1
   }   
    
  loghistMean.X <- list()  
  for(i in 1:length(histMean.X)){
     loghistMean.X[[i]] <-log(unlist(histMean.X[[i]]))
   }  
  
  #total animal population 
  totpop  <- sum(read.table(animalsPerDistrict, sep=";",header=TRUE,allowEscapes = TRUE, fill = TRUE,)$horse  )   

  model.X <- list()  
  for(i in 1:length(all.baseline.counts.X)){
      model.X[[i]] <- make.fit(all.baseline.counts.X[[i]], loghistMean.X[[i]], tx, week.of.year, log(totpop),time.unit = "week", distribution="poisson")
  }

  all.grid.counts.X <- list()  
  for(j in 1:length(projects)){
      all.grid.counts <- c()
      for(i in 1:length(my.gridinfo)){
          all.grid.counts  <- c(all.grid.counts, unlist(eval(parse(text = paste("all.week.counts.", projects[[j]], sep="")))))
      }
      all.grid.counts.X[[j]] <- as.vector(unlist(all.grid.counts))
  }

  #Why just for Neuro?
  all.grid.logpop <- c()  
  for(i in 1:length(my.gridinfo)){  
     all.grid.logpop   <-c(all.grid.logpop,rep(log(gridpop[[i]]),length(eval(parse(text=paste("all.week.counts.", projects[[1]],  sep=""))))))
  }

 
  all.grid.seasonal.X <- list()
  for(j in 1:length(projects)){
     all.grid.seasonal <- c()
     for(i in 1:length(my.gridinfo)){
         all.grid.seasonal <- c(all.grid.seasonal,  model.X[[j]][[2]]$fit)
     }
     all.grid.seasonal.X[[j]] <- all.grid.seasonal
  }
     
  ###
  all.grid.t <- c(); all.grid.week.of.year <- c(); all.grid.gridindex <- c()  
  for(i in 1:length(my.gridinfo)){  
    all.grid.t <- c(all.grid.t, tx)
    all.grid.week.of.year <- c(all.grid.week.of.year,week.of.year) 
    all.grid.gridindex <- c(all.grid.gridindex,rep(i, length(tx)))
  }

  ###  NOT SOLVING all.grid.seasonal.X is wrong length
  grid.fit.X <- list(); all.grid.exp.mean.X <- list()  
  for(j in 1:length(projects)){
      grid.fit.X[[j]]  <- glm(all.grid.counts.X[[j]]  ~ all.grid.logpop +vlog(all.grid.seasonal.X[[j]]), family="poisson")
      all.grid.exp.mean.X[[j]] <- predict.glm(grid.fit.X[[j]] ,type="response", se.fit=TRUE,newdata=as.data.frame(all.grid.seasonal.X[[j]]))
  }
  
  ###################################################################
  ###### Fit grid counts to zero truncated geometric distribution ###
  ###################################################################
  #prepare truncated vectors
  trunk.counts.X <- list()  
  for(j in 1:length(projects)){  
    trunk.counts.X[[j]] <- all.grid.counts.X[[j]][which(all.grid.counts.X[[j]] != 0)]
  }
    
  trunk.seasonal.X <- list()  
  for(j in 1:length(projects)){  
    trunk.seasonal.X[[j]] <- all.grid.seasonal.X[[j]][which(all.grid.seasonal.X[[j]] != 0)]
  }

  trunk.exp.mean.X <- list()  
  for(j in 1:length(projects)){  
    trunk.exp.mean.X[[j]] <- all.grid.exp.mean.X[[j]][which(all.grid.X[[j]] != 0)]
  }
    
  # fit truncated geometrical distribution to non-zero counts and corresponding expected mean of gridpoint
  
  trunk.model.geom.X <- list(); geom.coeff.X <- list()
  for(j in 1:length(projects)){  
    trunk.model.geom.X[[j]] <- vglm(trunk.counts.X[[j]] ~ trunk.exp.mean.X[[j]], link=log,family = truncgeometric(100))
    geom.coeff.X[[j]]       <- coef(trunk.model.geom.X[[j]], matrix = TRUE)
  }
  
  #extract vector with p parameter for zero trunkated geometrical distribution.
  geom.param.for.trunk.X <- list(); all.grid.geom.param.X <- list()
  for(j in 1:length(projects)){    
     geom.param.for.trunk.X[[j]] <- exp(as.numeric(geom.coeff.X[[j]][1])+ as.numeric(geom.coeff.X[[j]][2])* trunk.exp.mean.X[[j]])
     all.grid.geom.param.X[[j]]  <- exp(as.numeric(geom.coeff.X[[j]][1])+ as.numeric(geom.coeff.X[[j]][2])* all.grid.exp.mean.X[[j]]$fit)
  }

  all.grid.geom.mean.X <- list();  all.grid.geom.nonzero.X <- list()  
  for(j in 1:length(projects)){
     geom.param.for.trunk.X[[j]]  <- trunk.mean(geom.param.for.trunk.X[[j]])
     all.grid.geom.mean.X[[j]]    <- trunk.mean(all.grid.geom.param.X[[j]])
     all.grid.geom.nonzero.X[[j]] <- all.grid.exp.mean.X[[j]]$fit/all.grid.geom.mean.X[[j]]

    #insert in my gridinfo
    for( i in 1:length(my.gridinfo)){
        for(j in 1:length(projects)){
           name1 <- paste(projects[[j]], ".geom.p", sep="")
           name2 <- paste(projects[[j]], ".nonzero", sep="")
           my.gridinfo[[i]][[name1]]  <- all.grid.geom.param.X[[j]][which(all.grid.gridindex == i)]
           my.gridinfo[[i]][[name2]]  <- all.grid.geom.nonzero.X[[j]][which(all.grid.gridindex == i)]
         } 
     }
  }
  #################################################
  ############# make plots for geom distributions
  #################################################
  
  ## plot(c(),c(),xlim=c(0,0.15),ylim=c(0,2),xlab="grid mean",ylab="")
  ## points(all.grid.exp.mean.neuro$fit,all.grid.neuro.geom.mean,pch=".")
  ## mtext("mean of trunk geom",side=2,padj=-3)
  ## par(new = T)
  ## plot(all.grid.exp.mean.neuro$fit,all.grid.neuro.geom.nonzero,pch=".",xlim=c(0,0.15),axes=F,col=2,xlab="",ylab="")
  ## axis(side=4,padj=-1.4,cex.axis=0.8,tck=-0.01,col=2)
  ## mtext("p nonzero",side=4,padj=1.5,col=2)
  
  ## #plot(seq(1,10),log10(trunk.geom(seq(1,10),0.66)))
  ## plot(c(),c(),xlim=c(0,2.5),ylim=c(0,3),xlab="grid mean",ylab="")
  ## points(all.grid.exp.mean.respiratory$fit,all.grid.resp.geom.mean,pch=".")
  ## mtext("mean of trunk geom",side=2,padj=-3)
  ## par(new = T)
  ## plot(all.grid.exp.mean.respiratory$fit,all.grid.resp.geom.nonzero,pch=".",xlim=c(0,2.5),axes=F,col=2,xlab="",ylab="")
  ## axis(side=4,padj=-1.4,cex.axis=0.8,tck=-0.01,col=2)
  ## mtext("p nonzero",side=4,padj=1.5,col=2)
  
  ## i<-173
  ## nmax<-40
  ## test<-my.calculate.v(i,nmax)
  
  ## #calculate.v <- FALSE#  #"base", "outbreak.ID"
  ##### This is slow
  v.to.make <- c(outbreak.IDs,"base")
  
  if(vStuff == TRUE){
   for(i2 in 1:length(v.to.make)){
     outbreak.ID <- v.to.make[i2] 
     if(outbreak.ID=="base"){
       for(j in 1:length(projects)){  
         for(i in 1:length(my.gridinfo)){
           nmax<-40
           my.return<- my.calculate.v(i,nmax)
           my.gridinfo[[i]][[paste("v.", projects[[j]], sep="")]] <- my.return[[paste("v.", projects[[j]], sep="")]]
           my.gridinfo[[i]][["v.tot"]]   <-my.return[["v.tot"]]
         }
       }
     }else{
          for(j in 1:length(projects)){  
           v1[[j]] <- paste(outbreak.ID,paste("_v",projects[[j]],sep=""))
           v3 <- paste(outbreak.ID,"_vtot",sep="")
           for(i in 1:length(my.gridinfo)){
             nmax<-60
             #i<-1
             my.return <- my.calculate.v(i,nmax,outbreak.ID)
             my.gridinfo[[i]][[v1[[j]]]]<-my.return[[paste("v.", projects[[j]], sep="")]]
             my.gridinfo[[i]][[v3]]<-my.return[["v.tot"]]
            }
          }
         }
   }
   
  save(my.gridinfo, file=outFile)
  }else print("vStuff not run")

} #End script 3    
    
########################################
# End Script 3
########################################
  
########################################
# Script 4
########################################

if(All1 == TRUE || Script4 == TRUE) {
    if(All1 == TRUE)  print("Running all four scripts, in Script 4") else print("Running ONLY script four")  

  
  #rm(list=ls())
  
  ##############
  # GA comments Aug 25
  # We should not "hard code" the variable names for the outbreak scenarios..
  
  inFile     <- paste("Output/", shapeFileLayer, "GridDataPlusOutbreaksSurv_Script3.rda", sep="")
  outFile    <- paste("Output/", shapeFileLayer, "GridDataPlusOutbreaksSurvJRSFile_Script4.rda", sep="")
  
   x.try7 <-read.csv2("Data/IntroRisks/P_release2.csv")
        if('try-error' %in% class(x.try3))
            {print(paste("Can't open Data/IntroRisks/P_release2.csv"))} else {           
              release3 <- read.csv2("Data/IntroRisks/P_release2.csv")
          }
  
  #previously done
  gridShapeFile <- readOGR(dsn=gridFileName,layer=gridFileLayer)
  grid_no       <- slot(gridShapeFile,"data")$Grid_Code # to det the index in "my.gridinfo"
  gridindex     <- as.numeric(sapply(slot(gridShapeFile,"polygons"),function(x) slot(x,"ID")))+1
  coding        <- cbind(grid_no,gridindex)
  release3      <- merge(release3, coding, by=c("grid_no"))
  
  ###calculate p. release at the grid level and on a weekly scale (4 weeks per month):
  release3$Hrel50_ee  <- 1-(1-release3$Hrel50_ee)^(1/(release3$tot_grid*4))
  release3$Hrel5_ee   <- 1-(1-release3$Hrel5_ee)^(1/(release3$tot_grid*4))
  release3$Hrel95_ee  <- 1-(1-release3$Hrel95_ee)^(1/(release3$tot_grid*4))
  
  release3$Vrel50_ee  <- 1-(1-release3$Vrel50_ee)^(1/(release3$tot_grid*4))
  release3$Vrel5_ee   <- 1-(1-release3$Vrel5_ee)^(1/(release3$tot_grid*4))
  release3$Vrel95_ee  <- 1-(1-release3$Vrel95_ee)^(1/(release3$tot_grid*4))
  
  ###There are P_release=0 because no importation or no risk in low risk countries.
  #To avoid output= 0, Replace P_release=0 by "lowest":
  low0   <- subset(release3$Hrel50_ee, release3$Hrel50_ee!="0")
  lowest <- min(low0)*10/100 #"lowest" = 10% of the lowest probability
  release3$Hrel50_ee <- replace(release3$Hrel50_ee, release3$Hrel50_ee=="0", lowest)
  
  low0   <- subset(release3$Hrel5_ee, release3$Hrel5_ee!="0")
  lowest <- min(low0)*10/100 #"lowest" = 10% of the lowest probability
  release3$Hrel5_ee <- replace(release3$Hrel5_ee, release3$Hrel5_ee=="0", lowest)
  
  low0   <- subset(release3$Hrel95_ee, release3$Hrel95_ee!="0")
  lowest <- min(low0)*10/100 #"lowest" = 10% of the lowest probability
  release3$Hrel95_ee <- replace(release3$Hrel95_ee, release3$Hrel95_ee=="0", lowest)
  
  low0   <- subset(release3$Vrel50_ee, release3$Vrel50_ee!="0")
  lowest <- min(low0)*10/100 #"lowest" = 10% of the lowest probability
  release3$Vrel50_ee <- replace(release3$Vrel50_ee, release3$Vrel50_ee=="0", lowest)
  
  low0   <- subset(release3$Vrel5_ee, release3$Vrel5_ee!="0")
  lowest <- min(low0)*10/100 #"lowest" = 10% of the lowest probability
  release3$Vrel5_ee <- replace(release3$Vrel5_ee, release3$Vrel5_ee=="0", lowest)
  
  low0   <- subset(release3$Vrel95_ee, release3$Vrel95_ee!="0")
  lowest <- min(low0)*10/100 #"lowest" = 10% of the lowest probability
  release3$Vrel95_ee <- replace(release3$Vrel95_ee, release3$Vrel95_ee=="0", lowest)
  
  ###merge with timeData
  timeData               <- unique(read.csv2("Data/RefTimelineVice/ref_timeline_VICE.csv")[,c(3,5,7)])
  colnames(timeData)[1]  <- "month"
  release4<-merge(release3, timeData, by=c("month","year"))
  head(release4)
  
  load(inFile)
  #JRSfile<- paste("DATA/grid_gridinfo_JRS_france.Mar20.",radius,".rda",sep="")
  
  for(i in 1:length(my.gridinfo)){
    #i<-500
    temp <- release4[which(release4$grid_no == coding[i,1]),][order(release4[which(release4$grid_no == coding[i,1]),]$week_abs),] 
    variables <- names(temp)[6:11]
  
    for(j in 1:length(variables)){
      #j<-1
      temp1 <- rep(NA,length(my.gridinfo[[i]]$v.tot))
      temp2 <- rep(NA,length(my.gridinfo[[i]]$v.tot))
    
      for (k in 1:length(temp[,1])){
       #k<-1
        temp1[temp$week_abs[k]]  <-eval(parse(text=paste("temp$",variables[j],"[",k,"]",sep="")))
        temp2[temp$week_abs[k]] <- log10( temp1[temp$week_abs[k]]  )
      }
     eval(parse(text=paste("my.gridinfo[[",i,"]]$",variables[j],"<-temp2",sep="")))
    }
  }
  
  for(i in 1:length(my.gridinfo)){
    #i<-870  #860 has peaks 870 is Aquitaine, 140 has Normandie, 138, 139 has random peaks
    loghost.p <-eval(parse(text=paste("my.gridinfo[[",i,"]]$",variables[1],sep="")))
    logvector.p <-eval(parse(text=paste("my.gridinfo[[",i,"]]$",variables[4],sep="")))
    p.prior <- 10^loghost.p + 10^logvector.p -((10^loghost.p)*(10^logvector.p))
    logO.prior <- log10(p.prior/(1-p.prior))
    v.tot.base <-eval(parse(text=paste("my.gridinfo[[",i,"]]$","v.tot",sep="")))
    v.tot.normandie <-eval(parse(text=paste("my.gridinfo[[",i,"]]$","Normandie_vtot" ,sep="")))
    v.tot.aquitaine <-eval(parse(text=paste("my.gridinfo[[",i,"]]$","Aquitaine_vtot",sep="")))
    logJRS.base <- v.tot.base + logO.prior
    logJRS.Normandie <- v.tot.normandie + logO.prior
    logJRS.Aquitaine <- v.tot.aquitaine + logO.prior
    eval(parse(text=paste("my.gridinfo[[",i,"]]$","logO.prior","<-logO.prior",sep="")))
    eval(parse(text=paste("my.gridinfo[[",i,"]]$","logJRS.base","<-logJRS.base",sep="")))
    eval(parse(text=paste("my.gridinfo[[",i,"]]$","logJRS.Normandie","<-logJRS.Normandie",sep="")))
    eval(parse(text=paste("my.gridinfo[[",i,"]]$","logJRS.Aquitaine","<-logJRS.Aquitaine",sep="")))
    
   # plot(c(),c(),xlim=c(200,400),ylim=c(-12,4))
    #lines(seq(1:length(p.prior)),logO.prior,ylim=c(-12,4),col=1)
    #lines(seq(1:length(p.prior)),logJRS.base,col=2)
    #lines(seq(1:length(p.prior)),logJRS.Normandie+0.1,col=3)
    #lines(seq(1:length(p.prior)),logJRS.Aquitaine+0.2,col=4)
  }
  
  save(my.gridinfo , file=outFile)
}  
  

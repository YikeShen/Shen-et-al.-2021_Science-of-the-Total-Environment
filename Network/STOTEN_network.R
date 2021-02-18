#Correlation screening
# Figure 6 calculation
# Figure 6 used this R output files and input into Gephi for graphics

rm(list=ls())
setwd("/Users/yikeshen/Desktop/Salmonella_Paper")

#Bacterial Community antibiotics
OTUtableraw <- read_excel("Salmonella_OTU_Silva_Dec.xlsx", sheet = "OTU_all")
OTUtable <- OTUtableraw %>% filter(!grepl("Archaea",OTUID)) %>% filter(!grepl("Eukaryota",OTUID))
OTU <- OTUtable %>% dplyr::select(-"OTUID")
ROWNAMES <- OTU$FeatureID 
OTU = OTU[,mixedsort(colnames(OTU))]
row.names(OTU) <- OTU$FeatureID 
OTU <- OTU %>% select(-FeatureID) %>% 
  select(-c(P1:P12),-c(P25:P36),-c(S1:S12),-c(S25:S36)) %>% 
  select(-c(P16:P18),-c(P40:P42),-c(S16:S18),-c(S40:S42))
OTU <- OTU %>% select(-S21)#didn't amplify
OTU[is.na(OTU)] <- 0
OTU <- OTU %>% as.matrix()
row.names(OTU) <- ROWNAMES


map <- read_excel("Salmonella_OTU_Silva_Dec.xlsx", sheet = "ENV")
map <- map[!grepl("NO", map$Antibiotics),] %>% as.data.frame()
map <- map[-33,]#didn't amplifyS21
row.names(map) <- map$SampleID
map <- map %>% t() %>% as.data.frame()
map <- map %>% select(-c(P16:P18),-c(P40:P42),-c(S16:S18),-c(S40:S42)) %>% t() %>% as.data.frame()
SAMPLEIDROWNAMES <- map$SampleID
map <- map %>% as.matrix() %>% as.data.frame()
row.names(map) <- SAMPLEIDROWNAMES
ROANAMESLISTSAMPLESID <- row.names(map)

Taxonomy <- OTUtable %>% 
  dplyr::select("OTUID")

Taxonomy <- Taxonomy %>%
  separate("OTUID",into=c("Domain","Phylum","Class","Order","Family","Genus","Species"),sep=";")

TAX <- Taxonomy %>% as.data.frame() %>% 
  mutate(Domain=gsub("D_0__","",Domain)) %>% 
  mutate(Phylum=gsub("D_1__","",Phylum)) %>% 
  mutate(Class=gsub("D_2__","",Class)) %>% 
  mutate(Order=gsub("D_3__","",Order)) %>% 
  mutate(Family=gsub("D_4__","",Family)) %>% 
  mutate(Genus=gsub("D_5__","",Genus)) %>% 
  mutate(Species=gsub("D_6__","",Species)) %>% 
  as.matrix() %>% as.data.frame()

TAX[is.na(TAX)] <- NA
TAX <- apply(TAX, 2, function(x) gsub("^$|^ $", NA, x))
TAX <- apply(TAX, 2, function(x) gsub("__", NA, x))
TAX <- apply(TAX, 2, function(x) gsub("uncultured", NA, x))
TAX <- TAX %>% as.matrix()
row.names(TAX) <- ROWNAMES


OTU=otu_table(OTU,taxa_are_rows = TRUE)
rownames(TAX) <- row.names(OTU)
TAX=tax_table(TAX)
MAP <-  sample_data(map)
physeq = merge_phyloseq(OTU, TAX, MAP)
Tree=rtree(ntaxa(physeq),rooted=TRUE,tip.label = taxa_names(physeq))

physeq = phyloseq(OTU, TAX, MAP,Tree)

## Top 10 families ##
family.sum = tapply(taxa_sums(physeq), tax_table(physeq)[, "Family"], sum, na.rm=TRUE)

top10family = names(sort(family.sum, TRUE))[1:10]
physeqfamily10 = prune_taxa((tax_table(physeq)[, "Family"] %in% top10family), physeq)

phyloGlom = tax_glom(physeqfamily10, "Family")
glomTax = tax_table(phyloGlom)[,"Family"]
glomOTU = otu_table(phyloGlom)
glomTable = merge(glomOTU,glomTax,by=0,all=TRUE) 
rownames(glomTable) = glomTable[,"Family"]
glomTable$Row.names = NULL
glomTable$Family = NULL

glomTable2 = glomTable / rep(colSums(glomTable), each = nrow(glomTable))
glomTable3 <- glomTable2 %>% t()

Sample <- rownames(glomTable3)
glomTable3 <- cbind(Sample,glomTable3) %>% as.data.frame()
glomTable3 <- glomTable3[,-1]

indx <- sapply(glomTable3, is.factor)
glomTable3[indx] <- lapply(glomTable3[indx], function(x) as.numeric(as.character(x)))
glomTable3[glomTable3 == 0.000000e+00]<-NA

### Antibiotics ##
Antibioticsraw <- read_excel("Salmonella_OTU_Silva_Dec.xlsx", sheet = "ENV")
rownames(Antibioticsraw) <- Antibioticsraw$SampleID
Antibioticsraw <- Antibioticsraw %>% t() %>% as.data.frame() %>% 
  select(-c(P1:P12),-c(P25:P36),-c(S1:S12),-c(S25:S36)) %>% 
  select(-c(P16:P18),-c(P40:P42),-c(S16:S18),-c(S40:S42))
Antibioticsraw <- Antibioticsraw %>% select(-S21)#didn't amplify
Antibioticsraw <- Antibioticsraw %>% t() %>% as.data.frame()
Antibioticsraw <- Antibioticsraw[,-c(1:4)] %>% select(-Cephalexin)

## Antibiotic resistance genes ##
#Read raw WaferGen Data
#Only positive control for 16S gene is tested
FinalARG19 <- read_excel("SalmonellaARGs.xlsx",sheet = "All_Cts")

#Select data columns that will be used in the pipe
FinalARG19_Tidy <- FinalARG19 %>% 
  dplyr::select(-'Row Labels',-'Note')

FinalARG19_Tidy <- FinalARG19_Tidy[-2,]

#Pull out non numeric Assay column
Assay<-FinalARG19_Tidy %>% 
  dplyr::select(Gene,Description,Type)

#Numeric genes
FinalARG_process<-FinalARG19_Tidy %>% 
  dplyr::select(-Gene,-Description,-Type)

num_row<-nrow(FinalARG_process)
num_column<-ncol(FinalARG_process)

#Cutoff Ct=30. ForCt>30, return NA
for (i in 1:num_row){
  for (j in 1:num_column){
    if (is.na.data.frame(FinalARG_process[i,j])[1]){
      next
    }
    if (FinalARG_process[i,j] > 30){
      FinalARG_process[i,j] <- NA
    }
  }
}

FinalARG_process <- FinalARG_process[,-c(1:6)]
FinalARG_process <- FinalARG_process %>% select(-(S49:S63))#remove the nematode samples
FinalARG_process = FinalARG_process[,mixedsort(colnames(FinalARG_process))]

num_row1<-nrow(FinalARG_process)
num_column1<-ncol(FinalARG_process)

#If 2 in 3 WaferGen triplicates is NA, all of them is NA
for (i in 1:num_row1){
  for (j in 1:(num_column1/3)){
    Process_vector<-FinalARG_process[i, c(j*3-2,j*3-1,j*3)]
    if (sum(is.na.data.frame(Process_vector)) ==2){
      FinalARG_process[i, c(j*3-2,j*3-1,j*3)] <- NA
    }    
  }
}

FinalARG_process1 <- cbind(Assay,FinalARG_process)
#Filter rows with all NAs
FinalARG_process_filtered <- FinalARG_process1[
  rowSums(is.na.data.frame(FinalARG_process1)) < 72, 
]

Assay_filtered<-FinalARG_process_filtered %>% 
  dplyr::select(Gene,Description,Type)

FinalARG_process_filtered <- FinalARG_process_filtered[,-1]
dataframeAVESTDCOPIESRA <- FinalARG_process_filtered %>% select(-Description,-Type)

num_row2<-nrow(dataframeAVESTDCOPIESRA)
num_column2<-ncol(dataframeAVESTDCOPIESRA)

#Calculate Estimated Copies using equation 10^((30-average)/3.3333
estimated_copies<-data.frame(matrix(NA,nrow=19,ncol=72))

for (i in 1:num_row2){
  for (j in 1:num_column2){
    if (is.na.data.frame(dataframeAVESTDCOPIESRA[i,j])){
      next
    }
    estimated_copies[i,j] <- (10^((30-dataframeAVESTDCOPIESRA[i,j])/3.3333))
  }
}

#Calculate Relative Abundance normalized to 16s rRNA gene
relative_abundance<-data.frame(matrix(NA,nrow=19,ncol=72))
for (i in 1:num_row2){
  relative_abundance[i,] = estimated_copies[i,]/estimated_copies[1,]
}

#Add colnames to calculated dataframe
colnames(estimated_copies) <- colnames(dataframeAVESTDCOPIESRA)
colnames(relative_abundance) <- colnames(dataframeAVESTDCOPIESRA)

row.names(Assay_filtered) <- seq(from =1, to =19)

estimatedcopies_downstream <- cbind(Assay_filtered,estimated_copies)
relativeabundance_downstream <- cbind(Assay_filtered,relative_abundance)

relativeabundance_downstream <- relativeabundance_downstream %>% select(-Description,-Type)
rownames(relativeabundance_downstream) <- relativeabundance_downstream$Gene
relativeabundance_downstream <- relativeabundance_downstream[-1,] %>% select(-Gene)
relativeabundance_downstream <- relativeabundance_downstream %>%
  select(-c(P1:P12),-c(P25:P36),-c(S1:S12),-c(S25:S36)) #remove control without antibiotics
relativeabundance_downstream <- relativeabundance_downstream %>% select(-S21)#didn't amplify

#Generated Dataframes
Antibioticsraw
relativeabundance_downstream
glomTable3

#Filter and select Genes appear >=4 in our 7 Phar.samples
relativeabundance_downstream1 <- relativeabundance_downstream[
  rowSums(is.na.data.frame(relativeabundance_downstream)) < 32, 
] %>% 
  t() %>% 
  log2()

Correlation_ARGsRA_Antibiotics_Genes <- cbind(Antibioticsraw,relativeabundance_downstream1,glomTable3) 
correlationscreening_antibiotics <- Correlation_ARGsRA_Antibiotics_Genes
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
outputcorhackARG_Antibiotics<-rcorr(as.matrix(correlationscreening_antibiotics), type="spearman")
OutputtableARG_Antibiotics_p_r<-flattenCorrMatrix(outputcorhackARG_Antibiotics$r,outputcorhackARG_Antibiotics$P)
outputcorhackARG_Antibiotics <- outputcorhackARG_Antibiotics %>% as.matrix()

#Get the strong coorelation and significant output table
Outputtable_significant<-OutputtableARG_Antibiotics_p_r %>% 
  dplyr::rename(pvalue="p", coefficientvalue="cor")

CORRELATION<-
  filter(Outputtable_significant, coefficientvalue < -0.6 | coefficientvalue> 0.6)

CORRELATION1 <- filter(CORRELATION,pvalue<0.05)
#correlation1 contain ARGs, antibiotics, and bacterial community
write.table(CORRELATION1, "CORRELATION1.txt", sep="\t")


#################
rm(list=ls())

#Bacterial Community wihtout antibiotics
OTUtableraw <- read_excel("Salmonella_OTU_Silva_Dec.xlsx", sheet = "OTU_all")
OTUtable <- OTUtableraw %>% filter(!grepl("Archaea",OTUID)) %>% filter(!grepl("Eukaryota",OTUID))
OTU <- OTUtable %>% dplyr::select(-"OTUID")
ROWNAMES <- OTU$FeatureID 
OTU = OTU[,mixedsort(colnames(OTU))]
row.names(OTU) <- OTU$FeatureID 
OTU <- OTU %>% select(-FeatureID) %>% 
  #select(-c(P1:P12),-c(P25:P36),-c(S1:S12),-c(S25:S36)) %>% 
  select(-c(P4:P6),-c(S4:S6),-c(P28:P30),-c(S28:S30)) %>% 
  select(-c(P16:P18),-c(P40:P42),-c(S16:S18),-c(S40:S42))
OTU <- OTU %>% select(-S21)#didn't amplify
OTU[is.na(OTU)] <- 0
OTU <- OTU %>% as.matrix()
row.names(OTU) <- ROWNAMES

map <- read_excel("Salmonella_OTU_Silva_Dec.xlsx", sheet = "ENV")
map <- map[-69,]#didn't amplifyS21
row.names(map) <- map$SampleID
map <- map %>% t() %>% as.data.frame()
map <- map %>% 
  select(-c(P4:P6),-c(S4:S6),-c(P28:P30),-c(S28:S30)) %>%
  select(-c(P16:P18),-c(P40:P42),-c(S16:S18),-c(S40:S42)) %>% t() %>% as.data.frame()
SAMPLEIDROWNAMES <- map$SampleID
map <- map %>% as.matrix() %>% as.data.frame()
row.names(map) <- SAMPLEIDROWNAMES
ROANAMESLISTSAMPLESID <- row.names(map)


Taxonomy <- OTUtable %>% 
  dplyr::select("OTUID")

Taxonomy <- Taxonomy %>%
  separate("OTUID",into=c("Domain","Phylum","Class","Order","Family","Genus","Species"),sep=";")

TAX <- Taxonomy %>% as.data.frame() %>% 
  mutate(Domain=gsub("D_0__","",Domain)) %>% 
  mutate(Phylum=gsub("D_1__","",Phylum)) %>% 
  mutate(Class=gsub("D_2__","",Class)) %>% 
  mutate(Order=gsub("D_3__","",Order)) %>% 
  mutate(Family=gsub("D_4__","",Family)) %>% 
  mutate(Genus=gsub("D_5__","",Genus)) %>% 
  mutate(Species=gsub("D_6__","",Species)) %>% 
  as.matrix() %>% as.data.frame()

TAX[is.na(TAX)] <- NA
TAX <- apply(TAX, 2, function(x) gsub("^$|^ $", NA, x))
TAX <- apply(TAX, 2, function(x) gsub("__", NA, x))
TAX <- apply(TAX, 2, function(x) gsub("uncultured", NA, x))

OTU=otu_table(OTU,taxa_are_rows = TRUE)
rownames(TAX) <- row.names(OTU)
TAX=tax_table(TAX)
MAP <-  sample_data(map)
physeq = merge_phyloseq(OTU, TAX, MAP)
Tree=rtree(ntaxa(physeq),rooted=TRUE,tip.label = taxa_names(physeq))

physeq = phyloseq(OTU, TAX, MAP,Tree)

#Top10 families
family.sum = tapply(taxa_sums(physeq), tax_table(physeq)[, "Family"], sum, na.rm=TRUE)

top10family = names(sort(family.sum, TRUE))[1:10]
physeqfamily10 = prune_taxa((tax_table(physeq)[, "Family"] %in% top10family), physeq)

phyloGlom = tax_glom(physeqfamily10, "Family")
glomTax = tax_table(phyloGlom)[,"Family"]
glomOTU = otu_table(phyloGlom)
glomTable = merge(glomOTU,glomTax,by=0,all=TRUE) 
rownames(glomTable) = glomTable[,"Family"]
glomTable$Row.names = NULL
glomTable$Family = NULL

glomTable2 = glomTable / rep(colSums(glomTable), each = nrow(glomTable))
glomTable3 <- glomTable2 %>% t()

Sample <- rownames(glomTable3)
glomTable3 <- cbind(Sample,glomTable3) %>% as.data.frame()
glomTable3 <- glomTable3[,-1]


#Antibiotic resistance genes
#Read raw WaferGen Data
#Only positive control for 16S gene is tested
FinalARG19 <- read_excel("SalmonellaARGs.xlsx",sheet = "All_Cts")

#Select data columns that will be used in the pipe
FinalARG19_Tidy <- FinalARG19 %>% 
  dplyr::select(-'Row Labels',-'Note')

FinalARG19_Tidy <- FinalARG19_Tidy[-2,]


#Pull out non numeric Assay column
Assay<-FinalARG19_Tidy %>% 
  dplyr::select(Gene,Description,Type)

#Numeric genes
FinalARG_process<-FinalARG19_Tidy %>% 
  dplyr::select(-Gene,-Description,-Type)

num_row<-nrow(FinalARG_process)
num_column<-ncol(FinalARG_process)

#Cutoff Ct=30. ForCt>30, return NA
for (i in 1:num_row){
  for (j in 1:num_column){
    if (is.na.data.frame(FinalARG_process[i,j])[1]){
      next
    }
    if (FinalARG_process[i,j] > 30){
      FinalARG_process[i,j] <- NA
    }
  }
}

FinalARG_process <- FinalARG_process[,-c(1:6)]
FinalARG_process <- FinalARG_process %>% select(-(S49:S63))#remove the nematode samples
FinalARG_process = FinalARG_process[,mixedsort(colnames(FinalARG_process))]

num_row1<-nrow(FinalARG_process)
num_column1<-ncol(FinalARG_process)

#If 2 in 3 WaferGen triplicates is NA, all of them is NA
for (i in 1:num_row1){
  for (j in 1:(num_column1/3)){
    Process_vector<-FinalARG_process[i, c(j*3-2,j*3-1,j*3)]
    if (sum(is.na.data.frame(Process_vector)) ==2){
      FinalARG_process[i, c(j*3-2,j*3-1,j*3)] <- NA
    }    
  }
}
#Assay <- as.matrix(Assay)
#row.names(FinalARG_process) <- Assay[,1]

FinalARG_process1 <- cbind(Assay,FinalARG_process)
#Filter rows with all NAs
FinalARG_process_filtered <- FinalARG_process1[
  rowSums(is.na.data.frame(FinalARG_process1)) < 72, 
]

Assay_filtered<-FinalARG_process_filtered %>% 
  dplyr::select(Gene,Description,Type)

FinalARG_process_filtered <- FinalARG_process_filtered[,-1]
dataframeAVESTDCOPIESRA <- FinalARG_process_filtered %>% select(-Description,-Type)

num_row2<-nrow(dataframeAVESTDCOPIESRA)
num_column2<-ncol(dataframeAVESTDCOPIESRA)

#Calculate Estimated Copies using equation 10^((30-average)/3.3333
estimated_copies<-data.frame(matrix(NA,nrow=19,ncol=72))

for (i in 1:num_row2){
  for (j in 1:num_column2){
    if (is.na.data.frame(dataframeAVESTDCOPIESRA[i,j])){
      next
    }
    estimated_copies[i,j] <- (10^((30-dataframeAVESTDCOPIESRA[i,j])/3.3333))
  }
}

#Calculate Relative Abundance normalized to 16s rRNA gene
relative_abundance<-data.frame(matrix(NA,nrow=19,ncol=72))
for (i in 1:num_row2){
  relative_abundance[i,] = estimated_copies[i,]/estimated_copies[1,]
}

#Add colnames to calculated dataframe
colnames(estimated_copies) <- colnames(dataframeAVESTDCOPIESRA)
colnames(relative_abundance) <- colnames(dataframeAVESTDCOPIESRA)

row.names(Assay_filtered) <- seq(from =1, to =19)
#row.names(Assay_filtered) <- Assay_filtered

estimatedcopies_downstream <- cbind(Assay_filtered,estimated_copies)
relativeabundance_downstream <- cbind(Assay_filtered,relative_abundance)

relativeabundance_downstream <- relativeabundance_downstream %>% select(-Description,-Type)
rownames(relativeabundance_downstream) <- relativeabundance_downstream$Gene
relativeabundance_downstream <- relativeabundance_downstream[-1,] %>% select(-Gene)
#relativeabundance_downstream <- relativeabundance_downstream %>%
#select(-c(P1:P12),-c(P25:P36),-c(S1:S12),-c(S25:S36)) #remove control without antibiotics
relativeabundance_downstream <- relativeabundance_downstream %>% select(-S21)#didn't amplify
#relativeabundance_downstream <- relativeabundance_downstream

#Antibioticsraw
relativeabundance_downstream
glomTable3

#Filter and select Genes appear >=4 in our 7 Phar.samples
#relativeabundance_downstream <- relativeabundance_downstream %>% t()
relativeabundance_downstream1 <- relativeabundance_downstream[
  rowSums(is.na.data.frame(relativeabundance_downstream)) < 69, 
] %>% 
  t() %>% 
  log2()

#glomTable3 <- glomTable3 %>% as.matrix()
#glomTablex <- as.numeric(glomTable3)

indx <- sapply(glomTable3, is.factor)
glomTable3[indx] <- lapply(glomTable3[indx], function(x) as.numeric(as.character(x)))
glomTable3[glomTable3 == 0.000000e+00]<-NA
glomTable3[1:3,] <- NA

Correlation_ARGsRA_Antibiotics_Genes <- cbind(relativeabundance_downstream1,glomTable3) 


correlationscreening_antibiotics <- Correlation_ARGsRA_Antibiotics_Genes
correlationscreening_antibiotics[is.na(correlationscreening_antibiotics)] <- NA

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}


outputcorhackARG_Antibiotics<-rcorr(as.matrix(correlationscreening_antibiotics), type="spearman")
OutputtableARG_Antibiotics_p_r<-flattenCorrMatrix(outputcorhackARG_Antibiotics$r,outputcorhackARG_Antibiotics$P)
outputcorhackARG_Antibiotics <- outputcorhackARG_Antibiotics %>% as.matrix()

#Get the strong coorelation and significant output table
Outputtable_significant<-OutputtableARG_Antibiotics_p_r %>% 
  dplyr::rename(pvalue="p", coefficientvalue="cor")

CORRELATION<-
  filter(Outputtable_significant, coefficientvalue < -0.6 | coefficientvalue> 0.6)

CORRELATION2 <- filter(CORRELATION,pvalue<0.05)
write.table(CORRELATION2, "CORRELATION2.txt", sep="\t")

#Find Node, the ARGs, bacteria, and antibiotics
FINDNODE <- rbind(CORRELATION, CORRELATION2) %>% as.data.frame() %>% dplyr::rename(variable1 = "row") %>% dplyr::rename (variable2 ="column")

a <- unique(FINDNODE$variable1) %>% as.matrix()
b <- unique(FINDNODE$variable2) %>% as.matrix()
c <- c(a,b)
Uniquenode <- unique(c)
write.csv(Uniquenode, file="FINDNODES.csv") 

###Note: This is for full correlation screening for all samples#####
#### to get individual compartment, we need to slice the dataset and run individual compartment related correlations!##

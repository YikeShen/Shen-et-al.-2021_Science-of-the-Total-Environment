#Taxonomy profile
#Figure S1

rm(list=ls())
setwd("/Users/yikeshen/Desktop/Salmonella_Paper")

#########Turn OTU table and Env into Phyloseq object#######
OTUtableraw <- read_excel("Salmonella_OTU_Silva_Dec.xlsx", sheet = "OTU_all")
OTUtable <- OTUtableraw %>% filter(!grepl("Archaea",OTUID)) %>% filter(!grepl("Eukaryota",OTUID))
OTU <- OTUtable %>% dplyr::select(-"OTUID")
ROWNAMES <- OTU$FeatureID 
OTU = OTU[,mixedsort(colnames(OTU))]

dataframeAVESTD <- OTU %>% select(-FeatureID)
dataframeAVESTD[dataframeAVESTD == 0]<-NA
num_col1 <- ncol(dataframeAVESTD)
num_row1 <- nrow(dataframeAVESTD)
Average_df<-data.frame(matrix(0,nrow=3137,ncol=32))
std_df<-data.frame(matrix(0,nrow=3137,ncol=32))

for (i in 1:num_row1){
  for (j in 1:((num_col1)/3)){
    Process_vector<-as.numeric(unlist(dataframeAVESTD[i, c(j*3-2,j*3-1,j*3)]))
    
    if (sum(is.na.data.frame(Process_vector)) ==3){
      Average_df[i,j] = NA
      std_df[i,j] = NA      
    }
    else{
      average <- mean(na.omit(Process_vector))
      std_val <- sd(na.omit(Process_vector))
      Average_df[i,j] = average
      std_df[i,j] = std_val      
    }
  }
}

columnnamessalmonella <- c("D25-Shoots-CK","D27-Shoots-CK","D30-Shoots-CK","D35-Shoots-CK","D25-Shoots-Anti","D27-Shoots-Anti","D30-Shoots-Anti","D35-Shoots-Anti","D25-Roots-CK","D27-Roots-CK","D30-Roots-CK","D35-Roots-CK","D25-Roots-Anti","D27-Roots-Anti","D30-Roots-Anti","D35-Roots-Anti","D25-Soil-CK","D27-Soil-CK","D30-Soil-CK","D35-Soil-CK","D25-Soil-Anti","D27-Soil-Anti","D30-Soil-Anti","D35-Soil-Anti","D25-Rhizo-CK","D27-Rhizo-CK","D30-Rhizo-CK","D35-Rhizo-CK","D25-Rhizo-Anti","D27-Rhizo-Anti","D30-Rhizo-Anti","D35-Rhizo-Anti")

colnames(Average_df) <- columnnamessalmonella
colnames(std_df) <- columnnamessalmonella

#delete lines with all NAs, just screening, there was no NAs
taxainfo <- OTUtable %>% dplyr::select("OTUID","FeatureID")
OTUtablenew <- cbind(taxainfo,Average_df)
OTUtablenew <- OTUtablenew[
  rowSums(is.na.data.frame(OTUtablenew)) < 32, 
]


OTU <- OTUtablenew
row.names(OTU) <- OTU$FeatureID 
OTU <- OTU %>% select(-FeatureID,-OTUID)
OTU[is.na(OTU)] <- 0
OTU <- OTU %>% as.matrix()
row.names(OTU) <- ROWNAMES


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

map <- read_excel("Salmonella_OTU_Silva_Dec.xlsx", sheet = "ENV_avename")
SAMPLEIDROWNAMES <- map$SampleID
map <- map %>% as.matrix() %>% as.data.frame()
row.names(map) <- SAMPLEIDROWNAMES
ROANAMESLISTSAMPLESID <- row.names(map)

OTU=otu_table(OTU,taxa_are_rows = TRUE)
rownames(TAX) <- row.names(OTU)
TAX=tax_table(TAX)
MAP <-  sample_data(map)
physeq = merge_phyloseq(OTU, TAX, MAP)
Tree=rtree(ntaxa(physeq),rooted=TRUE,tip.label = taxa_names(physeq))

physeq = phyloseq(OTU, TAX, MAP,Tree)


########Figure S1, phylum level microbial community structure########
phylum.sum = tapply(taxa_sums(physeq), tax_table(physeq)[, "Phylum"], sum, na.rm=TRUE)
top10phyla = names(sort(phylum.sum, TRUE))[1:10]
physeqphylum10 = prune_taxa((tax_table(physeq)[, "Phylum"] %in% top10phyla), physeq)

phyloGlom = tax_glom(physeqphylum10, "Phylum")
glomTax = tax_table(phyloGlom)[,"Phylum"]
glomOTU = otu_table(phyloGlom)
glomTable = merge(glomOTU,glomTax,by=0,all=TRUE) 
rownames(glomTable) = glomTable[,"Phylum"]
glomTable$Row.names = NULL
glomTable$Phylum = NULL

glomTable2 = glomTable / rep(colSums(glomTable), each = nrow(glomTable))
glomTable3 <- glomTable2 %>% t()

Com <- map$Compartment
Sample <- rownames(glomTable3)
glomTable3 <- cbind(Com,Sample,glomTable3) %>% as.data.frame()

MicrobialCommunity2 <-glomTable3  %>% 
  gather(3:12,key="Bacteria_Phylum",value ="Percentage")

MicrobialCommunity2$Percentage <- as.numeric(as.character(MicrobialCommunity2$Percentage))
MicrobialCommunity2$Com <- factor(x = MicrobialCommunity2$Com,levels = c("Soil","Rhizosphere","Roots","Shoots"))
MicrobialCommunity2$Sample <- factor(x = MicrobialCommunity2$Sample,levels = c("D25-Shoots-CK","D27-Shoots-CK","D30-Shoots-CK","D35-Shoots-CK","D25-Shoots-Anti","D27-Shoots-Anti","D30-Shoots-Anti","D35-Shoots-Anti","D25-Roots-CK","D27-Roots-CK","D30-Roots-CK","D35-Roots-CK","D25-Roots-Anti","D27-Roots-Anti","D30-Roots-Anti","D35-Roots-Anti","D25-Soil-CK","D27-Soil-CK","D30-Soil-CK","D35-Soil-CK","D25-Soil-Anti","D27-Soil-Anti","D30-Soil-Anti","D35-Soil-Anti","D25-Rhizo-CK","D27-Rhizo-CK","D30-Rhizo-CK","D35-Rhizo-CK","D25-Rhizo-Anti","D27-Rhizo-Anti","D30-Rhizo-Anti","D35-Rhizo-Anti"))

Structurefigure2<- ggplot()+ theme_bw()+
  geom_bar(aes(x=Sample, y= Percentage*100, fill=Bacteria_Phylum), data = MicrobialCommunity2,stat="identity") +
  labs(x="Sample Name", y="Percentage (%)") +
  scale_y_continuous(labels = dollar_format(suffix = "", prefix = "")) +
  ggtitle("Microbial Community Structure (Phylum)") +
  theme(legend.position="top",axis.text=element_text(size=16),
        axis.title=element_text(size=16),plot.title = element_text(size=22, hjust=0.5),text = element_text(size = 16),
        axis.text.x = element_text(angle = 90))+
  facet_wrap(~Com,1,scales = "free_x")+scale_fill_d3()
Structurefigure2

########Figure S1, family level microbial community structure########
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

Com <- map$Compartment
Sample <- rownames(glomTable3)
glomTable3 <- cbind(Com,Sample,glomTable3) %>% as.data.frame()
View(glomTable3)

MicrobialCommunity <-glomTable3  %>% 
  gather(3:12,key="Bacteria_Family",value ="Percentage")

MicrobialCommunity$Percentage <- as.numeric(as.character(MicrobialCommunity$Percentage))
MicrobialCommunity$Com <- factor(x = MicrobialCommunity2$Com,levels = c("Soil","Rhizosphere","Roots","Shoots"))
MicrobialCommunity$Sample <- factor(x = MicrobialCommunity2$Sample,levels = c("D25-Shoots-CK","D27-Shoots-CK","D30-Shoots-CK","D35-Shoots-CK","D25-Shoots-Anti","D27-Shoots-Anti","D30-Shoots-Anti","D35-Shoots-Anti","D25-Roots-CK","D27-Roots-CK","D30-Roots-CK","D35-Roots-CK","D25-Roots-Anti","D27-Roots-Anti","D30-Roots-Anti","D35-Roots-Anti","D25-Soil-CK","D27-Soil-CK","D30-Soil-CK","D35-Soil-CK","D25-Soil-Anti","D27-Soil-Anti","D30-Soil-Anti","D35-Soil-Anti","D25-Rhizo-CK","D27-Rhizo-CK","D30-Rhizo-CK","D35-Rhizo-CK","D25-Rhizo-Anti","D27-Rhizo-Anti","D30-Rhizo-Anti","D35-Rhizo-Anti"))

Structurefigure<- ggplot()+ theme_bw()+
  geom_bar(aes(x=Sample, y= Percentage*100, fill=Bacteria_Family), data = MicrobialCommunity,stat="identity") +
  labs(x="Sample Name", y="Percentage (%)") +
  scale_y_continuous(labels = dollar_format(suffix = "", prefix = "")) +
  ggtitle("Microbial Community Structure (Family)") +
  theme(legend.position="top",axis.text=element_text(size=16),
        axis.title=element_text(size=16),plot.title = element_text(size=22, hjust=0.5),text = element_text(size = 16),
        axis.text.x = element_text(angle = 90))+
  facet_wrap(~Com,1,scales = "free_x")+scale_fill_d3()
Structurefigure

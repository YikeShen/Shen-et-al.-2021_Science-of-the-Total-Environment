rm(list=ls())
#Raw LEFSE was calculated from online galaxy platform developed at Huttenhower lab: https://huttenhower.sph.harvard.edu/lefse/ and cleaned up using unix command.

library(ggplot2)
library(scales)

LEFSE_raw <- read_excel("Salmonella_OTU_Silva_Dec.xlsx", sheet = "LEFSE")
Taxonomy <- LEFSE_raw %>% as.data.frame() %>% 
  separate(Taxon,into=c("Domain","Phylum","Class","Order","Family","Genus","Species"),sep="([\\.])")

Taxonomy <-  Taxonomy %>% dplyr::select(-Genus,-Species,-Domain)

#shoot
Shoot <- Taxonomy[73:75,] %>% dplyr::select(-log10,-p_kruskalwallis)
Shoot <- Shoot[order(Shoot$Phylum),]
Shoot <- Shoot[order(Shoot$Treatment),]
#soil
Soil <- Taxonomy[1:72,]
Soil <- Soil[order(Soil$Phylum),]
Soil <- Soil%>% dplyr::select(-log10,-p_kruskalwallis) %>% 
  filter(!Class=="Subgroup6") %>% 
  filter(LDA>3.5) %>% 
  distinct() %>% 
  na.omit(Family)
Soil<-Soil[!duplicated(Soil$Family),]  
Soil <- Soil[order(Soil$Treatment),]

#rhizo
Rhizo <- Taxonomy[76:113,]
Rhizo <- Rhizo[order(Rhizo$Phylum),]
Rhizo <- Rhizo%>% dplyr::select(-log10,-p_kruskalwallis) %>% 
  filter(!Class=="Subgroup6") %>% 
  filter(LDA>3.5) %>% 
  distinct() %>% 
  na.omit(Family)
Rhizo<-Rhizo[!duplicated(Rhizo$Family),]  
Rhizo <- Rhizo[order(Rhizo$Treatment),]

#root
roots <- Taxonomy[114:141,]
roots <- roots[order(roots$Phylum),]
roots <- roots%>% dplyr::select(-log10,-p_kruskalwallis) %>% 
  filter(!Class=="Subgroup6") %>% 
  filter(LDA>3.5) %>% 
  distinct() %>% 
  na.omit(Family)
roots<-roots[!duplicated(roots$Family),]  
roots <- roots[order(roots$Treatment),]

Taxonomy <- rbind(Soil,Rhizo,roots,Shoot)
row.names(Taxonomy) <- NULL
Taxonomy$Phylum <- paste(Taxonomy$Phylum, ".", Taxonomy$Class,".",Taxonomy$Order,".",Taxonomy$Family)
theme_set(theme_bw())
Taxonomy$Compartment <- factor(x = Taxonomy$Compartment,levels = c("Soil","Rhizosphere","Roots","Shoots"))

LDAPLOT <- Taxonomy %>% ggplot()+geom_bar(aes(x=Phylum,y=LDA,fill=Treatment),color="black",stat="identity")+scale_fill_manual(breaks = c("Antibiotics", "Control"),values=c("red", "blue"))+facet_wrap(~Compartment,ncol=1,scales="free_y")+labs(x="Taxon", y=expression("LDA Score (log10)"))+ggtitle("LDA analysis")+theme(legend.position="top",axis.text=element_text(size=12),axis.title=element_text(size=16),plot.title = element_text(size=22, hjust=0.5),text = element_text(size = 16),axis.text.x = element_text(angle = 0))+coord_flip()
LDAPLOT
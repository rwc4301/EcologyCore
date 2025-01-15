#ACHTUNG: From the RStudio menu, click on "Session" and then "Set Working Directory" to "To Source File Location"
#Script for visualisation of null modelling (beta NTI, Raup-Crick Beta-Diversity, and Elements of Metacommunity)
#Reference: http://uu.diva-portal.org/smash/get/diva2:1373632/DATASET09.txt
#Authors: Umer, Anna, and Simon
#Versions: 1.2 (fixed ordering issue)

library(ggplot2)
library(yarrr)

#PARAMETERS ###########################
label="Hypothesis1"
PairwisebNTI<-read.csv("PairwisebNTI_Hypothesis1.csv",header=T,row.names=1)
PairwiseRC<-read.csv("PairwiseRC_Hypothesis1.csv",header=T,row.names=1)
RC<-read.csv("RC_Hypothesis1.csv",header=T,row.names=1)
Coherence<-read.csv("Coherence_Hypothesis1.csv",header=T,row.names=1)
BoundaryClump<-read.csv("Boundary_Hypothesis1.csv",header=T,row.names=1)
Turnover<-read.csv("Turnover_Hypothesis1.csv",header=T,row.names=1)
ordering<-c(
  "COV",
  "CH"
)
colours <- c(
  "red",
  "blue",
  #Next colors are for lines mainly used in the PCoA script
  "#000080","#4876FF","#CAE1FF","#9FB6CD","#1E90FF","#00F5FF","#00C957",grey.colors(1000));
RC_width=8
RC_height=8
QPE_width=8
QPE_height=8
EMS_width=8
EMS_height=8
#PARAMETERS ###########################

QPE_table<-cbind(Var1=as.character(PairwiseRC$Var1), Var2=as.character(PairwiseRC$Var2),bNTI=PairwisebNTI[,"value",drop=F], RC=PairwiseRC[,"value",drop=F], Groups=PairwiseRC[,"Groups",drop=F])
#Now we want to get rid of any rows that have NA there to get the comparisons from
# N x N to N (N-1)/2. A way to do this to use complete.cases
QPE_table<-QPE_table[complete.cases(QPE_table),]
names(QPE_table)<-c("Var1","Var2","bNTI","RC","Groups")
QPE_table$Groups<-factor(as.character(QPE_table$Groups))

QPE_df<-NULL

for(i in levels(QPE_table$Groups)){
  tmp<-QPE_table[QPE_table$Groups==i,]
  sp_mask<-abs(tmp$bNTI)>2
  hs_count<-sum(tmp$bNTI[sp_mask]<2)
  vs_count<-sum(tmp$bNTI[sp_mask]>2)
  sig_count<-sum(sp_mask)
  total_count<-nrow(tmp)
  nonsig_count<-nrow(tmp[!sp_mask,])
  dl_count<-sum(tmp[!sp_mask,"RC"]>0.95)
  hd_count<-sum(tmp[!sp_mask,"RC"]<(-0.95))
  ed_count<-nonsig_count-dl_count-hd_count
  tmp2<-data.frame(measure="Homogeneous Selection",value=(hs_count/total_count*100),Groups=i)
  tmp2<-rbind(tmp2,data.frame(measure="Variable Selection",value=(vs_count/total_count*100),Groups=i))
  tmp2<-rbind(tmp2,data.frame(measure="Dispersal Limitation",value=(dl_count/total_count*100),Groups=i))
  tmp2<-rbind(tmp2,data.frame(measure="Undominated",value=(ed_count/total_count*100),Groups=i))
  tmp2<-rbind(tmp2,data.frame(measure="Homogenizing Dispersal",value=(hd_count/total_count*100),Groups=i))
  if(is.null(QPE_df)){QPE_df<-tmp2} else {QPE_df<-rbind(QPE_df,tmp2)}
}

#Change the orderign of the Groups
QPE_df$Groups<-factor(as.character(QPE_df$Groups),levels=ordering)
pdf(paste("QPE_",label,".pdf",sep=""),width=QPE_width,height=QPE_height)
p<-ggplot(QPE_df, aes(x=Groups, y=value, fill=Groups)) 
p<-p+geom_bar(stat="identity")+theme_minimal()
p<-p+geom_text(aes(label=sprintf("%0.2f", round(value, digits = 2))), vjust=-0.3, size=3.5)
p<-p+facet_wrap(~measure, strip.position="left", ncol=1,scales="free_y")
p<-p+scale_fill_manual(values=colours)
p<-p+ylim(0,110)
p<-p+ylab("% Assembly Processes")
p<-p+theme(strip.background = element_rect(fill="white"))+theme(panel.spacing = unit(2, "lines"),
                                                                axis.text.x = element_text(angle = 90, hjust = 1))
print(p)
dev.off()

#We move on to EMS (Elements of Metacommunity Structure)
#UMER: Not confident about interpretation, MAY BE BUGGY!
collated_community_types<-NULL
for(i in rownames(Coherence)){
  community_type="Random"
  if(Coherence[i,"p"]<0.05){
    #Top left figure
    if(Coherence[i,"z"]<(-1.96)){
      community_type="Checkerboard"
    } else if(Coherence[i,"z"]>1.96) {
      #Top middle left  
      if(Turnover[i,"z"]<(-1.96)){
        if(BoundaryClump[i,"p"]<0.05){
          if(BoundaryClump[i,"index"]<0){
            community_type="Nested Hyperdispersed species loss"
          } else {
            community_type="Nested Clumped species loss"
          }
        }
        else{
          community_type="Nested Random species loss"
        }
      }
      #Top middle right 
      else if(Turnover[i,"z"]>1.96)
      {
        if(BoundaryClump[i,"p"]<0.05){
          if(BoundaryClump[i,"index"]<0){
            community_type="Evenly spaced"
          } else {
            community_type="Clementsian"
          }
        }
        else{
          community_type="Gleasonian"
        }
      }
      if(Turnover[i,"p"]>0.05){
        community_type=paste("Quasi-structure",community_type)
      }
    }
    
  }
  #Collate all the information together
  if(is.null(collated_community_types)){collated_community_types<-community_type} else {collated_community_types<-c(collated_community_types,community_type)}
}

EMS<-cbind(Coherence,Turnover,BoundaryClump,Metacommunity=collated_community_types)
names(EMS)<-c("Coherence_Abs","Coherence_z","Coherence_p","Coherence_simMean","Coherence_simVariance",
              "Turnover_turn","Turnover_z","Turnover_p","Turnover_simMean","Turnover_simVariance",
              "Clumping_index","Clumping_p","Clumping_df","Metacommunity")
write.csv(EMS,file=paste("EMS_",label,".csv",sep=""))

EMS<-cbind(EMS,Groups=rownames(EMS))
#Change the ordering of the Groups
EMS$Groups<-factor(as.character(EMS$Groups),levels=ordering)
pdf(paste("EMS_",label,"_DONOTUSE.pdf",sep=""),width=EMS_width,height=EMS_height)
p<-ggplot(EMS,aes(Groups,Coherence_z,color=Turnover_z,size=Clumping_index,shape=Metacommunity))
p<-p+geom_point()
p <- p + geom_hline(yintercept = 1.96,linetype="dotted")
p <-p +geom_text(aes(x=1,y=1.96, label="1.96\n"), colour="blue",size=2) 

p <- p + geom_hline(yintercept = -1.96,linetype="dotted")
p <-p +geom_text(aes(x=1,y=-1.96, label="\n-1.96"), colour="blue",size=2) 

p <- p + ylab("Coherence (z-value)")
p <- p + scale_color_continuous("Turnover (z-value)")
p <- p + scale_size_continuous("Boundary clumping (Morisita's index)")
p<-p+theme(strip.background = element_rect(fill="white"))+theme(panel.spacing = unit(2, "lines"),
                                                                axis.text.x = element_text(angle = 90, hjust = 1))

p<-p+theme_minimal()
print(p)
dev.off()

RC<-cbind(RC,Groups=rownames(RC))
#Change the ordering of the Groups
RC$Groups<-factor(as.character(RC$Groups),levels=ordering)
pdf(paste("RC_",label,".pdf",sep=""),width=RC_width,height=RC_height)
p<-ggplot(RC,aes(Groups,value,colour=Groups))
p<-p+geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1, lty=1)+
  geom_point(size=5)

p <- p + geom_hline(yintercept = 0,linetype="dotted")
p <-p +geom_text(aes(x=1,y=0, label="0\n"), colour="blue",size=2) 

p <- p + geom_hline(yintercept = 1,linetype="dotted")
p <-p +geom_text(aes(x=1,y=1, label="+1\n"), colour="blue",size=2) 

p <- p + geom_hline(yintercept = -1,linetype="dotted")
p <-p +geom_text(aes(x=1,y=-1, label="\n-1"), colour="blue",size=2) 

p<-p+theme_minimal()
p<-p+scale_colour_manual(values=colours)
p<-p+ylab("Incidence-based beta-diversity (±SE)")
p<-p+theme(strip.background = element_rect(fill="white"))+theme(panel.spacing = unit(2, "lines"),
                                                                axis.text.x = element_text(angle = 90, hjust = 1))

print(p)
dev.off()


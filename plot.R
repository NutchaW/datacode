result1<-read.csv("./results_1000/result1.csv")
result2<-read.csv("./results_1000/result2.csv")
result3<-read.csv("./results_1000/result3.csv")
result4<-read.csv("./results_1000/result4.csv")
result5<-read.csv("./results_1000/result5.csv")
result6<-read.csv("./results_1000/result6.csv")
result7<-read.csv("./results_1000/result7.csv")
result8<-read.csv("./results_1000/result8.csv")
result9<-read.csv("./results_1000/result9.csv")
result10<-read.csv("./results_1000/result10.csv")
result11<-read.csv("./results_1000/result11.csv")
result12<-read.csv("./results_1000/result12.csv")
result13<-read.csv("./results_1000/result13.csv")
result14<-read.csv("./results_1000/result14.csv")
result15<-read.csv("./results_1000/result15.csv")
result16<-read.csv("./results_1000/result16.csv")

# transform data
humid<-read.csv("./AbsHumidity.csv")
states<-names(humid)[2:17]

plotDat <-function(data,state,statenum){
library(tidyverse)
check<-data %>%
  select(-"V1001")%>%
  gather(ensemble,value,V1:V1000,factor_key=TRUE) %>%
  filter(!(value>30))
check2<-data %>%
  select("X","V1001")%>%
  gather(ensemble,value,V1001,factor_key=TRUE)
ggplot()+
  geom_line(data=check,aes(x=X,y=value,color=ensemble),alpha=0.3,size=0.08,show.legend=FALSE)+
  geom_line(data=check2,aes(x=X,y=value),color="black",show.legend=FALSE)+
  geom_point(data=check2,aes(x=X,y=value),color="black",size=1,show.legend=FALSE)+
  labs(x="Step Ahead",y="Peak wILI", title=state[statenum])+
  theme(plot.title = element_text(size=10))
}
 
#### plot
loc1<-plotDat(result1,states,1)
loc2<-plotDat(result1,states,2)
loc3<-plotDat(result1,states,3)
loc4<-plotDat(result1,states,4)
loc5<-plotDat(result1,states,5)
loc6<-plotDat(result1,states,6)
loc7<-plotDat(result1,states,7)
loc8<-plotDat(result1,states,8)
loc9<-plotDat(result1,states,9)
loc10<-plotDat(result1,states,10)
loc11<-plotDat(result1,states,11)
loc12<-plotDat(result1,states,12)
loc13<-plotDat(result1,states,13)
loc14<-plotDat(result1,states,14)
loc15<-plotDat(result1,states,15)
loc16<-plotDat(result1,states,16)

library(gridExtra)
grid.arrange(loc1,loc2,loc3,loc4,nrow=2)
grid.arrange(loc5,loc6,loc7,loc8,nrow=2)
grid.arrange(loc9,loc10,loc11,loc12,nrow=2)
grid.arrange(loc13,loc14,loc15,loc16,nrow=2)


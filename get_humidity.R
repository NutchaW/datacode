library(surveillance)
data("campyDE")

abs_humid <- campyDE$hum
abs_humid <- (abs_humid[157:521]/1000)
hist(abs_humid)
AbsHumidity <- data.frame(matrix(rep(abs_humid,16),ncol=16))
location<-read.csv("/Users/ahctun_woon/Desktop/SPRING_19/Sem/german-flu-forecasting/data/GER_states_adjacency.csv")
names(AbsHumidity)<-names(location)[2:17]
write.csv(AbsHumidity, file = "AbsHumidity.csv")

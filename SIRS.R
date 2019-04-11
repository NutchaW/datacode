library(matrixStats)
library(abind)

SIRS <- function(x,ts,dt,tmstep,N,AH,seed,discrete) {
  # x=[S,I,total incidence,R0max,R0min,L,D]
  if (discrete==0) {
    AH=c(AH,AH)
    num_ens=dim(x)[2]
    BT1=matrix(0,length(AH),num_ens)
    BETA=matrix(0,length(AH),num_ens)
    for (i in 1:num_ens) {
      # b=log(x[4,i]-x[5,i])
      if ((x[4,i]-x[5,i])>=0) {
        b=log(x[4,i]-x[5,i])
      }
      else {
        b=complex(real=log(-(x[4,i]-x[5,i])),imaginary = pi)
      }
      a=-180
      BT1[,i]=exp(a*AH+b)+x[5,i]
      BETA[,i]=BT1[,i]/x[7,i]
    }
    # Sr=array(0,c(3,num_ens,tmstep+1))
    # Incidence=matrix(0,(tmstep+1),num_ens)
    Sr=array(0,c(3,num_ens,abs(tmstep)+1))
    Incidence=matrix(0,(abs(tmstep)+1),num_ens)
    Sr[,,1]=x[1:3,]
    L=x[6,]
    D=x[7,]
    # start integration
    tcnt=0
    for (t in seq(ts+dt,ts+tmstep,by = dt)) {
      tcnt=tcnt+1
      Eimmloss=dt*((N*rep(1,num_ens)-Sr[1,,tcnt]-Sr[2,,tcnt])/L)
      tempind=Mod(dt*(BETA[t,]*Sr[1,,tcnt]*Sr[2,,tcnt]/N))>Mod(Sr[1,,tcnt])
      Einf=dt*(BETA[t,]*Sr[1,,tcnt]*Sr[2,,tcnt]/N)
      Einf[tempind]=Sr[1,,tcnt][tempind]
      tempind2=Mod(dt*(Sr[2,,tcnt]/D))>Mod(Sr[2,,tcnt])
      Erecov=dt*(Sr[2,,tcnt]/D)
      Erecov[tempind2]=Sr[2,,tcnt][tempind2]
      Eimmloss[Mod(Eimmloss)<0]=0
      Einf[Mod(Einf)<0]=0
      Erecov[Mod(Erecov)<0]=0
      sk1=Eimmloss-Einf
      ik1=Einf-Erecov
      ik1i=Einf
      Ts1=Sr[1,,tcnt]+round(sk1/2)
      Ti1=Sr[2,,tcnt]+round(ik1/2)
      Eimmloss=dt*((N*rep(1,num_ens)-Ts1-Ti1)/L)
      tempind=Mod(dt*(BETA[t,]*Ts1*Ti1/N))>Mod(Ts1)
      Einf=dt*(BETA[t,]*Ts1*Ti1/N)
      Einf[tempind]=Ts1[tempind]
      tempind2=Mod(dt*(Ti1/D))>Mod(Ti1)
      Erecov=dt*(Ti1/D)
      Erecov[tempind2]=Ti1[tempind2]
      Eimmloss[Mod(Eimmloss)<0]=0
      Einf[Mod(Einf)<0]=0
      Erecov[Mod(Erecov)<0]=0
      sk2=Eimmloss-Einf
      ik2=Einf-Erecov
      ik2i=Einf
      Ts2=Sr[1,,tcnt]+round(sk2/2)
      Ti2=Sr[2,,tcnt]+round(ik2/2)
      Eimmloss=dt*((N*rep(1,num_ens)-Ts2-Ti2)/L)
      tempind=Mod(dt*(BETA[t,]*Ts2*Ti2/N))>Mod(Ts2)
      Einf=dt*(BETA[t,]*Ts2*Ti2/N)
      Einf[tempind]=Ts2[tempind]
      tempind2=Mod(dt*(Ti2/D))>Mod(Ti2)
      Erecov=dt*(Ti2/D)
      Erecov[tempind2]=Ti2[tempind2]
      Eimmloss[Mod(Eimmloss)<0]=0
      Einf[Mod(Einf)<0]=0
      Erecov[Mod(Erecov)<0]=0
      sk3=Eimmloss-Einf
      ik3=Einf-Erecov
      ik3i=Einf
      Ts3=Sr[1,,tcnt]+round(sk3)
      Ti3=Sr[2,,tcnt]+round(ik3)
      Eimmloss=dt*((N*rep(1,num_ens)-Ts3-Ti3)/L)
      tempind=Mod(dt*(BETA[t,]*Ts3*Ti3/N))>Mod(Ts3)
      Einf=dt*(BETA[t,]*Ts3*Ti3/N)
      Einf[tempind]=Ts3[tempind]
      tempind2=Mod(dt*(Ti3/D))>Mod(Ti3)
      Erecov=dt*(Ti3/D)
      Erecov[tempind2]=Ti3[tempind2]
      Eimmloss[Mod(Eimmloss)<0]=0
      Einf[Mod(Einf)<0]=0
      Erecov[Mod(Erecov)<0]=0
      sk4=Eimmloss-Einf
      ik4=Einf-Erecov
      ik4i=Einf
      Sr[1,,tcnt+1]=Sr[1,,tcnt]+round(sk1/6+sk2/3+sk3/3+sk4/6)-seed
      Sr[2,,tcnt+1]=Sr[2,,tcnt]+round(ik1/6+ik2/3+ik3/3+ik4/6)+seed
      Incidence[tcnt+1,]=round(ik1i/6+ik2i/3+ik3i/3+ik4i/6)+seed
    }
    x[1,]=Sr[1,,tcnt+1]
    x[2,]=Sr[2,,tcnt+1]
    x[3,]=colSums(Incidence)/(N/100000)
  } 
  
  if (discrete==1) {
    AH=c(AH,AH)
    num_ens=dim(x)[2]
    BT1=matrix(0, length(AH),num_ens)
    BETA=matrix(0,length(AH),num_ens)
    for (i in 1:num_ens) {
      # b=log(x[4,i]-x[5,i])
      if ((x[4,i]-x[5,i])>=0) {
        b=log(x[4,i]-x[5,i])
      }
      else {
        b=complex(real=log(-(x[4,i]-x[5,i])),imaginary = pi)
      }
      a=-180
      BT1[,i]=exp(a*AH+b)+x[5,i]
      BETA[,i]=BT1[,i]/x[7,i]
    }
    Sr=array(0,c(3,num_ens,abs(tmstep)+1))
    Incidence=matrix(0,abs(tmstep)+1,num_ens)
    Sr[,,1]=x[1:3,]
    L=x[6,]
    D=x[7,]
    # start integration
    tcnt=0
    for (t in seq(ts+dt,ts+tmstep,by = dt)) {
      tcnt=tcnt+1
      Eimmloss=dt*((N*rep(1,num_ens)-Sr[1,,tcnt]-Sr[2,,tcnt])/L)
      tempind=Mod(dt*(BETA[t,]*Sr[1,,tcnt]*Sr[2,,tcnt]/N))>Mod(Sr[1,,tcnt])
      Einf=dt*(BETA[t,]*Sr[1,,tcnt]*Sr[2,,tcnt]/N)
      Einf[tempind]=Sr[1,,tcnt][tempind]
      tempind2=Mod(dt*(Sr[2,,tcnt]/D))>Mod(Sr[2,,tcnt])
      Erecov=dt*(Sr[2,,tcnt]/D)
      Erecov[tempind2]=Sr[2,,tcnt][tempind2]
      Eimmloss[Mod(Eimmloss)<0]=0
      Einf[Mod(Einf)<0]=0
      Erecov[Mod(Erecov)<0]=0
      l1=rbind(Eimmloss,Einf,Erecov)
      l <- rpois(length(l1), l1)
      dim(l) <- dim(l1)
      sk1=l[1,]-l[2,]
      ik1=l[2,]-l[3,]
      ik1i=l[2,]
      Ts1=Sr[1,,tcnt]+round(sk1/2)
      Ti1=Sr[2,,tcnt]+round(ik1/2)
      Eimmloss=dt*((N*rep(1,num_ens)-Ts1-Ti1)/L)
      tempind=Mod(dt*(BETA[t,]*Ts1*Ti1/N))>Mod(Ts1)
      Einf=dt*(BETA[t,]*Ts1*Ti1/N)
      Einf[tempind]=Ts1[tempind]
      tempind2=Mod(dt*(Ti1/D))>Mod(Ti1)
      Erecov=dt*(Ti1/D)
      Erecov[tempind2]=Ti1[tempind2]
      Eimmloss[Mod(Eimmloss)<0]=0
      Einf[Mod(Einf)<0]=0
      Erecov[Mod(Erecov)<0]=0
      l1=rbind(Eimmloss,Einf,Erecov)
      l <- rpois(length(l1), l1)
      dim(l) <- dim(l1)
      sk2=l[1,]-l[2,]
      ik2=l[2,]-l[3,]
      ik2i=l[2,]
      Ts2=Sr[1,,tcnt]+round(sk2/2)
      Ti2=Sr[2,,tcnt]+round(ik2/2)
      Eimmloss=dt*((N*rep(1,num_ens)-Ts2-Ti2)/L)
      tempind=Mod(dt*(BETA[t,]*Ts2*Ti2/N))>Mod(Ts2)
      Einf=dt*(BETA[t,]*Ts2*Ti2/N)
      Einf[tempind]=Ts2[tempind]
      tempind2=Mod(dt*(Ti2/D))>Mod(Ti2)
      Erecov=dt*(Ti2/D)
      Erecov[tempind2]=Ti2[tempind2]
      Eimmloss[Mod(Eimmloss)<0]=0
      Einf[Mod(Einf)<0]=0
      Erecov[Mod(Erecov)<0]=0
      l1=rbind(Eimmloss,Einf,Erecov)
      l <- rpois(length(l1), l1)
      dim(l) <- dim(l1)
      sk3=l[1,]-l[2,]
      ik3=l[2,]-l[3,]
      ik3i=l[2,]
      Ts3=Sr[1,,tcnt]+round(sk3)
      Ti3=Sr[2,,tcnt]+round(ik3)
      Eimmloss=dt*((N*rep(1,num_ens)-Ts3-Ti3)/L)
      tempind=Mod(dt*(BETA[t,]*Ts3*Ti3/N))>Mod(Ts3)
      Einf=dt*(BETA[t,]*Ts3*Ti3/N)
      Einf[tempind]=Ts3[tempind]
      tempind2=Mod(dt*(Ti3/D))>Mod(Ti3)
      Erecov=dt*(Ti3/D)
      Erecov[tempind2]=Ti3[tempind2]
      Eimmloss[Mod(Eimmloss)<0]=0
      Einf[Mod(Einf)<0]=0
      Erecov[Mod(Erecov)<0]=0
      l1=rbind(Eimmloss,Einf,Erecov)
      l <- rpois(length(l1), l1)
      dim(l) <- dim(l1)
      sk4=l[1,]-l[2,]
      ik4=l[2,]-l[3,]
      ik4i=l[2,]
      travel=rpois(1,seed)
      Sr[1,,tcnt+1]=Sr[1,,tcnt]+round(sk1/6+sk2/3+sk3/3+sk4/6)-travel
      Sr[2,,tcnt+1]=Sr[2,,tcnt]+round(ik1/6+ik2/3+ik3/3+ik4/6)+travel
      Incidence[tcnt+1,]=round(ik1i/6+ik2i/3+ik3i/3+ik4i/6)+travel
    }

    x[1,]=Sr[1,,tcnt+1]
    x[2,]=Sr[2,,tcnt+1]
    x[3,]=colSums(Incidence)/(N/100000)

  }
  
  return (x)
}

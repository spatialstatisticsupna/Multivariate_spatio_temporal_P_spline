################################################################################
## Title: A novel multivariate spatio-temporal approach based on splines to   ##
##        analyse different crimes against women                              ##
##                                                                            ##
## Authors: Vicente, G. - Goicoa, T.- Ugarte, M.D.                            ##
##                                                                            ##
## doi:                                                                       ##
##                                                                            ##
################################################################################
rm(list=ls())

## libraries
library(sp); library(tmap); library(INLA); library(grid); library(RColorBrewer)
library(ggplot2); library(ggpubr)

### save figures
if(!file.exists("figures")) {dir.create("figures")}
################################################################################
## data loading and organization                                              ##
################################################################################
## Load data and Maharashtra SpatialPolygonsDataFrame
load("data_Psplines.RData")

## ID
ID<- data.frame(dist=carto$dist,ID_area=carto$ID_area)

## t.from and t.to 
t.from<- min(data$year)
t.to<- max(data$year)
x<- t.from:t.to

## colors 
selected_colors<-c(rgb(154,192,205,alpha=150, maxColorValue=255),
                   rgb(69,139,116, alpha=150, maxColorValue=255),
                   rgb(255,160,122,alpha=150, maxColorValue=255),
                   rgb(147,112,219,alpha=150, maxColorValue=255) )


## Number of areas and number of time periods
n<- length(unique(data$ID_area))
t<- length(unique(data$ID_year))

## selected crimes
crimes<- c("rape", "assault", "cruelty", "kidnapping")
e_crimes<- paste0("e_",crimes)
crime_names<-c("Rape", "Assault", "Cruelty","Kidnapping")

## Parameter of the model (number of crimes)
nc <- length(crimes)

## data.frame ID_disease
data<- data[order(data$ID_year, data$ID_area),]
d.data_frame<- NULL
for(i in 1:t){
  for(j in 1:length(crimes)){
    d1<- data[data$ID_year==i, c(crimes[j], e_crimes[j], "ID_area", "ID_year" )]
    d1$ID_disease<- j
    colnames(d1)<- c("OBS", "EXP", "ID_area", "ID_year", "ID_disease")
    d.data_frame<- rbind(d.data_frame, d1)
    rm(d1)
  }
}
rm(list = c("i","j"))

## orden: disease - year - area
d.data_frame<- d.data_frame[order(d.data_frame$ID_disease,d.data_frame$ID_year,d.data_frame$ID_area),]
################################################################################
## Load results                                                               ##
################################################################################
## Selected strategy
strategy <- "simplified.laplace" # see 'control.inla' (?control.inla)

## selected type of interaction
type<- 2 # 1: Type I; 2: Type II; 3: Type III; 4: Type IV 
########################
## Centered Psplines
########################
for(ms in 1:2){
  for(mt in 1:2){
    resulta<-list()
    ## with temporal correlation
    load(paste0("resul/resul_Centered_Psplines_ps",ms,"_pt",mt,"_type",type,"_",gsub("\\.", "_", strategy),".RData"))
    resulta[[1]]<- model
    rm(model)
    ## without temporal correlation
    load(paste0("resul/resul_Centered_Psplines_ps",ms,"_pt",mt,"_type",type,"_Indep_",gsub("\\.", "_", strategy),".RData"))
    resulta[[2]]<- model
    rm(model)
    names(resulta)<- c("dep.temp", "indep.temp")
    eval(parse(text = paste0("resulta.ps",ms,".pt",mt,"<-resulta") ))
    rm(resulta)
  }
}

rm(list = c("ms","mt"))
################################################################################
## Tables                                                                     ##
################################################################################
########################
## Table 2. Model selection criteria, DIC, WAIC and LS, for multivariate models 
########################
criterios<- c()
for(ms in c(1:2)){
  for(mt in c(1:2)){
    eval(parse(text= paste0("resulta<- resulta.ps", ms, ".pt", mt)  ))
    aux.tabla<- matrix(NA, nrow=length(resulta), ncol=5)
    for(j in 1:length(resulta)){
      aux<- resulta[[j]]
      aux.tabla[j,]<-c(aux$dic$mean.deviance, aux$dic$p.eff, aux$dic$dic, 
                       aux$waic$waic, -mean(log(aux$cpo$cpo),na.rm=T))
      rm(aux)
    }
    colnames(aux.tabla)<-c("Mean Post D","pD","DIC","WAIC","LS")
    aux.tabla<- as.data.frame(aux.tabla)
    aux.tabla$inter<- c("Type II")
    aux.tabla$model<- names(resulta)
    aux.tabla$ps<- ms
    aux.tabla$pt<- mt
    criterios<-rbind(criterios, aux.tabla)
    rm(list=c("aux.tabla", "resulta","j"))
  }
}
rm(list = c("ms","mt"))
criterios<- criterios[order(criterios$inter, criterios$model), ]
criterios<- criterios[,c("inter","model","ps","pt","Mean Post D","pD","DIC","WAIC","LS")]
criterios

#######################
## Table 3. Estimated correlations (posterior medians and 95% credible intervals)
##          between the spatial P-spline coefficients (below main diagonal) and
##          between the temporal P-spline coefficients (above main diagonal)
#######################
model<- resulta.ps1.pt1[[1]]
aux.name<-c() 
for(r.c in 2:nc){ aux<- r.c:nc; aux.name<- c(aux.name,aux) }
rm(list = c("aux","r.c"))

## Estimated correlations between the spatial P-spline coefficients
corre.s<- matrix(NA, nrow = (nc*(nc-1)/2), ncol=3)
colnames(corre.s)<- c("median","q.025", "q.975")
for(l in c(1:(nc*(nc-1)/2))){
  eval(parse( text=  paste0("marg<-inla.tmarginal(fun= function(x){(2* exp(x))/(1+exp(x)) -1}, model$marginals.hyperpar$`Theta",nc+l," for idx`)") ))
  corre.s[l,]<- c(round(inla.qmarginal(c(0.5, 0.025, 0.975), marg), 3))
  rm(marg)
}
corre.s<- as.data.frame(corre.s)
corre.s$param<- paste0("corre.s.", rep(1:(nc-1),times=c((nc-1):1)),".",aux.name)
corre.s<- corre.s[c("param","median", "q.025", "q.975" )]
corre.s


## Estimated correlations between the temporal P-spline coefficients
corre.t<- matrix(NA, nrow = (nc*(nc-1)/2), ncol=3)
colnames(corre.t)<- c("median","q.025", "q.975")
for(l in c(1:(nc*(nc-1)/2))){
  eval(parse( text=  paste0("marg<-inla.tmarginal(fun= function(x){(2* exp(x))/(1+exp(x)) -1}, model$marginals.hyperpar$`Theta",nc+l," for idy`)") ))
  corre.t[l,]<- c(round(inla.qmarginal(c(0.5, 0.025, 0.975), marg), 3))
  rm(marg)
}
corre.t<- as.data.frame(corre.t)
corre.t$param<- paste0("corre.t", rep(1:(nc-1),times=c((nc-1):1)),".",aux.name)
corre.t<- corre.t[c("param","median", "q.025", "q.975" )]
corre.t


## rm
rm(list = c("model","aux.name","l"))
################################################################################
## Figures                                                                    ##
################################################################################
########################
## Figure 4.  Posterior median of the district-specific spatial risk for 
##            rapes (top left), assault (top right), cruelty (bottom left), 
##            and kidnapping (bottom right)
########################
model<- resulta.ps1.pt1[[1]]
marg<- lapply(model$marginals.lincomb.derived[1:(n*nc)], inla.tmarginal, 
              fun=function(x){exp(x)})
spatial<- c(unlist(lapply(marg, function(x){inla.qmarginal(0.5,x)} ))) 

carto_use<- carto
spatial_matrix<- matrix(spatial,nrow=n, ncol=nc)
colnames(spatial_matrix)<- paste0("spatial_",crimes[1:nc])
spatial_data_frame <- data.frame(dist=carto$dist, spatial_matrix)
attr(carto_use, "data") = data.frame(attr(carto_use,"data"), spatial_data_frame)

paleta<- brewer.pal(8,"YlOrRd")
minimo<- min(spatial_data_frame[,2:(1+nc)])-0.01
maximo<- max(spatial_data_frame[,2:(1+nc)])+0.01
values<- c(round(seq(minimo, 1, length.out = 5),2), 
           round(seq(1, maximo, length.out = 5),2)[-1] )

Rates.spatial<- list()
for(j in 1:nc){
  Rates.spatial[[j]]<- tm_shape(carto_use) +
    tm_polygons(col=paste("spatial_",crimes[j],sep=""), palette=paleta, title="",
                legend.show=T, legend.reverse=T, style="fixed", breaks=values, 
                interval.closure="left") +
    tm_layout(main.title="", main.title.position="center", legend.text.size=1,
              panel.labels=crime_names[j], panel.show=T, panel.label.size=1.2, 
              panel.label.bg.color="lightskyblue3", bg.color="aliceblue", 
              inner.margins=c(.04,.03, .02, .01), earth.boundary = TRUE,
              space.color="grey90")
  }

## pdf 
pdf("figures/Figure_4.pdf", height=10, width=10, onefile=FALSE)
pushViewport(viewport(layout=grid.layout(2,2)))
print(Rates.spatial[[1]], vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(Rates.spatial[[2]], vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(Rates.spatial[[3]], vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(Rates.spatial[[4]], vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
dev.off()

## rm
rm(list = c("model","marg","spatial","carto_use","spatial_matrix","paleta","j",
            "spatial_data_frame","minimo","maximo","values","Rates.spatial"))

########################
## Figure 5.  Temporal pattern of incidence risks for rape, assault, cruelty,
##            and kidnapping
########################
model<- resulta.ps1.pt1[[1]]
marg<- lapply(model$marginals.lincomb.derived[n*nc+1:(t*nc)],inla.tmarginal, 
              fun=function(x){exp(x)})
temp<- matrix(unlist(lapply(marg, function(x){inla.qmarginal(0.5,x)})),
              ncol=nc, byrow = FALSE)
q1<- matrix(unlist(lapply(marg, function(x){inla.qmarginal(0.025,x)})),
            ncol=nc, byrow = FALSE)
q2<- matrix(unlist(lapply(marg, function(x){inla.qmarginal(0.975,x)})),
            ncol=nc, byrow = FALSE)

df.final<- df.media<- c()
for(i in 1:nc){
  aux<-data.frame(X.Vec=c(x, tail(x, 1), rev(x), x[1]),
                  Y.Vec=c(q1[,i], tail(q2[,i], 1), rev(q2[,i]), q1[1,i]), 
                  crime=crime_names[i])
  df.final<- rbind(df.final,aux)
  rm(aux)
  aux<- data.frame(media=temp[,i], year=t.from:t.to, crime=crime_names[i])
  df.media<- rbind(df.media,aux)
  rm(aux)
}
rm(list = c("model","marg","temp","q1","q2","i"))

## pdf
for(i in 1:nc){
  p<- ggplot(data=df.final[df.final$crime==crime_names[i],]) +
    geom_polygon(data=df.final[df.final$crime==crime_names[i],],
                 aes(X.Vec, Y.Vec, fill=crime, group=crime), alpha = 0.7) +
    geom_line(data=df.media[df.media$crime==crime_names[i],],
              aes(year,media, group=crime), color='black',lwd=.6) +
    facet_wrap(. ~ crime, ncol = 2, scales = "fixed") +
    xlab("Year") +
    ylab(expression(exp(f[j](x[t])))) +
    ylim(c(0.6,2.25)) +
    geom_hline(yintercept=1, linetype="dashed", color = "black") +
    theme(legend.position = c("none"),
          strip.text.x = element_text(size = 13)) +
    scale_fill_manual(values = selected_colors[i])
  eval(parse(text =  paste0("p.",i,"<-p") ))
  rm(p)
}

## pdf
pdf("figures/Figure_5.pdf", height=8, width=10, onefile=FALSE)
ggpubr::ggarrange(p.1,p.2,p.3,p.4, nrow = 2, ncol = 2)
dev.off()

## rm
rm(list = c("i",paste0("p.",1:nc)))

########################
## Figure 6.  Relative risk evolution (posterior median of $R_{ijt}$) for three 
##            selected districts: Aurangabad, Garhchiroli, and Greater Bombay
########################
id.area<- ID[ID$dist %in% c("Aurangabad","Garhchiroli","Greater Bombay"),"ID_area"]

model<- resulta.ps1.pt1[[1]]
d.data_frame_use<- d.data_frame
d.data_frame_use$risks.mean<- model$summary.fitted.values$`0.5quant`[1:(nc*t*n)]
d.data_frame_use$risks.q025<- model$summary.fitted.values$`0.025quant`[1:(nc*t*n)]
d.data_frame_use$risks.q975<- model$summary.fitted.values$`0.975quant`[1:(nc*t*n)]

df.final<- df.media<-c()
for(i in id.area){
  df<- c()
  df.med<- c()
  for(j in 1:length(crimes)){
    X.Vec <- c(x, tail(x, 1), rev(x), x[1])
    Y.Vec <- c(d.data_frame_use[d.data_frame_use$ID_area==i & d.data_frame_use$ID_disease==j, c("risks.q025")],
               tail(d.data_frame_use[d.data_frame_use$ID_area==i & d.data_frame_use$ID_disease==j, c("risks.q975")], 1), 
               rev(d.data_frame_use[d.data_frame_use$ID_area==i & d.data_frame_use$ID_disease==j, c("risks.q975")]),
               d.data_frame_use[d.data_frame_use$ID_area==i & d.data_frame_use$ID_disease==j, c("risks.q025")][1])
    df.aux<- data.frame(X.Vec=X.Vec,
                        Y.Vec=Y.Vec,
                        ID=ID[ID$ID_area==i,"ID_area"],
                        District= ID[ID$ID_area==i,"dist"],
                        crime=crime_names[j])
    df<- rbind(df,df.aux)
    df.med.aux<- data.frame(year=t.from:t.to,
                            media=d.data_frame_use[d.data_frame_use$ID_area==i & d.data_frame_use$ID_disease==j, c("risks.mean")],
                            ID=ID[ID$ID_area==i,"ID_area"],
                            District= ID[ID$ID_area==i,"dist"],
                            crime=crime_names[j]
    )
    df.med<- rbind(df.med, df.med.aux)
    rm(list = c("X.Vec","Y.Vec","df.aux","df.med.aux"))
  }
  rm(j)
  df.final<- rbind(df.final, df)
  df.media<- rbind(df.media, df.med)
  rm(list=c("df","df.med"))
}
df.final$label<- paste0(df.final$District,"_",df.final$crime)
df.media$label<- paste0(df.media$District,"_",df.media$crime)

## pdf
p<- ggplot(data=df.final) +
  geom_polygon(data=df.final,
               aes(X.Vec,Y.Vec, fill=crime, group=District),alpha = 0.7) +
  geom_line(data=df.media,
            aes(year,media, group=District), color='black',lwd=.6) +
  facet_wrap(District ~  crime , ncol = 4, scales = "fixed",dir="h") +
  xlab("Year") +
  ylab(expression(R[itj])) +
  geom_hline(yintercept=1, linetype="dashed", color = "black") +
  theme(legend.position = c("top"),
        legend.title = element_blank(),
        legend.text = element_text(color = "black", size = 13),
        strip.text.x = element_text(size = 13)
        ) +
  scale_fill_manual(values = selected_colors)

## pdf
pdf("figures/Figure_6.pdf", height=8, width=14, onefile=FALSE)
print(p)
dev.off()

## rm
rm(list = c("id.area","model","d.data_frame_use","df.final","df.media","i","p"))

########################
## Figure 7.  Functional boxplots of the final relative risk trends for rape 
##            (top left), assault (top right), cruelty (bottom left), and 
##            kidnapping (bottom right)
########################
library(fda)
## 'fbplot_modified.R' is a modification of the 'fbplot' function of the 'fda' library.
source('fbplot/fbplot_modified.R')

risk.array<- array(resulta.ps1.pt1[[1]]$summary.fitted.values$`0.5quant`[1:(nc*t*n)],
                   dim=c(n,t,nc))
inf<- min(risk.array)-0.01
top<- round(max(risk.array)+0.51,2)

## without arrows
pdf("figures/Figure_7.pdf", height=8, width=10, onefile=FALSE)
par(mfrow=c(2,2))
for(i in 1:nc){
  risk.mat<- risk.array[,,i]
  t.risk.mat<- as.data.frame(t(risk.mat))
  fbplot(t.risk.mat,method='MBD',ylim=c(inf,top), fullout=T, 
         barcol="slateblue4", outliercol="red", xlab="",ylab="", yaxt="n",
         col=rgb(255,48,48,alpha=200, maxColorValue=255))
  axis(side=2, at=seq(0,4,1), labels = TRUE, lty = 0)
  title(ylab=expression(R[ijt]), xlab="Year", cex.lab=1.3, line=2.0)
  rect(par("usr")[1],top-0.5,par("usr")[2],par("usr")[4], 
       col=" gray83", border = "gray83")
  legend(7,top-0.5, crime_names[i], box.col = "transparent", 
         bg = "transparent", adj = 0.2, xjust = 0.5, yjust = 0.1, cex = 1.2)
  rm(list = c("risk.mat","t.risk.mat"))
}
dev.off()

rm(list = c("risk.array","inf","top","i"))

## with arrows
pdf("figures/Figure_7.pdf", height=8, width=10, onefile=FALSE)
par(mfrow=c(2,2))

fbplot(as.data.frame(t( risk.array[,,1])), method='MBD', ylim=c(inf,top), fullout=T, 
       barcol="slateblue4", outliercol="red",xlab="",ylab="", yaxt="n", 
       col=rgb(255,48,48,alpha=200, maxColorValue=255)) 
axis(side=2, at=seq(0,4,1), labels = TRUE, lty = 0)
title(ylab=expression(R[ijt]), xlab="Year", cex.lab=1.3, line=2.0)
rect(par("usr")[1],top-0.5,par("usr")[2],par("usr")[4],
     col = "gray83", border = "gray83")
legend(7,top-0.5, crime_names[i], box.col = "transparent",
       bg = "transparent", adj = 0.2, xjust = 0.5, yjust = 0.1, cex = 1.2)
arrows(10,2.5, 10,1.9, length=0.2, col=1, angle = 30)
text(10,2.7,ID$dist[32],col=1)

fbplot(as.data.frame(t( risk.array[,,2])), method='MBD',ylim=c(inf,top), fullout=T, 
       barcol="slateblue4", outliercol="red",xlab="",ylab="", yaxt="n",
       col=rgb(255,48,48,alpha=200, maxColorValue=255))
axis(side=2, at=seq(0,4,1), labels = TRUE, lty = 0)
title(ylab=expression(R[ijt]), xlab="Year", cex.lab=1.3, line=2.0)
rect(par("usr")[1],top-0.5,par("usr")[2],par("usr")[4],
     col = "gray83", border = "gray83")
legend(7,top-0.5, crime_names[2], box.col = "transparent",
       bg = "transparent", adj = 0.2, xjust = 0.4, yjust = 0.1, cex = 1.2)
arrows(8,3.1, 8,2.5, length=0.2, col=1, angle = 35)
text(8,3.3,ID$dist[32],col=1)

fbplot(as.data.frame(t( risk.array[,,3])), method='MBD',ylim=c(inf,top), fullout=T, 
       barcol="slateblue4", outliercol="red",xlab="",ylab="", yaxt="n",
       col=rgb(255,48,48,alpha=200, maxColorValue=255))
axis(side=2, at=seq(0,4,1), labels = TRUE, lty = 0)
title(ylab=expression(R[ijt]), xlab="Year", cex.lab=1.3, line=2.0)
rect(par("usr")[1],top-0.5,par("usr")[2],par("usr")[4],
     col = "gray83", border = "gray83")
legend(7,top-0.5, crime_names[3], box.col = "transparent",
       bg = "transparent", adj = 0.2, xjust = 0.38, yjust = 0.1, cex = 1.2)

fbplot(as.data.frame(t( risk.array[,,4])), method='MBD',ylim=c(inf,top), fullout=T, 
       barcol="slateblue4", outliercol="red", xlab="",ylab="", yaxt="n",
       col=rgb(255,48,48,alpha=200, maxColorValue=255))
axis(side=2, at=seq(0,4,1), labels = TRUE, lty = 0)
title(ylab=expression(R[ijt]), xlab="Year", cex.lab=1.3, line=2.0)
rect(par("usr")[1],top-0.5,par("usr")[2],par("usr")[4],
     col = "gray83", border = "gray83")
legend(7,top-0.5, crime_names[4], box.col = "transparent", 
       bg = "transparent", adj = 0.2, xjust = 0.38, yjust = 0.1, cex = 1.2)
dev.off()

rm(list = c("risk.array","inf","top","i"))

########################
## Figure A1. Exceedance probabilities for rapes (top left), assault (top right),
##            cruelty (bottom left), and kidnapping (bottom right).
########################
model<- resulta.ps1.pt1[[1]]
marg<- lapply(model$marginals.lincomb.derived[1:(n*nc)], inla.tmarginal,
              fun=function(x){exp(x)})
prob<- c(unlist(lapply(marg, function(x){1-inla.pmarginal(1,x)})))

carto_use<- carto
prob_matrix<- matrix(prob,nrow=n, ncol=nc)
colnames(prob_matrix)<- paste0("prob_",crimes[1:nc])
spatial_data_frame <- data.frame(dist=carto$dist,prob_matrix)
attr(carto_use, "data") = data.frame(attr(carto_use,"data"), spatial_data_frame)

paleta2 <- brewer.pal(5,"PuBu")
values2 <- c(0,0.1,0.2,0.8,0.9,1)

Prob.spatial<- list()
for(j in 1:nc){
  Prob.spatial[[j]]<- tm_shape(carto_use) +
    tm_polygons(col=paste("prob_",crimes[j],sep=""), palette=paleta2, title="",
                legend.show=T, legend.reverse=T, style="fixed", breaks=values2,
                interval.closure="left") +
    tm_layout(main.title="", main.title.position="center", legend.text.size=1,
              panel.labels=crime_names[j], panel.show=T, panel.label.size=1.2, 
              panel.label.bg.color="lightskyblue3", bg.color="aliceblue",
              inner.margins=c(.04,.03, .02, .01), earth.boundary = TRUE,
              space.color="grey90")
  }

## pdf
pdf("figures/Figure_A1.pdf", height=10, width=10, onefile=FALSE)
pushViewport(viewport(layout=grid.layout(2,2)))
print(Prob.spatial[[1]], vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(Prob.spatial[[2]], vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(Prob.spatial[[3]], vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(Prob.spatial[[4]], vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
dev.off()

## rm
rm(list = c("model","marg","prob","carto_use","prob_matrix","paleta2","j",
            "spatial_data_frame","values2","Prob.spatial"))

########################
## Figure A2. Specific temporal trends (posterior median of exp(delta_ijt) for
##            three selected districts: Aurangabad, Garhchiroli, and Greater Bombay
########################
id.area<- ID[ID$dist %in% c("Aurangabad","Garhchiroli","Greater Bombay"),"ID_area"]

model<- resulta.ps1.pt1[[1]]
q.5<- q.025<- q.975<- c()
for(i in 1:nc){
  eval(parse(text = paste0("aux<- lapply(model$marginals.random$idxy.",i,",
                           inla.tmarginal, fun=function(x){exp(x)})") ))
  q.5<- c(q.5, unlist(lapply(aux, function(x){inla.qmarginal(0.5,x)} )))
  q.025<- c(q.025, unlist(lapply(aux, function(x){inla.qmarginal(0.025,x)})))
  q.975<- c(q.975, unlist(lapply(aux, function(x){inla.qmarginal(0.975,x)})))
  rm(aux)
}
rm(i)
d.data_frame_use<- d.data_frame
d.data_frame_use$risks.mean<- q.5
d.data_frame_use$risks.q025<- q.025
d.data_frame_use$risks.q975<- q.975
rm(list = c("q.5","q.025","q.975"))

df.final<- df.media<-c()
for(i in id.area){
  df<- df.med<- c()
  for(j in 1:length(crimes)){
    X.Vec <- c(x, tail(x, 1), rev(x), x[1])
    Y.Vec <- c(d.data_frame_use[d.data_frame_use$ID_area==i & 
                                  d.data_frame_use$ID_disease==j, c("risks.q025")],
               tail(d.data_frame_use[d.data_frame_use$ID_area==i & 
                                       d.data_frame_use$ID_disease==j, c("risks.q975")], 1), 
               rev(d.data_frame_use[d.data_frame_use$ID_area==i &
                                      d.data_frame_use$ID_disease==j, c("risks.q975")]),
               d.data_frame_use[d.data_frame_use$ID_area==i &
                                  d.data_frame_use$ID_disease==j, c("risks.q025")][1])
    df.aux<- data.frame(X.Vec=X.Vec, Y.Vec=Y.Vec, ID=ID[ID$ID_area==i,"ID_area"],
                        District= ID[ID$ID_area==i,"dist"], crime=crime_names[j])
    df<- rbind(df,df.aux)
    df.med.aux<- data.frame(year=t.from:t.to,ID=ID[ID$ID_area==i,"ID_area"],
                          media=d.data_frame_use[d.data_frame_use$ID_area==i &
                                      d.data_frame_use$ID_disease==j, c("risks.mean")],
                            District= ID[ID$ID_area==i,"dist"],crime=crime_names[j]
                          )
    df.med<- rbind(df.med, df.med.aux)
    rm(list = c("X.Vec","Y.Vec","df.aux","df.med.aux"))
  }
  df.final<- rbind(df.final, df)
  df.media<- rbind(df.media, df.med)
  rm(list=c("df","df.med"))
}
rm(list = c("i","j"))

## pdf
p<- ggplot(data=df.final) +
  geom_polygon(data=df.final,
               aes(X.Vec, Y.Vec, fill=crime, group=District), alpha = 0.7) + # , fill=selected_colors
  geom_line(data=df.media,
            aes(year,media, group=District), color='black',lwd=.6) +
  facet_wrap(District ~  crime , ncol = 4, scales = "fixed",dir="h") +
  xlab("Year") +
  ylab(expression(exp(delta[ijt]))) +
  geom_hline(yintercept=1, linetype="dashed", color = "black") +
  theme(legend.position = c("top"),
        legend.justification = c("center", "top"),
        # legend.box.just = "center",
        # legend.margin = margin(4, 4, 4, 4),
        legend.title = element_blank(),
        legend.text = element_text(color = "black", size = 12),
        strip.text.x = element_text(size = 13)
  ) +
  scale_fill_manual(values = selected_colors)

pdf("figures/Figure_A2.pdf", height=8, width=14, onefile=TRUE)
print(p)
dev.off()

## rm
rm(list = c("id.area","model","d.data_frame_use","df.final","df.media","p"))

########################
## Figure A3: A6. Map of estimated incidence risks for rape/assault/cruelty/kidnapping (top) 
##                and posterior probabilities that the relative risk is greater than one 
##                in Maharashtra between 2001 and 2013.
########################
model<- resulta.ps1.pt1[[1]]
d.data_frame_use<-d.data_frame
d.data_frame_use$risks<- model$summary.fitted.values$`0.5quant`[1:(nc*t*n)]
d.data_frame_use$prp<- 1-model$summary.linear.predictor[,"0 cdf"][1:(nc*t*n)]

risks<- list()
prob<- list()
for(j in 1:nc){
  risks[[j]]<- data.frame(matrix(d.data_frame_use[d.data_frame_use$ID_disease==j,c("risks")], nrow=n, ncol=t, byrow=FALSE))
  colnames(risks[[j]])<- paste0("risk_", crimes[j],"_",t.from:t.to)
  
  prob[[j]]<- data.frame(matrix(d.data_frame_use[d.data_frame_use$ID_disease==j,c("prp")], nrow=n, ncol=t, byrow=FALSE))
  colnames(prob[[j]])<- paste0("prob_", crimes[j],"_",t.from:t.to)
}

carto_use<- carto
risks_data_frame <- data.frame(dist=carto$dist, risks, prob)
attr(carto_use, "data") = data.frame(attr(carto_use,"data"), risks_data_frame)

## risks
paleta <- RColorBrewer::brewer.pal(8,"YlOrRd")
minimo<- min(d.data_frame$risks)-0.01
maximo<- max(d.data_frame$risks)+0.01
values<- c(round(seq(minimo, 1, length.out = 5),2), round(seq(1, maximo, length.out = 5),2)[-1] )
Rates.Risks<- list()
for(j in 1:nc){
  Rates.Risks[[j]]<- tm_shape(carto_use) + 
    tm_polygons(col=paste("risk_",crimes[j],"_",seq(t.from,t.to), sep=""),
                palette=paleta, title="", legend.show=T, legend.reverse=T, 
                style="fixed", breaks=values, interval.closure="left") +
    tm_layout(main.title=crime_names[j], main.title.position="center",
              legend.text.size=1, panel.show=T, panel.label.size=2, 
              panel.labels=as.character(round(seq(t.from,t.to,length.out=t))),
              panel.label.bg.color="lightskyblue3", bg.color="aliceblue",
              inner.margins=c(.04,.03,.02,.01), earth.boundary = TRUE, 
              space.color="grey90", legend.outside=T, legend.outside.position="right",
              legend.frame=F, legend.outside.size=0.25, outer.margins=c(0.02,0.01,0.02,0.01)) +
    tm_facets(ncol=5, nrow=3)
}

## prob risks
paleta2 <- brewer.pal(5,"PuBu")
values2 <- c(0,0.1,0.2,0.8,0.9,1)
Prob.Risks<- list()
for(j in 1:nc){
  Prob.Risks[[j]]<- tm_shape(carto_use) + 
    tm_polygons(col=paste("prob_",crimes[j],"_",seq(t.from,t.to), sep=""),
                palette=paleta2, title="",legend.show=T, legend.reverse=T,
                style="fixed", breaks=values2, interval.closure="left") +
    tm_layout(main.title=crime_names[j] , main.title.position="center",
              legend.text.size=1, panel.show=T, panel.label.size=2,
              panel.labels=as.character(round(seq(t.from,t.to,length.out=t))),
              panel.label.bg.color="lightskyblue3", bg.color="aliceblue",
              inner.margins=c(.04,.03,.02,.01), earth.boundary = TRUE, space.color="grey90",
              legend.outside=T, legend.outside.position="right", legend.frame=F,
              legend.outside.size=0.25, outer.margins=c(0.02,0.01,0.02,0.01)) +
    tm_facets(ncol=5, nrow=3)
}

## pdf
for(i in 1:nc){
  pdf(paste0("figures/Figure_A",i+2,".pdf"), height=10, width=7.5, onefile=TRUE)
  pushViewport(viewport(layout=grid.layout(2,1)))
  print(Rates.Risks[[i]], vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(Prob.Risks[[i]],  vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
  dev.off()
}

## rm
rm(list = c("model","d.data_frame_use","risks","prob","carto_use",
            "risks_data_frame","paleta","minimo","maximo","values",
            "Rates.Risks","paleta2","values2","Prob.Risks","j"))

################################################################################
################################################################################
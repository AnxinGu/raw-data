if (T) {
  dir.create("PDFs")
  dir.create("PDFs/图片")
  dir.create("files")
  dir.create("files/文件")
  dir.create("origin_datas")
  dir.create("origin_datas/TCGA")
  dir.create("origin_datas/GEO")
}

## 无虚线
ggplotKMCox=function(dat,title='Groups',labs=NULL,add_text=NULL,pal=NULL){
  library(ggplot2)
  library(ggsci)
  colnames(dat)=c('time','status','groups')
  #sdf<-survdiff(Surv(time,status) ~ groups,data=dat)
  #print((sdf))
  #summary(sdf)
  #p<-pchisq(sdf$chisq,length(sdf$n)-1,lower.tail=FALSE)
  sf<-survfit(Surv(time,status) ~ groups,data=dat)
  if(is.null(pal)){
    surv=survminer::ggsurvplot(sf, data = dat
                               , palette =pal_lancet()(8)[c(1,2)]
                               ,pval = TRUE, surv.median.line='hv'
                               ,conf.int = T
                               ,linetype = "strata"
                               ,conf.int.style ='step'
                               , pval.coord=c(0, 0.2), #Add p-value 
                               risk.table = TRUE, 
                               legend.title = title
                               ,legend.labs = labs
    )
    
  }else{
    if(pal=='npg'){
      surv=survminer::ggsurvplot(sf, data = dat
                                 , palette = pal_npg()(9)[c(2,1,3,4:9)] #jco palette 
                                 ,pval = TRUE, surv.median.line='hv'
                                 # ,conf.int = T
                                 # ,linetype = "strata"
                                 ,conf.int.style ='step'
                                 , pval.coord=c(0, 0.2), #Add p-value 
                                 risk.table = TRUE, 
                                 legend.title = title
                                 ,legend.labs = labs
      )
    }else{
      surv=survminer::ggsurvplot(sf, data = dat
                                 , palette = pal #jco palette 
                                 ,pval = TRUE, surv.median.line='hv'
                                 # ,conf.int = T
                                 # ,linetype = "strata"
                                 ,conf.int.style ='step'
                                 , pval.coord=c(0, 0.2), #Add p-value 
                                 risk.table = TRUE, 
                                 legend.title = title
                                 ,legend.labs = labs
      )
    }
    
  }
  
  p1=surv$plot+theme_bw()+theme(axis.text.y=element_text(family="Times",face="plain")
                                ,axis.text.x=element_blank()
                                ,axis.title.x=element_blank()
                                ,plot.margin=unit(c(0.2, 0.2, 0, 0.1), "inches")
                                #,axis.title.y=element_blank()
                                ,legend.position=c(1,1), legend.justification=c(1,1)
                                ,legend.background = element_rect(fill = NA, colour = NA)
                                ,legend.title = element_text(family="Times",face="plain")
                                ,legend.text = element_text(family="Times",face="plain"))
  #p1=p1+text()
  #tms=data.frame(Group=tms.gp,value=tms.tps,Attribute=rep(data_m[1,1],length(tms.gp))
  #               ,ymax=rep(max(ylim),length(tms.gp)))
  #p4=p4+geom_text(data=tms,aes(x=Group, y=ymax, label=value),color="yellow")
  if(!is.null(add_text)){
    text.tb=surv$data.survplot[1,]
    text.tb[1,1]=0
    text.tb[1,5]=0
    text.tb$Text=add_text
    p1=p1+geom_text(data=text.tb,aes(x=time, y=surv, label=Text),color="yellow",hjust =0)
  }
  
  p2=surv$table+theme_bw()+theme(axis.text.y=element_text(family="Times",face="plain")
                                 #,axis.text.x=element_blank()
                                 #,axis.title.x=element_blank()
                                 #,axis.title.y=element_blank()
                                 ,plot.margin=unit(c(0, 0.2, 0.2, 0.1), "inches")
                                 ,plot.title=element_blank()
                                 ,legend.position=c(1,1), legend.justification=c(1,1)
                                 #,legend.background = element_rect(fill = NA, colour = NA)
                                 ,legend.title = element_text(family="Times",face="plain")
                                 ,legend.text = element_text(family="Times",face="plain"))
  p1=surv$plot
  p2=surv$table
  g2=ggpubr::ggarrange(p1,p2, ncol = 1, nrow = 2,heights = c(0.9,0.3),align = "v")
  return(g2)
}

################ timeROC
ggplotTimeROC_use=function(time,status,score,mks=c(1,3,5)){
  #time=g.os
  #status=g.ev
  #score=as.numeric(cpm.score)
  #cx=coxRun(data.frame(time,status,score))
  #if(cx[1]<=1){
  #  score=-1*score
  #}
  roc.tm=mg_surv_pROC(time,status,score,mks)
  print('roc.tm')
  print((roc.tm))
  library(survival)
  library(ggplot2)
  mks=mg_predict_time_ymd(time,mks)
  print(mks)  
  ROC.DSST=timeROC::timeROC(T=time,
                            delta=status
                            ,marker=score,
                            cause=1,weighting="marginal",
                            times=mks,
                            iid=TRUE)
  print(ROC.DSST)
  mks=mks[which(!is.na(ROC.DSST$AUC)&ROC.DSST$AUC>0)]
  print(mks)
  if(length(mks)>0){
    if(max(ROC.DSST$AUC)<0.5){
      score=-1*score
    }
    ROC.DSST=timeROC::timeROC(T=time,
                              delta=status
                              ,marker=score,
                              cause=1,weighting="marginal",
                              times=mks,
                              iid=TRUE)
    print(ROC.DSST$times)
    if(max(ROC.DSST$times)<20){
      lb=paste0(ROC.DSST$times,'-Years')
    }else if(max(ROC.DSST$times)<365){
      lb=paste0(round(ROC.DSST$times/12,0),'-Years')
    }else{
      lb=paste0(round(ROC.DSST$times/365,0),'-Years')
    }
    
    lbs=paste0(lb,',AUC=',round(ROC.DSST$AUC,2),',95%CI(',paste0(round(confint(ROC.DSST,level = 0.95,na.rm=T)$CI_AUC[,1]/100,2),'-',
                                                                 round(confint(ROC.DSST,level = 0.95,na.rm=T)$CI_AUC[,2]/100,2)),')')
    #roc.tm=ROC.DSST$times[which(ROC.DSST$times>0)]
    
    #p.dat=rbind()
    #for(i in which(ROC.DSST$times>0)){
    #los=lowess(ROC.DSST$FP[,i], y=ROC.DSST$TP[,i], f = 1/3, iter = 100)
    #los$x=c(0,los$x,1)
    #los$y=c(0,los$y,1)
    # p.dat=rbind(p.dat,data.frame(los$x, y=los$y,rep(lbs[i],length(los$y)),stringsAsFactors = F))
    #}
    
    p.dat=rbind()
    print(length(roc.tm))
    for(i in 1:length(roc.tm)){
      #print(i)
      r1=roc.tm[[i]]
      x1=1-r1$specificities
      y1=r1$sensitivities
      #print(cbind(1-r1$specificities,r1$sensitivities))
      nx1=unique(x1)
      ny1=c()
      for(x in unique(x1)){
        x.inds=which(x1==x)
        if(length(x.inds)>0&x<0.5){
          ny1=c(ny1,min(y1[x.inds]))
        }else if(length(x.inds)>0){
          ny1=c(ny1,max(y1[x.inds]))
        }else{
          ny1=c(ny1,y1[x.inds][1])
        }
      }
      #print(cbind(nx1,ny1))
      p.dat=rbind(p.dat,data.frame(x=nx1, y=ny1,rep(lbs[i],length(nx1)),stringsAsFactors = F))
    }
    colnames(p.dat)=c('V1','V2','Type')
    p.dat=as.data.frame(p.dat)
    
    p1=ggplot(p.dat, aes(x=V1,y=V2, fill=Type))
    p1=p1+geom_line(aes(colour=Type),lwd=1.1)+theme_bw()+xlab('False positive fraction')+ylab('True positive fraction')
    # p1=p1+stat_smooth(aes(colour=Type),se = FALSE, size = 1)+theme_bw()+xlab('False positive fraction')+ylab('True positive fraction')
    
    p1=p1+scale_colour_manual(values = c(pal_npg(alpha =0.8)(9)[c(8,4,3,5)]))
    
    p1=p1+theme(axis.text.y=element_text(family="Times",face="plain"),axis.text.x=element_text(family="Times",face="plain")
                ,axis.title.x=element_text(family="Times",face="plain"),axis.title.y=element_text(family="Times",face="plain")
                ,plot.title=element_blank()
                ,plot.margin=unit(c(0.1, 0.1, 0.1, 0.1), "inches")
                ,legend.position=c(1,0)
                ,legend.justification=c(1,0)
                ,legend.background = element_rect(fill = NA, colour = NA)
                ,legend.title = element_text(family="Times",face="plain")
                ,legend.text = element_text(family="Times",face="plain"))
    return(p1)
  }else{
    return(mg_getplot_bank('No data plot by ROC!'))
  }
}

plotRiskScoreModel_use=function(riskScore,dat,time,event,cutoff,hetTitle='z-score of expression',hetColor=c('green','black','red')){
  srt.inds=order(riskScore)
  dat=dat[srt.inds,]
  time=time[srt.inds]
  event=event[srt.inds]
  riskScore=riskScore[srt.inds]
  library(ggplot2)
  dt1=data.frame(V1=1:length(riskScore),V2=riskScore,RiskType=ifelse(riskScore>cutoff,'High','Low')) 
  # p1=ggplot(dt1, aes(x = V1, y = V2, colour = RiskType,fill=RiskType)) +geom_bar(stat = 'identity', position = 'dodge')+ggsci::scale_fill_npg()+theme_bw()
  # p1=ggplot(dt1, aes(x = V1, y = V2, colour = RiskType,fill=RiskType)) +geom_point(stat = 'identity', position = 'identity')+ggsci::scale_fill_npg()+theme_bw()
  p1=ggplot(dt1, aes(x = V1, y = V2, colour = RiskType,fill=RiskType)) +geom_point(stat = 'identity', position = 'identity')+scale_colour_manual(values = ggsci::pal_lancet()(9)[c(2,1)],aesthetics = c("colour", "fill"))+theme_bw()
  p1=p1+geom_hline(aes(yintercept=cutoff),colour='black',linetype="dashed")
  p1=p1+ylab('RiskScore')+theme(axis.text.y=element_text(family="Times",face="plain"),axis.text.x=element_blank()
                                ,axis.title.x=element_blank(),legend.position=c(1,0), legend.justification=c(1,0)
                                ,legend.background = element_rect(fill = NA, colour = NA)
                                ,plot.margin=unit(c(0.1, 0.1, 0, 0.1), "inches")
                                ,legend.title = element_text(family="Times",face="plain")
                                ,legend.text = element_text(family="Times",face="plain"))
  
  dt2=data.frame(V1=1:length(riskScore),V2=time,Status=ifelse(event==1,'Dead','Alive'))  
  p2=ggplot(dt2, aes(x = V1, y = V2, colour = Status,shape =Status)) +geom_point()+ggsci::scale_fill_npg()+theme_bw()
  p2=p2+ylab('Time')+theme(axis.text.y=element_text(family="Times",face="plain"),axis.text.x=element_blank()
                           ,axis.title.x=element_blank(),legend.position=c(1,1), legend.justification=c(1,1)
                           ,legend.background = element_rect(fill = NA, colour = NA)
                           ,plot.margin=unit(c(0, 0.1, 0, 0.1), "inches")
                           ,legend.title = element_text(family="Times",face="plain")
                           ,legend.text = element_text(family="Times",face="plain"))
  
  data=as.data.frame(scale(dat))
  hc.r = hclust(dist(t(data)))
  data=data[,hc.r$order]
  data$ID <- 1:nrow(dat)
  #colnames(data)
  data_m <- reshape2::melt(data, id.vars=c("ID"))
  colnames(data_m)=c('ID','V1','V2')
  data_m$V2[which(data_m$V2>mean(data_m$V2)+3*sd(data_m$V2))]=mean(data_m$V2)+3*sd(data_m$V2)
  data_m$V2[which(data_m$V2<mean(data_m$V2)-3*sd(data_m$V2))]=mean(data_m$V2)-3*sd(data_m$V2)
  
  data_m$V1=mg_str_outline(data_m$V1,isCut = T,n=50)
  #print(data_m$V1)
  #print(head(data_m))
  #data_m[1:20,]
  p3 <- ggplot(data_m, aes(x=ID,y=V1)) 
  p3 <- p3 + geom_tile(aes(fill=V2))
  p3=p3+scale_fill_gradient2(low = hetColor[1],mid=hetColor[2], high = hetColor[3])
  p3=p3+theme_bw()
  p3=p3+ labs(fill=hetTitle) 
  #p3=p3+guides(fill=guide_legend(title="New Legend Title"))
  p3=p3+xlab('Samples')
  #p3 <- p3 + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
  p3=p3+theme(axis.text.y=element_text(family="Times",face="plain")
              ,axis.text.x=element_blank()
              #,axis.title.x=element_blank()
              ,axis.title.y=element_blank()
              ,legend.position='bottom'
              #,legend.justification=c(1,1)
              #,legend.background = element_rect(fill = NA, colour = NA)
              ,plot.margin=unit(c(0, 0.1, 0.1, 0.1), "inches")
              ,legend.title = element_text(family="Times",face="plain")
              ,legend.text = element_text(family="Times",face="plain"))
  
  g1=ggpubr::ggarrange(p1,p2,p3, ncol = 1, nrow = 3,heights = c(0.5,0.5,1),align = "v")
  return(g1)
}

################ KM
plotCoxModel_Batch_use=function(riskScore,dat,time,event,cutoff,title='Groups',hetTitle='z-score of expression',hetColor=c('green','black','red'),mks=c(1,3,5)){
  g1=plotRiskScoreModel_use(riskScore,dat,time,event,cutoff,hetTitle,hetColor)
  dat=data.frame(time,event,ifelse(riskScore>cutoff,'high','low'),stringsAsFactors = F)
  dat=dat[order(dat[,3]),]
  cx=coxRun(dat = data.frame(time,event,riskScore))
  txt=paste0('HR=',round(cx[2],2),' 95CI%(',round(cx[3],2),'-',round(cx[4],2),')')
  dat[,3]=factor(dat[,3],levels = c('low','high'))
  print(unique(dat[,3]))
  g3=ggplotKMCox(dat,title=title,add_text = txt,pal = NULL)
  g2=ggplotTimeROC_use(time,event,riskScore,mks)
  g23=ggpubr::ggarrange(g2,g3, ncol = 1, nrow = 2,heights = c(1,1)
                        #,align = "hv"
                        ,labels = toupper(letters)[2:3])
  g23.row=ggpubr::ggarrange(g2,g3, ncol = 2, nrow = 1,heights = c(1,1)
                            #,align = "hv"
                            ,labels = toupper(letters)[2:3])
  
  gal=ggpubr::ggarrange(g1,g23, ncol = 2, nrow = 1,widths =  c(1.5,1)
                        #,align = "hv"
                        ,labels = c(toupper(letters)[1],''))
  return(list(gal,g23,g23.row))
}


#################
mg_violin_use=function(data,group.col=NULL,xangle=0,ylab='value',xlab='',leg.title='Group',test_method='anova',cmp_test_method='t.test',legend.pos='r',melt=F,jitter=T,ylim=NULL,show_compare=NULL,point_size=NULL){
  library(ggplot2)
  if(melt){
    data_m=data
    colnames(data_m)=c('Group','value')
  }else{
    data_m=reshape2::melt(data)
    colnames(data_m)=c('Group','value')
  }
  data_m=data_m[which(!is.na(data_m[,1])),]
  if(xangle==0){
    tx=element_text(colour="black",family="Times")
  }else{
    tx=element_text(angle=xangle,hjust = 1,colour="black",family="Times")
  }
  
  pos='right'
  if(is.null(legend.pos)){
    pos='none'
  }else if(legend.pos=='tr'){
    pos=c(1,1)
  }else if(legend.pos=='br'){
    pos=c(1,0)
  }else if(legend.pos=='tl'){
    pos=c(0,1)
  }else if(legend.pos=='bl'){
    pos=c(0,0)
  }else if(legend.pos=='t'){
    pos='top'
  }else if(legend.pos=='r'){
    pos='right'
  }else if(legend.pos=='b'){
    pos='bottom'
  }
  uni.group=unique(data_m[,1])
  ct=length(uni.group)
  if(!is.null(ylim)){
    data_m=data_m[data_m[,2]<=max(ylim),]
    data_m=data_m[which(!is.na(data_m[,1])),]
    ylim[2]=1.2*ylim[2]
  }
  p1<-ggplot(data_m,aes(x=Group,y=value))+geom_violin(alpha=0.7)
  if(ct<=10){
    # p1=p1+ggsci::scale_fill_npg(name=leg.title)
    if(!is.null(group.col)){
      p1=p1+scale_fill_manual(values = group.col)
    }else{
      p1=p1+scale_fill_manual(values = c(pal_lancet('lanonc',alpha =0.8)(9)[c(7,1,4,3,5:6)]))
    }
    
  }else if(ct<=20){
    p1=p1+ggsci::scale_fill_d3(palette = "category20",name=leg.title)
  }else if(ct<=30){
    cbPalette=c(ggsci::pal_npg("nrc", alpha = 0.6)(10),ggsci::pal_d3("category20", alpha = 0.6)(20))
    p1=p1+scale_fill_manual(values=cbPalette[1:ct])
  }else if(ct<=38){
    cbPalette=c(ggsci::pal_npg("nrc", alpha = 0.6)(10)
                ,ggsci::pal_d3("category20", alpha = 0.6)(20)
                ,ggsci::pal_nejm("default", alpha = 0.6)(8))
    p1=p1+scale_fill_manual(values=cbPalette[1:ct])
  }
  
  if(jitter){
    if(is.null(point_size)){
      p1<-p1+geom_jitter(alpha=0.3,col='black',show.legend=FALSE,width = 0.2)
    }else{
      p1<-p1+geom_jitter(alpha=0.3,col='black',show.legend=FALSE,width = 0.2,size=point_size)
    }
  }
  
  p1=p1+theme_bw()+geom_boxplot(width=0.2,aes(fill=Group),outlier.shape = NA)
  p1=p1+theme(axis.text.x=tx, #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
              axis.text.y=element_text(family="Times",face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
              axis.title.y=element_text(family="Times",face="plain"), #设置y轴标题的字体属性
              #panel.border = element_blank(),axis.line = element_line(colour = "black"), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
              legend.text=element_text(face="plain", family="Times", colour="black"  #设置图例的子标题的字体属性
              ),
              legend.title=element_text(face="plain", family="Times", colour="black" #设置图例的总标题的字体属性
              ),
              legend.justification=pos, legend.position=pos
              ,legend.background = element_rect(fill = NA, colour = NA)
              #,panel.grid.major = element_blank(),   #不显示网格线
              #panel.grid.minor = element_blank()
  )+ylab(ylab)+xlab(xlab)
  til=''
  if(test_method=='anova'){
    if(length(unique(data_m[,1]))<3){
      x1=data_m[,2][which(data_m[,1]==unique(data_m[,1])[1])]
      x2=data_m[,2][which(data_m[,1]==unique(data_m[,1])[2])]
      pv=t.test(x1,x2)$p.value
      til=paste0('t-tests p=',signif(pv,2))
    }else{
      fit <- aov(value~Group, data = data_m)
      pv=summary(fit)[[1]][5][[1]]
      fv=summary(fit)[[1]][4][[1]]
      til=paste0('ANOVA tests p=',signif(pv,2))
    }
  }else{
    if(length(unique(data_m[,1]))<3){
      x1=data_m[,2][which(data_m[,1]==unique(data_m[,1])[1])]
      x2=data_m[,2][which(data_m[,1]==unique(data_m[,1])[2])]
      pv=wilcox.test(x1,x2)$p.value
      til=paste0('wilcox.tests p=',signif(pv,2))
    }else{
      fit=kruskal.test(value~Group, data = data_m)
      pv=fit$p.value
      til=paste0('Kruskal-Wallis test p=',signif(pv,2))
    }
  }
  p1=p1+ggtitle(til)
  if(!is.null(ylim)){
    p1=p1+ylim(ylim)
  }
  if(is.null(show_compare)){
    if(length(uni.group)>5){
      show_compare=F
    }else{
      show_compare=T
    }
  }
  if(show_compare){
    comps=list()
    for(i in 1:(length(uni.group)-1)){
      for(j in (i+1):length(uni.group)){
        comps=c(comps,list(c(uni.group[i],uni.group[j])))
      }
    }
    if(length(comps)<5){
      p1=p1+ggpubr::stat_compare_means(comparisons = comps,method = cmp_test_method,label= "p.signif")
    }
  }
  return(p1)
}

plotMutiBar_use=function(dat,ist=F,margin=T,xlb='',ylb='',lineCol='black',lineW=0.5,legTitle='Group',showValue=F,showLine=T,xangle=0,isAuto=T){
  library(ggplot2)
  #library(tidyverse)
  #library(reshape2)
  #library(optparse)
  if(ist){
    dat=t(dat)
  }
  lbc=colnames(dat)
  lbr=row.names(dat)
  bk_dat=dat
  if(margin){
    dat=dat%*%diag(1/c(apply(t(dat), 1, sum)))
  }
  row.names(dat)=paste0('R',1:(nrow(dat)))
  colnames(dat)=paste0('C',1:(ncol(dat)))
  row.names(bk_dat)=paste0('R',1:(nrow(bk_dat)))
  colnames(bk_dat)=paste0('C',1:(ncol(bk_dat)))
  #df=cbind(bg=paste0('R',1:nrow(dat)),dat)
  #colnames(df)=c('bg',paste0('C',1:(ncol(dat))))
  tp.dat=as.data.frame(cbind(bg=row.names(dat),dat))
  tp.dat[,1]=as.character(tp.dat[,1])
  for(i in 2:ncol(tp.dat)){
    tp.dat[,i]=as.numeric(as.character(tp.dat[,i]))
  }
  mt.df=reshape2::melt(tp.dat)
  colnames(mt.df)=c('bg','variable','value')
  pg=ggplot(mt.df, aes(x=variable, y=value, fill=bg))+geom_bar(stat = "identity", width=lineW, col=lineCol)
  if(showLine){
    for (i in 2:(ncol(tp.dat)-1)) {
      tmp=tp.dat[order(tp.dat[,1],decreasing = T),]
      tmp[,i]=base::cumsum(tmp[,i])
      tmp[,i+1]=base::cumsum(tmp[,i+1])
      colnames(tmp)[c(i,i+1)]=c('STY','ED')
      tmp1=cbind(tmp,STX=rep(i-1+lineW/2,nrow(tmp))
                 ,EDX=rep(i-lineW/2,nrow(tmp)))
      pg=pg+geom_segment(data=tmp1,aes(x=STX, xend=EDX, y=STY, yend=ED))
    }
  }
  
  if(showValue){
    pg=pg+geom_text(data=mt.df,aes(label=sprintf("%0.2f", round(value, digits = 2))),position=position_stack(vjust=0.5))
  }
  pg=pg+scale_x_discrete(breaks = paste0('C',1:(ncol(dat))),label = lbc)
  pg=pg+labs(x=xlb, y=ylb)+ggsci::scale_fill_d3()+theme(legend.position = "bottom")
  pg=pg+scale_fill_manual(values = c(pal_locuszoom(alpha =0.8)(7)[c(5,1,2,3,6,7)]),breaks = paste0('R',1:nrow(dat)),label = lbr,name=legTitle)
  if(xangle>0){
    pg=pg+theme(axis.text.x = element_text(angle = xangle, hjust = 1),legend.position = "bottom")
  }
  
  g.tb=matrix(0,nrow=ncol(dat),ncol=ncol(dat))
  for(i in 1:(ncol(dat))){
    for(j in 1:ncol(dat)){
      if(i!=j){
        g.tb[i,j]=round(-log10((chisq.test(bk_dat[,c(i,j)])$p.value)),2)
      }
    }
  }
  colnames(g.tb)=lbc
  row.names(g.tb)=lbc
  g.tb=reshape2::melt(g.tb)
  colnames(g.tb)=c('A1','A2','A3')
  g.tb$A4=paste0(g.tb[,3],ifelse(g.tb[,3]>-log10(0.05),'(*)',''))
  stable.p=ggplot(g.tb, aes(A1, A2)) + geom_tile(aes(fill = A3),colour = "white") +xlab('')+ylab('')+ scale_fill_gradient(low = "white",high = "steelblue")+geom_text(aes(x=A1,y=A2,label=A4))+theme(legend.position="none",axis.title.x=element_blank(),axis.text.x=element_blank())
  stable.p=stable.p+ggtitle('-log10(anova p value)')
  if(isAuto){
    g1=ggpubr::ggarrange(stable.p,pg, ncol = 1, nrow = 2,heights = c(0.5,1),align = "hv")
    return(g1)
  }else{
    return(list(Bar=pg,Table=stable.p))
  }
}

plotMutiBar_tmp=function(dat,ist=F,margin=T,xlb='',ylb='',lineCol='black',lineW=0.5,legTitle='Group',showValue=F,showLine=T,xangle=0,isAuto=T,legend.nrow=0){
  library(ggplot2)
  #library(tidyverse)
  #library(reshape2)
  #library(optparse)
  if(ist){
    dat=t(dat)
  }
  lbc=colnames(dat)
  lbr=row.names(dat)
  bk_dat=dat
  if(margin){
    dat=dat%*%diag(1/c(apply(t(dat), 1, sum)))
  }
  row.names(dat)=paste0('R',1:(nrow(dat)))
  colnames(dat)=paste0('C',1:(ncol(dat)))
  row.names(bk_dat)=paste0('R',1:(nrow(bk_dat)))
  colnames(bk_dat)=paste0('C',1:(ncol(bk_dat)))
  #df=cbind(bg=paste0('R',1:nrow(dat)),dat)
  #colnames(df)=c('bg',paste0('C',1:(ncol(dat))))
  tp.dat=as.data.frame(cbind(bg=row.names(dat),dat))
  tp.dat[,1]=as.character(tp.dat[,1])
  for(i in 2:ncol(tp.dat)){
    tp.dat[,i]=as.numeric(as.character(tp.dat[,i]))
  }
  mt.df=reshape2::melt(tp.dat)
  colnames(mt.df)=c('bg','variable','value')
  
  pg=ggplot(mt.df, aes(x=variable, y=value, fill=bg))+geom_bar(stat = "identity", width=lineW, col=lineCol) 
  if(showLine){
    for (i in 2:(ncol(tp.dat)-1)) {
      tmp=tp.dat[order(tp.dat[,1],decreasing = T),]
      tmp[,i]=base::cumsum(tmp[,i])
      tmp[,i+1]=base::cumsum(tmp[,i+1])
      colnames(tmp)[c(i,i+1)]=c('STY','ED')
      tmp1=cbind(tmp,STX=rep(i-1+lineW/2,nrow(tmp))
                 ,EDX=rep(i-lineW/2,nrow(tmp)))
      pg=pg+geom_segment(data=tmp1,aes(x=STX, xend=EDX, y=STY, yend=ED))
    }
  }
  
  if(showValue){
    pg=pg+geom_text(data=mt.df,aes(label=sprintf("%0.2f", round(value, digits = 2))),position=position_stack(vjust=0.5))
  }
  pg=pg+scale_x_discrete(breaks = paste0('C',1:(ncol(dat))),label = lbc)
  pg=pg+labs(x=xlb, y=ylb)+ggsci::scale_fill_d3()+theme(legend.position = "bottom",legend.key.size = unit(0.2, 'cm'))
  # pg=pg+ggsci::scale_fill_d3()+scale_fill_discrete(breaks = paste0('R',1:nrow(dat)),label = lbr,name=legTitle)
  pg=pg+scale_fill_manual(values = c(pal_locuszoom(alpha =0.8)(7)[c(5,1,2,3,6,7)]),breaks = paste0('R',1:nrow(dat)),label = lbr,name=legTitle)
  if(xangle>0){
    pg=pg+theme(axis.text.x = element_text(angle = xangle, hjust = 1),legend.position = "bottom")
  }
  if(legend.nrow>0){
    print('YES')
    print(legend.nrow)
    pg=pg+guides(fill = guide_legend(nrow = legend.nrow, byrow = TRUE))
  }
  g.tb=matrix(0,nrow=ncol(dat),ncol=ncol(dat))
  for(i in 1:(ncol(dat))){
    for(j in 1:ncol(dat)){
      if(i!=j){
        g.tb[i,j]=round(-log10((chisq.test(bk_dat[,c(i,j)])$p.value)),2)
      }
    }
  }
  colnames(g.tb)=lbc
  row.names(g.tb)=lbc
  g.tb=reshape2::melt(g.tb) 
  colnames(g.tb)=c('A1','A2','A3')
  g.tb$A4=paste0(g.tb[,3],ifelse(g.tb[,3]>-log10(0.05),'(*)',''))
  stable.p=ggplot(g.tb, aes(A1, A2)) + geom_tile(aes(fill = A3),colour = "white") +xlab('')+ylab('')+ scale_fill_gradient(low = "white",high = "steelblue")+geom_text(aes(x=A1,y=A2,label=A4))+theme(legend.position="none",axis.title.x=element_blank(),axis.text.x=element_blank())
  stable.p=stable.p+ggtitle('-log10(anova p value)')
  if(isAuto){
    g1=ggpubr::ggarrange(stable.p,pg, ncol = 1, nrow = 2,heights = c(0.5,1),align = "hv")
    return(g1)
  }else{
    return(list(Bar=pg,Table=stable.p))
  }
}

###################################################
################### GEO 数据集
# ## GSE37745 dataset
# # GSE37745=getGEOExpDataByCel('GSE37745')
# # save(GSE37745,file='origin_datas/GEO/GSE37745.exp.rma.RData')
# load('origin_datas/GEO/GSE37745.exp.rma.RData')
# gse37745.cli<-GSE37745$Sample
# gse37745.t.cli<-gse37745.cli[gse37745.cli$histology=='adeno',]
# gse37745.t.cli.tnm<-clean_TNMStage(ss = gse37745.cli$`tumor stage`
#                                    , sex = gse37745.cli$gender
#                                    , age = gse37745.cli$age
#                                    , age_cut = NULL)
# rownames(gse37745.t.cli.tnm)<-gse37745.cli$Acc
# 
# gse37745.t.cli.os<-data.frame(OS=gse37745.t.cli$dead
#                               ,OS.time=gse37745.t.cli$`days to determined death status`
#                               ,PFI=gse37745.t.cli$recurrence
#                               ,PFI.time=gse37745.t.cli$`days to recurrence / to last visit`)
# rownames(gse37745.t.cli.os)<-gse37745.t.cli$Acc
# gse37745.t.cli.os<-crbind2DataFrame(gse37745.t.cli.os)
# table(gse37745.t.cli.os$OS)
# gse37745.t.cli.os$OS[gse37745.t.cli.os$OS=='yes']<-1
# gse37745.t.cli.os$OS[gse37745.t.cli.os$OS=='no']<-0
# 
# table(gse37745.t.cli.os$PFI)
# gse37745.t.cli.os$PFI[gse37745.t.cli.os$PFI=='yes']<-1
# gse37745.t.cli.os$PFI[gse37745.t.cli.os$PFI=='no']<-0
# gse37745.t.cli.os$PFI[gse37745.t.cli.os$PFI=='not known']<-NA
# 
# gse37745.t.cli.os<-crbind2DataFrame(gse37745.t.cli.os)
# gse37745.t.cli.os=gse37745.t.cli.os[gse37745.t.cli.os$OS.time>0,]
# 
# # gse37745.exp<-GSE37745$Exp$GPL570_Data
# # gse37745.exp<-exp_probe2symbol_v2(gse37745.exp,GPL='GPL570')
# # save(gse37745.exp,file='origin_datas/GEO/gse37745.exp.RData')
# load('origin_datas/GEO/gse37745.exp.RData')
# 
# gse37745.t.exp<-gse37745.exp[,rownames(gse37745.t.cli.os)]
# dim(gse37745.t.exp)
# gse37745.t.cli.tnm<-gse37745.t.cli.tnm[rownames(gse37745.t.cli.os),]
# 
# writeMatrix(gse37745.t.exp,outpath = 'files/文件/GSE37745.t.exp.txt')
# writeMatrix(cbind(gse37745.t.cli.tnm,gse37745.t.cli.os),outpath = 'files/文件/GSE37745.t.exp.clini.txt')
# 

## GSE50081 dataset
# GSE50081=getGEOExpDataByCel('GSE50081')
# save(GSE50081,file='origin_datas/GEO/GSE50081.exp.rma.RData')
load('origin_datas/GEO/GSE50081.exp.rma.RData')
gse50081.cli<-GSE50081$Sample
gse50081.t.cli<-gse50081.cli[gse50081.cli$histology=='adenocarcinoma',]
gse50081.t.cli.tnm<-clean_TNMStage(ss = gse50081.cli$Stage
                                   ,st=gse50081.cli$`t-stage`
                                   ,sn=gse50081.cli$`n-stage`
                                   ,sm=gse50081.cli$`m-stage`
                                   , sex = gse50081.cli$Sex
                                   , age = gse50081.cli$age
                                   , age_cut = NULL)
rownames(gse50081.t.cli.tnm)<-gse50081.cli$Acc

gse50081.t.cli.os<-data.frame(OS=gse50081.t.cli$status
                              ,OS.time=gse50081.t.cli$`survival time`
                              ,PFI=gse50081.t.cli$recurrence
                              ,PFI.time=gse50081.t.cli$`disease-free survival time`)
rownames(gse50081.t.cli.os)<-gse50081.t.cli$Acc
gse50081.t.cli.os<-crbind2DataFrame(gse50081.t.cli.os)
table(gse50081.t.cli.os$OS)
gse50081.t.cli.os$OS[gse50081.t.cli.os$OS=='dead']<-1
gse50081.t.cli.os$OS[gse50081.t.cli.os$OS=='alive']<-0

table(gse50081.t.cli.os$PFI)
gse50081.t.cli.os$PFI[gse50081.t.cli.os$PFI=='Y']<-1
gse50081.t.cli.os$PFI[gse50081.t.cli.os$PFI=='N']<-0
gse50081.t.cli.os$PFI[gse50081.t.cli.os$PFI=='U']<-NA

gse50081.t.cli.os<-crbind2DataFrame(gse50081.t.cli.os)
gse50081.t.cli.os$OS.time=gse50081.t.cli.os$OS.time*365
gse50081.t.cli.os=gse50081.t.cli.os[gse50081.t.cli.os$OS.time>0,]

# gse50081.exp<-GSE50081$Exp$GPL570_Data
# gse50081.exp<-exp_probe2symbol_v2(gse50081.exp,GPL='GPL570')
# save(gse50081.exp,file='origin_datas/GEO/gse50081.exp.RData')
load('origin_datas/GEO/gse50081.exp.RData')

table(gse50081.t.cli.os$OS.time>30)
table(gse50081.t.cli.os$OS.time<365*15)

gse50081.t.exp<-gse50081.exp[,rownames(gse50081.t.cli.os)]
dim(gse50081.t.exp)
gse50081.t.cli.tnm<-gse50081.t.cli.tnm[rownames(gse50081.t.cli.os),]

writeMatrix(gse50081.t.exp,outpath = 'files/文件/GSE50081.t.exp.txt')
writeMatrix(cbind(gse50081.t.cli.tnm,gse50081.t.cli.os),outpath = 'files/文件/GSE50081.t.exp.clini.txt')

# ## GSE19188 dataset
# # GSE19188=getGEOExpDataByCel('GSE19188')
# # save(GSE19188,file='origin_datas/GEO/GSE19188.exp.rma.RData')
# load('origin_datas/GEO/GSE19188.exp.rma.RData')
# gse19188.cli<-GSE19188$Sample
# gse19188.t.cli<-gse19188.cli[gse19188.cli$`cell type`=='ADC',]
# gse19188.t.cli.tnm<-clean_TNMStage(sex = gse19188.cli$gender)
# rownames(gse19188.t.cli.tnm)<-gse19188.cli$Acc
# 
# gse19188.t.cli.os<-data.frame(OS=gse19188.t.cli$status
#                               ,OS.time=gse19188.t.cli$`overall survival`
#                               ,PFI=NA
#                               ,PFI.time=NA
# )
# rownames(gse19188.t.cli.os)<-gse19188.t.cli$Acc
# gse19188.t.cli.os<-crbind2DataFrame(gse19188.t.cli.os)
# table(gse19188.t.cli.os$OS)
# gse19188.t.cli.os$OS[gse19188.t.cli.os$OS=='deceased']<-1
# gse19188.t.cli.os$OS[gse19188.t.cli.os$OS=='alive']<-0
# gse19188.t.cli.os$OS[gse19188.t.cli.os$OS=='Not available']<-NA
# 
# gse19188.t.cli.os<-crbind2DataFrame(gse19188.t.cli.os)
# gse19188.t.cli.os$OS.time[gse19188.t.cli.os$OS.time=='Not available']=NA
# gse19188.t.cli.os<-crbind2DataFrame(gse19188.t.cli.os)
# gse19188.t.cli.os$OS.time=gse19188.t.cli.os$OS.time*30
# gse19188.t.cli.os=gse19188.t.cli.os[!is.na(gse19188.t.cli.os$OS.time),]
# gse19188.t.cli.os<-gse19188.t.cli.os[gse19188.t.cli.os$OS.time>0,]
# 
# # gse19188.exp<-GSE19188$Exp$GPL570_Data
# # gse19188.exp<-exp_probe2symbol_v2(gse19188.exp,GPL='GPL570')
# # save(gse19188.exp,file='origin_datas/GEO/gse19188.exp.RData')
# load('origin_datas/GEO/gse19188.exp.RData')
# 
# gse19188.t.exp<-gse19188.exp[,rownames(gse19188.t.cli.os)]
# dim(gse19188.t.exp)
# gse19188.t.cli.tnm<-gse19188.t.cli.tnm[rownames(gse19188.t.cli.os),]
# 
# writeMatrix(gse19188.t.exp,outpath = 'files/文件/GSE19188.t.exp.txt')
# writeMatrix(cbind(gse19188.t.cli.tnm,gse19188.t.cli.os),outpath = 'files/文件/GSE19188.t.exp.clini.txt')
# 
# ## GSE30219 dataset
# # GSE30219=getGEOExpDataByCel('GSE30219')
# # save(GSE30219,file='origin_datas/GEO/GSE30219.exp.rma.RData')
# load('origin_datas/GEO/GSE30219.exp.rma.RData')
# gse30219.cli<-GSE30219$Sample
# gse30219.t.cli<-gse30219.cli[gse30219.cli$histology=='ADC',]
# gse30219.t.cli.tnm<-clean_TNMStage(st = gse30219.cli$`pt stage`
#                                    ,sn=gse30219.t.cli$`pn stage`
#                                    ,sm=gse30219.t.cli$`pm stage`
#                                    , sex = gse30219.cli$gender
#                                    , age = gse30219.cli$`age at surgery`
#                                    , age_cut = NULL)
# rownames(gse30219.t.cli.tnm)<-gse30219.cli$Acc
# 
# gse30219.t.cli.os<-data.frame(OS=gse30219.t.cli$status
#                               ,OS.time=gse30219.t.cli$`follow-up time (months)`
#                               ,PFI=gse30219.t.cli$`relapse (event=1; no event=0)`
#                               ,PFI.time=gse30219.t.cli$`disease free survival in months`)
# rownames(gse30219.t.cli.os)<-gse30219.t.cli$Acc
# gse30219.t.cli.os<-crbind2DataFrame(gse30219.t.cli.os)
# table(gse30219.t.cli.os$OS)
# gse30219.t.cli.os$OS[gse30219.t.cli.os$OS=='DEAD']<-1
# gse30219.t.cli.os$OS[gse30219.t.cli.os$OS=='ALIVE']<-0
# 
# table(gse30219.t.cli.os$PFI)
# # gse30219.t.cli.os$PFI[gse30219.t.cli.os$PFI=='relapsed']<-1
# # gse30219.t.cli.os$PFI[gse30219.t.cli.os$PFI=='not relapsed']<-0
# 
# gse30219.t.cli.os<-crbind2DataFrame(gse30219.t.cli.os)
# gse30219.t.cli.os$OS.time=gse30219.t.cli.os$OS.time*30
# gse30219.t.cli.os$PFI.time=gse30219.t.cli.os$PFI.time*30
# gse30219.t.cli.os=gse30219.t.cli.os[gse30219.t.cli.os$OS.time>0,]
# 
# # gse30219.exp<-GSE30219$Exp$GPL570_Data
# # gse30219.exp<-exp_probe2symbol_v2(gse30219.exp,GPL='GPL570')
# # save(gse30219.exp,file='origin_datas/GEO/gse30219.exp.RData')
# load('origin_datas/GEO/gse30219.exp.RData')
# 
# gse30219.t.exp<-gse30219.exp[,rownames(gse30219.t.cli.os)]
# dim(gse30219.t.exp)
# gse30219.t.cli.tnm<-gse30219.t.cli.tnm[rownames(gse30219.t.cli.os),]
# 
# writeMatrix(gse30219.t.exp,outpath = 'files/文件/GSE30219.t.exp.txt')
# writeMatrix(cbind(gse30219.t.cli.tnm,gse30219.t.cli.os),outpath = 'files/文件/GSE30219.t.exp.clini.txt')

## GSE31210 dataset
# GSE31210=getGEOExpData('GSE31210')
# save(GSE31210,file='origin_datas/GEO/GSE31210.RData')
load('origin_datas/GEO/GSE31210.RData')
gse31210.cli<-GSE31210$Sample
gse31210.t.cli<-gse31210.cli[gse31210.cli$tissue=='primary lung tumor',]
gse31210.t.cli.tnm<-clean_TNMStage(ss = gse31210.cli$`pathological stage`
                                   , sex = gse31210.cli$gender
                                   , age = gse31210.cli$`age (years)`
                                   , age_cut = NULL)
rownames(gse31210.t.cli.tnm)<-gse31210.cli$Acc

gse31210.t.cli.os<-data.frame(OS=gse31210.t.cli$death
                              ,OS.time=gse31210.t.cli$`days before death/censor`
                              ,PFI=gse31210.t.cli$relapse
                              ,PFI.time=gse31210.t.cli$`days before relapse/censor`)
rownames(gse31210.t.cli.os)<-gse31210.t.cli$Acc
gse31210.t.cli.os<-crbind2DataFrame(gse31210.t.cli.os)
table(gse31210.t.cli.os$OS)
gse31210.t.cli.os$OS[gse31210.t.cli.os$OS=='dead']<-1
gse31210.t.cli.os$OS[gse31210.t.cli.os$OS=='alive']<-0

table(gse31210.t.cli.os$PFI)
gse31210.t.cli.os$PFI[gse31210.t.cli.os$PFI=='relapsed']<-1
gse31210.t.cli.os$PFI[gse31210.t.cli.os$PFI=='not relapsed']<-0

gse31210.t.cli.os<-crbind2DataFrame(gse31210.t.cli.os)
table(gse31210.t.cli.os$OS.time>30)

# gse31210.exp<-GSE31210$Exp$GPL570_54675_Data_col1
# gse31210.exp<-exp_probe2symbol_v2(gse31210.exp,GPL='GPL570')
# save(gse31210.exp,file='origin_datas/GEO/gse31210.exp.RData')
load('origin_datas/GEO/gse31210.exp.RData')
gse31210.t.cli.os=gse31210.t.cli.os[which(gse31210.t.cli.os$OS.time>0),]

gse31210.t.exp<-gse31210.exp[,rownames(gse31210.t.cli.os)]
dim(gse31210.t.exp)
gse31210.t.cli.tnm<-gse31210.t.cli.tnm[rownames(gse31210.t.cli.os),]
dim(gse31210.t.exp)
range(gse31210.t.exp)

writeMatrix(gse31210.t.exp,outpath = 'files/文件/GSE31210.t.exp.txt')
writeMatrix(cbind(gse31210.t.cli.tnm,gse31210.t.cli.os),outpath = 'files/文件/GSE31210.t.exp.clini.txt')

## GSE72094 dataset
# GSE72094=getGEOExpData('GSE72094')
# save(GSE72094,file='origin_datas/GEO/GSE72094.RData')
load('origin_datas/GEO/GSE72094.RData')
gse72094.cli<-GSE72094$Sample
rownames(gse72094.cli)=gse72094.cli$Acc
gse72094.t.cli<-gse72094.cli
colnames(gse72094.cli)
gse72094.t.cli.tnm<-clean_TNMStage(ss = gse72094.t.cli$Stage
                                   , sex = gse72094.t.cli$gender
                                   , age = gse72094.t.cli$age_at_diagnosis
                                   , age_cut = NULL)
rownames(gse72094.t.cli.tnm)<-gse72094.t.cli$Acc

gse72094.t.cli.os<-data.frame(OS=gse72094.t.cli$vital_status
                              ,OS.time=gse72094.t.cli$survival_time_in_days)
rownames(gse72094.t.cli.os)<-gse72094.t.cli$Acc
gse72094.t.cli.os<-crbind2DataFrame(gse72094.t.cli.os)
table(gse72094.t.cli.os$OS)
gse72094.t.cli.os$OS[gse72094.t.cli.os$OS=='Dead']<-1
gse72094.t.cli.os$OS[gse72094.t.cli.os$OS=='Alive']<-0
gse72094.t.cli.os$OS[gse72094.t.cli.os$OS=='NA']<-NA

gse72094.t.cli.os<-crbind2DataFrame(gse72094.t.cli.os)
gse72094.t.cli.os<-gse72094.t.cli.os[!is.na(gse72094.t.cli.os$OS),]
gse72094.t.cli.os<-gse72094.t.cli.os[!is.na(gse72094.t.cli.os$OS.time),]

gse72094.t.cli.os<-gse72094.t.cli.os[gse72094.t.cli.os$OS.time>0,]
dim(gse72094.t.cli.os)

# gse72094.exp<-GSE72094$Exp$GPL15048_60607_Data_col1
# gse72094.exp<-exp_probe2symbol_v2(gse72094.exp,anno = GSE72094$Anno$GPL15048[,c(1,4)])
# save(gse72094.exp,file='origin_datas/GEO/gse72094.exp.RData')
load('origin_datas/GEO/gse72094.exp.RData')

gse72094.t.cli.tnm=gse72094.t.cli.tnm[rownames(gse72094.t.cli.os),]
gse72094.t.exp=gse72094.exp[,rownames(gse72094.t.cli.os)]
dim(gse72094.t.exp)
range(gse72094.t.exp)

writeMatrix(gse72094.t.exp,outpath = 'files/文件/GSE72094.t.exp.txt')
writeMatrix(cbind(gse72094.t.cli.tnm,gse72094.t.cli.os),outpath = 'files/文件/GSE72094.t.exp.clini.txt')

############
ann.pcg=ann[which(ann$gene_type=='protein_coding'),]
dim(ann.pcg)

########### GDC 下载数据处理
luad.tcga.cli=readMatrix('origin_datas/TCGA/Merge_TCGA-LUAD_clinical.txt',row = F)
dim(luad.tcga.cli)
luad.tcga.cli$Tumor_Sample_Barcode=luad.tcga.cli$A0_Samples
luad.tcga.cli$os=luad.tcga.cli$A1_OS
luad.tcga.cli$os_status=luad.tcga.cli$A2_Event
dim(luad.tcga.cli)
luad.tcga.cli=luad.tcga.cli[,c(179:181,1:178)]
names(table(luad.tcga.cli$os_status))

luad.tcga.cli=data.frame(SampleID=paste0(luad.tcga.cli$Tumor_Sample_Barcode,"-01"),luad.tcga.cli,stringsAsFactors = F,check.names = F)
rownames(luad.tcga.cli)=luad.tcga.cli$SampleID

table(luad.tcga.cli$os_status)
luad.tcga.cli$os_status[luad.tcga.cli$os_status=='Alive']=0
luad.tcga.cli$os_status[luad.tcga.cli$os_status=='Dead']=1
luad.tcga.cli$os_status[luad.tcga.cli$os_status=='']=NA
table(luad.tcga.cli$os_status)
luad.tcga.cli=luad.tcga.cli[which(!is.na(luad.tcga.cli$os_status)),]
###TCGA-LUAD FPKM数据
# luad.tcga.exp=read.table("origin_datas/TCGA/Merge_TCGA-LUAD_FPKM.txt"
#                          ,row.names = 1,sep="\t",header=T,as.is=T,quote="\""
#                          ,fill=T,check.names = F,stringsAsFactors = F)
# rownames(luad.tcga.exp)=gsub("\\..*","",rownames(luad.tcga.exp))
# luad.tcga.exp=mg_FPKM2TPMs(luad.tcga.exp)
# luad.tcga.exp=log2(luad.tcga.exp+1)
# luad.tcga.exp=exp_ensg2symbol(luad.tcga.exp,method = 'median')
# save(luad.tcga.exp,file = 'origin_datas/TCGA/luad.tcga.exp.RData')
load('origin_datas/TCGA/luad.tcga.exp.RData')

dim(luad.tcga.exp)
range(luad.tcga.exp)
table(substr(colnames(luad.tcga.exp),14,15))

######### -06 转移肿瘤样本
luad.tcga.exp.t.inds=grep('-01[A-Z]?',colnames(luad.tcga.exp))
length(luad.tcga.exp.t.inds)

luad.tcga.t.exp=luad.tcga.exp[,luad.tcga.exp.t.inds]
dim(luad.tcga.t.exp)

##############################
luad.tcga.t.cli=luad.tcga.cli[grep("-01$",luad.tcga.cli$SampleID),]
dim(luad.tcga.t.cli)

luad.tcga.t.cli=dplyr::distinct(luad.tcga.t.cli,SampleID,.keep_all=TRUE)
rownames(luad.tcga.t.cli)=luad.tcga.t.cli$SampleID
dim(luad.tcga.t.cli)

##########
setdiff(colnames(luad.tcga.t.exp),luad.tcga.t.cli$SampleID)

luad.tcga.t.cli=luad.tcga.t.cli[colnames(luad.tcga.t.exp),]
dim(luad.tcga.t.cli)

luad.tcga.t.cli=dplyr::rename(luad.tcga.t.cli,c('OS.time'='A1_OS'
                                                ,'OS'='A2_Event'))
table(luad.tcga.t.cli$OS)
luad.tcga.t.cli$OS[luad.tcga.t.cli$OS=='Alive']=0
luad.tcga.t.cli$OS[luad.tcga.t.cli$OS=='Dead']=1
luad.tcga.t.cli$OS[luad.tcga.t.cli$OS=='']=NA
table(luad.tcga.t.cli$OS)

luad.tcga.t.cli=crbind2DataFrame(luad.tcga.t.cli)

table(luad.tcga.t.cli$OS)
table(is.na(luad.tcga.t.cli$OS))
table(is.na(luad.tcga.t.cli$OS.time))

luad.tcga.t.cli<-luad.tcga.t.cli[!is.na(luad.tcga.t.cli$OS.time),]
luad.tcga.t.cli<-luad.tcga.t.cli[!is.na(luad.tcga.t.cli$OS),]
dim(luad.tcga.t.cli)

range(luad.tcga.t.cli$OS.time)

table(luad.tcga.t.cli$OS.time<30)
table(luad.tcga.t.cli$OS.time>365*10)
table(luad.tcga.t.cli$OS.time>365*15)

###
luad.tcga.t.cli=luad.tcga.t.cli[intersect(luad.tcga.t.cli$SampleID,colnames(luad.tcga.t.exp)), ]
dim(luad.tcga.t.cli)
luad.tcga.t.exp=luad.tcga.t.exp[,luad.tcga.t.cli$SampleID]
dim(luad.tcga.t.exp)

## TCGA 生存数据
luad.tcga.t.exp.os=luad.tcga.t.cli[,c('OS','OS.time')]
luad.tcga.t.exp.os=luad.tcga.t.exp.os[which(luad.tcga.t.exp.os$OS.time>0),]

luad.tcga.t.cli=luad.tcga.t.cli[rownames(luad.tcga.t.exp.os),]
luad.tcga.t.exp=luad.tcga.t.exp[,rownames(luad.tcga.t.exp.os)]
dim(luad.tcga.t.cli)
dim(luad.tcga.t.exp)

writeMatrix(luad.tcga.t.cli,outpath = 'files/文件/luad.tcga.t.cli.txt')
writeMatrix(luad.tcga.t.exp,outpath = 'files/文件/luad.tcga.t.exp.txt')

## TCGA 临床数据
clin.selected=c('A17_Age','A18_Sex','A3_T','A4_N','A5_M','A6_Stage')
setdiff(clin.selected,colnames(luad.tcga.t.cli))

luad.tcga.t.exp.cli=luad.tcga.t.cli[,clin.selected]

luad.tcga.t.exp.cli.tnm=clean_TNMStage(st=luad.tcga.t.exp.cli$A3_T
                                       , sn = luad.tcga.t.exp.cli$A4_N
                                       , sm = luad.tcga.t.exp.cli$A5_M
                                       , ss = luad.tcga.t.exp.cli$A6_Stage
                                       , sex = luad.tcga.t.exp.cli$A18_Sex
                                       , age = luad.tcga.t.exp.cli$A17_Age
                                       , age_cut = NULL)
row.names(luad.tcga.t.exp.cli.tnm)=rownames(luad.tcga.t.cli)
dim(luad.tcga.t.exp.cli.tnm)

colnames(luad.tcga.t.exp.cli.tnm)

table(luad.tcga.t.exp.cli.tnm$Clinical_T)
table(luad.tcga.t.exp.cli.tnm$Clinical_N)
table(luad.tcga.t.exp.cli.tnm$Clinical_M)
table(luad.tcga.t.exp.cli.tnm$Clinical_Stage)

#################### Step 1 ###########################
APOBEC.genes<-read.table('origin_datas/APOBEC.gene.family.txt'
                         ,header = T,check.names = F,stringsAsFactors = F)
head(APOBEC.genes)
APOBEC.genes$Gene_symbol[which(APOBEC.genes$Gene_symbol=='AID')]='AICDA'

########### 正常组织中的基因表达分析
GTEx.exprs=c()
for(gene in APOBEC.genes$Gene_symbol){
  MGID=mg_mysql_getGene_Info(gene)$GeneID
  df_expr=get_pan_tissue_expression(gid = MGID)
  df_expr=data.frame(df_expr,Gene=gene)
  GTEx.exprs=rbind(GTEx.exprs,df_expr)
}
table(GTEx.exprs$primarySite,GTEx.exprs$Gene)
GTEx.exprs=GTEx.exprs[which(GTEx.exprs$primarySite!='<not provided>'),]
head(GTEx.exprs)
GTEx.exprs_a=reshape2::dcast(GTEx.exprs,primarySite ~ Gene,value.var='val',mean)
dim(GTEx.exprs_a)  
GTEx.exprs_a[1:4,1:5]
rownames(GTEx.exprs_a)=GTEx.exprs_a$primarySite
GTEx.exprs_a=GTEx.exprs_a[,-1]
GTEx.exprs_a=as.matrix(GTEx.exprs_a)
GTEx.exprs_a=log2(GTEx.exprs_a+1)
range(GTEx.exprs_a)
GTEx.exprs_a=t(GTEx.exprs_a)

min(GTEx.exprs_a)
max(GTEx.exprs_a)

library(ComplexHeatmap)
fig1a=Heatmap(as.matrix(GTEx.exprs_a)
        # , name = "NES"
        , rect_gp = gpar(col = 'grey')
        , col = circlize::colorRamp2(c(0, 8), c('white', 'red'))
        , border = TRUE
        # , show_column_names = T
        , show_column_dend = F
        , show_row_dend = F
        # , cluster_columns=T
        # , cluster_rows=F
        # , rect_gp = gpar(col = "white", lwd = 1)
        # , row_names_gp = gpar(fontsize = 10)
        ,heatmap_legend_param = list(title = 'TPM')
        ,jitter=TRUE
)

pdf('PDFs/Fig1A.pdf',height = 4,width = 8)
fig1a
dev.off()

########### TCGA数据库中 -01肿瘤组织中的基因表达分析
pan.cancer.exprs=c()
for(gene in APOBEC.genes$Gene_symbol){
  # MGID=mg_mysql_getGene_Info(gene)$GeneID
  df_expr=mg_get_pancancer_deg_exp(gene=gene,onlyTCGA = F)
  # df_expr=get_pan_cancer_expression(gid = MGID)
  df_expr=data.frame(df_expr,Gene=gene)
  pan.cancer.exprs=rbind(pan.cancer.exprs,df_expr)
}

table(gsub("(.*?)-.*","\\1",pan.cancer.exprs$SampleName))
length(unique(pan.cancer.exprs$CODE))

pan.cancer.exprs=pan.cancer.exprs[gsub("(.*?)-.*","\\1",pan.cancer.exprs$SampleName) %in% c('TCGA','TARGET'),]
table(pan.cancer.exprs$CODE,pan.cancer.exprs$Type)

pan.cancer.exprs=pan.cancer.exprs[which(pan.cancer.exprs$Type=='Tumor'),]

pan.cancer.exprs_m=reshape2::dcast(pan.cancer.exprs,Gene ~ CODE,value.var='Value',mean)
dim(pan.cancer.exprs_m)  
pan.cancer.exprs_m[1:4,1:5]
rownames(pan.cancer.exprs_m)=pan.cancer.exprs_m$Gene
pan.cancer.exprs_m=pan.cancer.exprs_m[,-1]
range(pan.cancer.exprs_m)
pan.cancer.exprs_m=as.matrix(pan.cancer.exprs_m)
range(pan.cancer.exprs_m)
min(pan.cancer.exprs_m)
max(pan.cancer.exprs_m)
pan.cancer.exprs_m=log2(pan.cancer.exprs_m+1)
dim(pan.cancer.exprs_m)

library(ComplexHeatmap)
fig1b=Heatmap(as.matrix(pan.cancer.exprs_m)
              # , name = "NES"
              ,rect_gp = gpar(col = 'grey')
              , col = circlize::colorRamp2(c(0, 8), c('white', 'red'))
              , border = TRUE
              # , show_column_names = T
              , show_column_dend = F
              , show_row_dend = F
              # , cluster_columns=T
              # , cluster_rows=F
              # , rect_gp = gpar(col = "white", lwd = 1)
              # , row_names_gp = gpar(fontsize = 10)
              ,heatmap_legend_param = list(title = 'TPM')
              ,jitter=TRUE
)

pdf('PDFs/Fig1B.pdf',height = 4,width = 8)
fig1b
dev.off()

############# 肿瘤-正常组织中的差异分析[TCGA, TARGET, GTEX]
pancancer.deg.exp=c()
for(gene in APOBEC.genes$Gene_symbol){
  df_expr=mg_get_pancancer_deg_exp(gene=gene,onlyTCGA = F)
  df_expr=data.frame(df_expr,Gene=gene)
  pancancer.deg.exp=rbind(pancancer.deg.exp,df_expr)
}

head(pancancer.deg.exp)
table(substr(pancancer.deg.exp$SampleName,1,4))

pancancer.deg.exp_m=reshape2::dcast(pancancer.deg.exp,Gene ~ SampleName,value.var='Value',mean)
rownames(pancancer.deg.exp_m)=pancancer.deg.exp_m$Gene
pancancer.deg.exp_m=pancancer.deg.exp_m[,-1]
pancancer.deg.exp_m[1:3,1:5]
range(pancancer.deg.exp_m)
pancancer.deg.exp_m=log2(pancancer.deg.exp_m+1)
dim(pancancer.deg.exp_m)
range(pancancer.deg.exp_m)

pancancer.deg.res=c()
for (disease in unique(pancancer.deg.exp$CODE)) {
  # disease='UCEC'
  smp=pancancer.deg.exp$SampleName[pancancer.deg.exp$CODE==disease]
  df_exp=pancancer.deg.exp_m[,smp]
  df_group=pancancer.deg.exp$Type[pancancer.deg.exp$CODE==disease]
  limma.res=mg_limma_DEG(exp=df_exp,group = df_group,ulab = 'Tumor',dlab = 'Normal')
  limma.res=limma.res$DEG
  limma.res=cbind(Gene=rownames(limma.res),limma.res,CODE=disease)
  pancancer.deg.res=rbind(pancancer.deg.res,limma.res)
}

pancancer.deg.res_m=reshape2::dcast(pancancer.deg.res,Gene ~ CODE,value.var='logFC',mean)
rownames(pancancer.deg.res_m)=pancancer.deg.res_m$Gene
pancancer.deg.res_m=pancancer.deg.res_m[,-1]

pancancer.deg.res_pval=reshape2::dcast(pancancer.deg.res,Gene ~ CODE,value.var='adj.P.Val',mean)
rownames(pancancer.deg.res_pval)=pancancer.deg.res_pval$Gene
pancancer.deg.res_pval=pancancer.deg.res_pval[,-1]
pancancer.deg.res_pval=apply(pancancer.deg.res_pval, c(1,2), function(x){mg_format_p_values(x)})
pancancer.deg.res_pval=pancancer.deg.res_pval[rownames(pancancer.deg.res_m),colnames(pancancer.deg.res_m)]

col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "tomato"))
fig1c=Heatmap(as.matrix(pancancer.deg.res_m), name = "logFC"
        , col = col_fun
        ,cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(pancancer.deg.res_pval[i, j], x, y, gp = gpar(fontsize = 8))}
        ,show_row_dend = F
        ,show_column_dend = F
        ,rect_gp = gpar(col = 'grey')
        ,border = TRUE
        )
pdf('PDFs/Fig1C.pdf',height = 4,width = 8)
fig1c
dev.off()

fig1d=mg_pan_cancer_deg_boxplot(gene = 'APOBEC3B',group_col=c("#4DBBD5", "#E64B35")
                                ,ylab = 'APOBEC3B expression'
                                ,legend.pos='tr')
savePDF('PDFs/Fig1D.pdf',fig1d,height = 4,width = 8)

########### 配对样本中的表达
get_tcga_exp_pair=function(gene,CODE){
  # gene='APOBEC3B'
  # CODE='LUAD'
  df_expr=mg_get_pancancer_deg_exp(gene=gene,onlyTCGA = T)
  df_expr=df_expr[which(df_expr$CODE==CODE),]
  df_expr$Patient=substr(df_expr$SampleName,1,12)
  pair.smp=names(which(table(df_expr$Patient)>1))
  df_expr=df_expr[df_expr$Patient %in% pair.smp,]
  df_expr=data.frame(df_expr,Gene=gene,check.names = F,stringsAsFactors = F)
  df_expr=df_expr[order(df_expr$Type,df_expr$Patient),]
  # all(df_expr$Patient[df_expr$Type=='Normal']==(df_expr$Patient[df_expr$Type=='Tumor']))
  return(df_expr)
}

luad.tcga.expr.pair=c()
for(gene in APOBEC.genes$Gene_symbol){
  df_expr=get_tcga_exp_pair(gene = gene,CODE = 'LUAD')
  luad.tcga.expr.pair=rbind(luad.tcga.expr.pair,df_expr)
}
range(luad.tcga.expr.pair$Value)
luad.tcga.expr.pair$Value=log2(luad.tcga.expr.pair$Value+1)
all(luad.tcga.expr.pair$Patient[which(luad.tcga.expr.pair$Type=='Tumor')]==luad.tcga.expr.pair$Patient[which(luad.tcga.expr.pair$Type=='Normal')])

library("ggpubr")
p.all=list()
for(gene in unique(luad.tcga.expr.pair$Gene)){
  luad.tcga.expr.pair_sub=luad.tcga.expr.pair[which(luad.tcga.expr.pair$Gene==gene),]
  luad.tcga.expr.pair_m=reshape2::dcast(luad.tcga.expr.pair_sub,Patient~Type,value.var='Value',mean)
  p=ggpaired(luad.tcga.expr.pair_m[,c(2,3)], cond1 = 'Normal', cond2 = 'Tumor', fill = "condition", palette = c("#4DBBD5", "#E64B35"),
             legend.title="Group",xlab="",ylab = paste(gene,"expression",sep=" "))+
    stat_compare_means(paired = TRUE, label = "p.format", label.x = 1.35)
    # stat_compare_means(paired = TRUE, symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif",label.x = 1.35)
  p.all=c(p.all,list(p))
}

fig1e=mg_merge_plot(p.all,nrow = 2,ncol = 6,labels = c('E'),common.legend = T)
savePDF('PDFs/Fig1E.pdf',fig1e,width = 16,height = 8)

############################# 02_SNV_CNV ######
#################################
#### pan-cancer analysis
#################################
# tcga.pancancer.GeneMut=read.table('origin_datas/PANCAN/mc3.v0.2.8.PUBLIC.nonsilentGene.xena/mc3.v0.2.8.PUBLIC.nonsilentGene.xena',sep='\t',header = T,fill = T,check.names = F,stringsAsFactors = F)
# tcga.pancancer.GeneMut=tcga.pancancer.GeneMut[which(!is.na(tcga.pancancer.GeneMut$sample)),]
# rownames(tcga.pancancer.GeneMut)=tcga.pancancer.GeneMut$sample
# tcga.pancancer.GeneMut=tcga.pancancer.GeneMut[,-1]
# save(tcga.pancancer.GeneMut,file = 'origin_datas/PANCAN/tcga.pancancer.GeneMut.RData')
load('origin_datas/PANCAN/tcga.pancancer.GeneMut.RData')
dim(tcga.pancancer.GeneMut)
tcga.pancancer.GeneMut[1:4,1:5]
table(substr(colnames(tcga.pancancer.GeneMut),14,15))

#################
tcga.pancancer.GeneMut.sub=tcga.pancancer.GeneMut[match(APOBEC.genes$Gene_symbol,rownames(tcga.pancancer.GeneMut)),]
tcga.pancancer.GeneMut.sub=tcga.pancancer.GeneMut.sub[,which(substr(colnames(tcga.pancancer.GeneMut.sub),14,15)=="01")]
dim(tcga.pancancer.GeneMut.sub)
colnames(tcga.pancancer.GeneMut.sub)

##################
#####################
tcga.cancer.type$英文名称=str_to_title(tcga.cancer.type$英文名称)
############
tcga.pancancer.phenotype=read.table('origin_datas/PANCAN/TCGA_phenotype_denseDataOnlyDownload.tsv/TCGA_phenotype_denseDataOnlyDownload.tsv',sep='\t',header = T,fill = T,check.names = F,stringsAsFactors = F)
tcga.pancancer.phenotype[1:4,1:4]
rownames(tcga.pancancer.phenotype)=tcga.pancancer.phenotype$sample

setdiff(colnames(tcga.pancancer.GeneMut.sub),tcga.pancancer.phenotype$sample)
tcga.pancancer.phenotype=tcga.pancancer.phenotype[colnames(tcga.pancancer.GeneMut.sub),]
dim(tcga.pancancer.phenotype)

dim(tcga.pancancer.phenotype)
head(tcga.pancancer.phenotype)

tcga.pancancer.phenotype$`_primary_disease`=str_to_title(tcga.pancancer.phenotype$`_primary_disease`)
tcga.pancancer.phenotype$`_primary_disease`=gsub("Cancer$","Carcinoma",tcga.pancancer.phenotype$`_primary_disease`)
tcga.pancancer.phenotype$`_primary_disease`[which(tcga.pancancer.phenotype$`_primary_disease`=='Cervical & Endocervical Carcinoma')]='Cervical squamous cell carcinoma and endocervical adenocarcinoma'
tcga.pancancer.phenotype$`_primary_disease`[which(tcga.pancancer.phenotype$`_primary_disease`=='Diffuse Large B-Cell Lymphoma')]='Lymphoid Neoplasm Diffuse Large B-cell Lymphoma'
tcga.pancancer.phenotype$`_primary_disease`[which(tcga.pancancer.phenotype$`_primary_disease`=='Head & Neck Squamous Cell Carcinoma')]='Head and Neck squamous cell carcinoma'
tcga.pancancer.phenotype$`_primary_disease`[which(tcga.pancancer.phenotype$`_primary_disease`=='Kidney Clear Cell Carcinoma')]='Kidney renal clear cell carcinoma'
tcga.pancancer.phenotype$`_primary_disease`[which(tcga.pancancer.phenotype$`_primary_disease`=='Kidney Papillary Cell Carcinoma')]='Kidney renal papillary cell carcinoma'
tcga.pancancer.phenotype$`_primary_disease`[which(tcga.pancancer.phenotype$`_primary_disease`=='Pheochromocytoma & Paraganglioma')]='Pheochromocytoma and Paraganglioma'
tcga.pancancer.phenotype$`_primary_disease`[which(tcga.pancancer.phenotype$`_primary_disease`=='Testicular Germ Cell Tumor')]='Testicular Germ Cell Tumors'
tcga.pancancer.phenotype$`_primary_disease`[which(tcga.pancancer.phenotype$`_primary_disease`=='Uterine Corpus Endometrioid Carcinoma')]='Uterine Corpus Endometrial Carcinoma'
tcga.pancancer.phenotype$`_primary_disease`=str_to_title(tcga.pancancer.phenotype$`_primary_disease`)
setdiff(tcga.pancancer.phenotype$`_primary_disease`,tcga.cancer.type$英文名称)

head(tcga.pancancer.phenotype)
tcga.pancancer.phenotype$CODE=tcga.cancer.type$Cohort[match(tcga.pancancer.phenotype$`_primary_disease`,tcga.cancer.type$英文名称)]
head(tcga.pancancer.phenotype)

#########################
table(tcga.pancancer.phenotype$CODE,tcga.pancancer.phenotype$sample_type)

#####################
tcga.pancancer.stat=table(tcga.pancancer.phenotype$CODE,tcga.pancancer.phenotype$sample_type)
head(tcga.pancancer.stat)
tcga.pancancer.stat

################
tcga.pancancer.GeneMut.sub[1:4,1:5]
tcga.pancancer.GeneMut.sub_use=data.frame(Gene=rownames(tcga.pancancer.GeneMut.sub),tcga.pancancer.GeneMut.sub,check.names = F)
tcga.pancancer.GeneMut.sub_use[1:4,1:5]

tcga.pancancer.GeneMut.sub_long=reshape2::melt(tcga.pancancer.GeneMut.sub_use,id.vars='Gene')
tcga.pancancer.GeneMut.sub_long=merge(tcga.pancancer.GeneMut.sub_long,tcga.pancancer.phenotype,by.x='variable',by.y='sample',all.x=TRUE,sort=F)
head(tcga.pancancer.GeneMut.sub_long)

tcga.pancancer.GeneMut.sub_m=reshape2::dcast(tcga.pancancer.GeneMut.sub_long,CODE~Gene,value.var='value',sum)
dim(tcga.pancancer.GeneMut.sub_m)
rownames(tcga.pancancer.GeneMut.sub_m)=tcga.pancancer.GeneMut.sub_m$CODE
tcga.pancancer.GeneMut.sub_m=tcga.pancancer.GeneMut.sub_m[,-1]

tcga.pancancer.MF=tcga.pancancer.GeneMut.sub_m/tcga.pancancer.stat[rownames(tcga.pancancer.stat),1]
tcga.pancancer.MF=as.matrix(tcga.pancancer.MF)
rownames(tcga.pancancer.MF)=gsub("TCGA-","",rownames(tcga.pancancer.MF))
range(tcga.pancancer.MF)
tcga.pancancer.MF=tcga.pancancer.MF*100
tcga.pancancer.MF=t(tcga.pancancer.MF)

#######################################
library("ComplexHeatmap")
library("circlize")

heat_colors2 <- colorRamp2(c(0, 2.5, 5), c( "white","red", "red4"))

pdf('PDFs/Fig2A.pdf',width =8,height = 4)
Heatmap(tcga.pancancer.MF,
        name="MF (%)",
        border='grey',
        rect_gp = gpar(col = 'grey'),
        cluster_rows = F,
        col=heat_colors2,
        color_space = "RGB",
        cluster_columns = T,
        show_column_dend=F,
        show_row_dend=F,
        row_order=NULL,
        column_order=NULL,
        show_column_names = T,
        show_row_names = T,
        row_names_gp = gpar(fontsize = 12),
        gap = unit(2, "mm"),
        column_title = "",
        column_title_gp = gpar(fontsize = 6),
        show_heatmap_legend = TRUE,
        # cell_fun = function(j, i, x, y, width, height, fill) {
        #   grid.text(signif(tcga.pancancer.MF,digits = 1)[i, j], x, y, gp = gpar(fontsize = 8))},
        heatmap_legend_param=list(labels_gp = gpar(fontsize = 6),title_gp = gpar(fontsize = 6, fontface = "bold")))
dev.off()

######################### 肺腺癌
library(maftools)
LUAD=getTCGAMAFByCode('LUAD')
LUAD_new<-LUAD@data
######### 核验过，4个基因没有
setdiff(APOBEC.genes$Gene_symbol,unique(LUAD_new$Hugo_Symbol))
grep("KA6",APOBEC.genes$Gene_symbol)

# LUAD_new$Hugo_Symbol[which(LUAD_new$Hugo_Symbol=='DFNA5')]='GSDME'

##################
LUAD=read.maf(maf=LUAD_new,isTCGA = T)
LUAD_new=LUAD@data

####### 突变的maf
LUAD_new_APOBEC<-LUAD_new[which(LUAD_new$Hugo_Symbol %in% APOBEC.genes$Gene_symbol),]

unique(LUAD_new_APOBEC$Tumor_Sample_Barcode)

luad.tcga.cli.forMut=luad.tcga.cli[,-1]
luad.tcga.cli.forMut$os=luad.tcga.cli.forMut$os/365

writeMatrix(luad.tcga.cli.forMut,'origin_datas/TCGA/luad.tcga.cli.txt',row = F)

######## APOBEC.genes对应的maf文件
LUAD_new_APOBEC_maf<-read.maf(maf=LUAD_new_APOBEC
                              ,clinicalData ='origin_datas/TCGA/luad.tcga.cli.txt')

table(substr(unique(LUAD_new_APOBEC_maf@data$Tumor_Sample_Barcode),14,15))

######## 整个maf文件
LUAD_new_maf<-read.maf(maf=LUAD_new
                       ,clinicalData ='origin_datas/TCGA/luad.tcga.cli.txt'
                       ,isTCGA=T)

pdf('PDFs/Fig2B.pdf',width =8 ,height = 4)
oncoplot(maf = LUAD,genes = APOBEC.genes$Gene_symbol)
dev.off()

##############
luad.titv = titv(maf = LUAD_new_APOBEC_maf, plot = FALSE, useSyn = TRUE)
# plot titv summary
pdf('PDFs/Fig2A-2.pdf',width =8 ,height = 6)
plotTiTv(res = luad.titv)
dev.off()

############ 生存分析
mut_gene<-LUAD_new_APOBEC$Hugo_Symbol[which(!duplicated(LUAD_new_APOBEC$Hugo_Symbol))]
mut_gene

if(!dir.exists('files/01_SNV_CNV/mut_gene_survival')){
  dir.create('files/01_SNV_CNV/mut_gene_survival',recursive =T)
}

####### Mut/nor 单个基因生存分析
## ADARB2 os和pfs均显著
mut_gene
for (i in 1:length(mut_gene))
{
  pdf(paste0('files/01_SNV_CNV/mut_gene_survival/',mut_gene[i],'.pdf'),width = 4,height = )
  mafSurvival(maf=LUAD_new_maf, genes=mut_gene[i], time="os", Status="os_status", isTCGA=TRUE)
  dev.off()
}

### 整个 APOBEC 的生存分析
#### all survival
pdf(paste0('files/01_SNV_CNV/APOBEC.all.pdf'),width = 5,height = 5)
mafSurvival(maf=LUAD_new_maf, genes=mut_gene, time="os", Status="os_status", isTCGA=TRUE,col=mg_colors[1:2])
dev.off()

pdf(paste0('PDFs/Fig2C.pdf'),width = 4,height = 4)
mafSurvival(maf=LUAD_new_maf, genes=mut_gene, time="os", Status="os_status", isTCGA=TRUE,col=mg_colors[1:2])
dev.off()

### 单个基因的结构信息
if(!dir.exists('files/01_SNV_CNV/mut_gene_loolipplot/')){
  dir.create('files/01_SNV_CNV/mut_gene_loolipplot/',recursive = T)
}

for (i in 1:length(mut_gene))
{
  if(mut_gene[i]=='GSDME'){
    next;
  }
  pdf(paste0('files/01_SNV_CNV/mut_gene_loolipplot/',mut_gene[i],'.pdf'),width = 8,height = 4)
  lollipopPlot(maf=LUAD_new_maf,gene=mut_gene[i])
  dev.off()
}

### 具有 APOBEC 突变信息的样本
write.mafSummary(maf=LUAD_new_APOBEC_maf, basename="files/01_SNV_CNV/LUAD_APOBEC") 
## 具有突变的样本
write.mafSummary(maf=LUAD_new_maf, basename="files/01_SNV_CNV/LUAD") 

#####################################################
########### -01 肿瘤样本中 突变和非突变的 GSEA
##### 具有 APOBEC 基因突变的样本
mut_sample<-read.table('files/01_SNV_CNV/LUAD_APOBEC_sampleSummary.txt',header = T,check.names = F,row.names = 1)
mut_sample<-rownames(mut_sample)
length(mut_sample)

dim(luad.tcga.t.exp)
range(luad.tcga.t.exp)

dim(luad.tcga.exp)
range(luad.tcga.exp)

exp_matrix<-luad.tcga.exp
dim(exp_matrix)
table(substr(colnames(exp_matrix),14,15))

##### -01 肿瘤样本的转录表达谱
exp_matrix_T<-exp_matrix[,which(substr(colnames(exp_matrix),14,15) == "01")]
colnames(exp_matrix_T)<-substr(colnames(exp_matrix_T),1,12)
dim(exp_matrix_T)

# ##### -06 转移样本的转录表达谱
# exp_matrix_M<-exp_matrix[,which(substr(colnames(exp_matrix),14,15) == "06")]
# colnames(exp_matrix_M)<-substr(colnames(exp_matrix_M),1,12)
# dim(exp_matrix_M)

# ##### -11 转移样本的转录表达谱
# exp_matrix_N<-exp_matrix[,which(substr(colnames(exp_matrix),14,15) == "11")]
# colnames(exp_matrix_N)
# colnames(exp_matrix_N)<-substr(colnames(exp_matrix_N),1,12)
# dim(exp_matrix_N)

table(substr(colnames(exp_matrix),14,15)) 

## -01 肿瘤样本中具有 APOBEC 突变和非突变的表达谱
dim(exp_matrix_T)
length(mut_sample)

exp_matrix_M_mut<-exp_matrix_T[,which(colnames(exp_matrix_T) %in% mut_sample)]
exp_matrix_M_nor<-exp_matrix_T[,which(!(colnames(exp_matrix_T) %in% mut_sample))]
dim(exp_matrix_M_mut)
dim(exp_matrix_M_nor)

exp_matrix_for_deg_mutnor<-cbind(exp_matrix_M_mut,exp_matrix_M_nor)
dim(exp_matrix_for_deg_mutnor)
dim(exp_matrix_T)

if(!dir.exists('files/01_SNV_CNV/GSEA')){
  dir.create('files/01_SNV_CNV/GSEA',recursive = T)
}

write.table(exp_matrix_for_deg_mutnor,'files/01_SNV_CNV/GSEA/luad.APOBEC.mut.nor.exp.for.deg.txt',sep = '\t', row.names=T,col.names=T,quote=F,)
range(exp_matrix_for_deg_mutnor)

degs_mut_nor=mg_limma_DEG(exp_matrix_for_deg_mutnor,
                          c(rep('Mut',times=ncol(exp_matrix_M_mut)),rep('Normal',times=ncol(exp_matrix_M_nor))),
                          ulab = 'Mut',
                          dlab = 'Normal')
table(degs_mut_nor$DEG$adj.P.Val<0.05)
table(degs_mut_nor$DEG$adj.P.Val<0.01)

library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)

## 未过滤任何基因
degs_mut_nor_sig<-degs_mut_nor$DEG

degs_mut_nor_sig_gene<-rownames(degs_mut_nor_sig)

degs_gene_entz=bitr(degs_mut_nor_sig_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
degs_gene_entz <- dplyr::distinct(degs_gene_entz,SYMBOL,.keep_all=TRUE)

########## 合并
gene_df <- data.frame(logFC=degs_mut_nor_sig$logFC,
                      SYMBOL = rownames(degs_mut_nor_sig))
gene_df <- merge(gene_df,degs_gene_entz,by="SYMBOL")
head(gene_df)

geneList<-gene_df$logFC
names(geneList)=gene_df$ENTREZID 
head(geneList)

geneList=sort(geneList,decreasing = T)
head(geneList) ## 降序排列

## KEGG 通路gmt文件
kegmt<-read.gmt("origin_datas/h.all.v7.5.1.entrez.gmt")

hallmark<-GSEA(geneList,TERM2GENE = kegmt)
dim(hallmark@result)
library(enrichplot)
library(ggplot2)
dotplot(hallmark,color="pvalue")
dotplot(hallmark)

table(hallmark@result$p.adjust<0.05)
View(hallmark@result)
### .sign was reserved for the sign of NES (activated for NES > 0 and suppressed for NES < 0). We can pass split parameter and then it will apply showCategory by splitting the results. The following example plot 30 activated and 30 suppressed enriched disease sets.
dotplot(hallmark,split=".sign",showCategory=20)+facet_grid(~.sign)
# writeMatrix(hallmark@result,'files/文件/tcga.Mut-WT.GSEA.res.txt')

###################################################
############## 拷贝数变异分析
###################################################
library("ComplexHeatmap")
library("circlize")
get_CNV_Preprocess=function(df_cnv){
  df_cnv$`Gene Symbol`=gsub("\\..*","",df_cnv$`Gene Symbol`)
  rownames(df_cnv)=df_cnv$`Gene Symbol`
  
  library(TCGAutils)
  aliquot_id_to_submitter_id=UUIDtoBarcode(colnames(df_cnv)[-c(1:3)]
                                           ,from_type = 'aliquot_ids')
  colnames(df_cnv)[-c(1:3)]=aliquot_id_to_submitter_id[,2]
  colnames(df_cnv)=substr(colnames(df_cnv),1,15)
  df_cnv=df_cnv[,-which(duplicated(colnames(df_cnv)))]
  
  df_cnv=dplyr::distinct(df_cnv,`Gene Symbol`,.keep_all=TRUE)
  
  ensg=df_cnv$`Gene Symbol`
  
  idmap=bitr(ensg,fromType="ENSEMBL",toType="SYMBOL",OrgDb="org.Hs.eg.db")
  idmap=dplyr::distinct(idmap,ENSEMBL,.keep_all=TRUE)## 一个基因匹配到多个geneid,随机取一个
  
  cnv.inds=which(!is.na(idmap$SYMBOL)) ## 去掉没有注释到gene symbol的行
  idmap=idmap[cnv.inds,]
  df_cnv=df_cnv[idmap$ENSEMBL,]
  df_cnv$`Gene Symbol`=idmap$SYMBOL
  
  df_cnv=dplyr::distinct(df_cnv,`Gene Symbol`,.keep_all=TRUE)
  dim(df_cnv)
  
  # rownames(df_cnv)=df_cnv$`Gene Symbol`
  df_cnv=df_cnv[,-c(2:3)]
  return(df_cnv)
}

abbr<-"LUAD"
cnv.all=read.csv('origin_datas/TCGA/Merge_TCGA-LUAD_GeneLevelCopyNumber.txt',sep = '\t'
                 ,stringsAsFactors = F,header = T,check.names = F)
cnv.all=get_CNV_Preprocess(cnv.all)
dim(cnv.all)
cnv.all[1:4,1:5]

table(substr(colnames(cnv.all)[-1],13,15))

names(table(substr(colnames(cnv.all)[-1],13,15)))

length(which(substr(colnames(cnv.all)[-1],13,15)=="-01"))

tcga.cnv.sp.selected=colnames(cnv.all)[-1][(which(substr(colnames(cnv.all)[-1],13,15)=="-01"))]

cnv.all=cnv.all[,c("Gene Symbol",tcga.cnv.sp.selected)]
dim(cnv.all)

table(substr(colnames(cnv.all)[-1],13,15))

##############################################################
################# 基因拷贝数变异 在 转移肿瘤中的频率统计
##############################################################
get_CNV_Freq=function(df_cnv,genes_custom,genes_type=NULL){
  df_cnv=reshape2::melt(df_cnv,id.vars='Gene Symbol',measure.vars=colnames(cnv.all)[-1]
                        ,variable.name='Sample',value.name='CopyNum')
  head(df_cnv)
  cnv_frq=as.data.frame(matrix(data=0,nrow=1,ncol=5))
  colnames(cnv_frq)<-c('type','gene','CNV_gain','CNV_loss','none_CNV')
  r<-1
  for(i in 1:length(genes_custom))
  {
    cnv_gene<-df_cnv[which(df_cnv$`Gene Symbol` == genes_custom[i]),]
    total_counts<-dim(cnv_gene)[1]
    if(is.null(genes_type)){
      cnv_frq[r,1]<-NA
    }else{
      cnv_frq[r,1]<-genes_type[i]
    }
    cnv_frq[r,2]<-genes_custom[i]
    
    cnv_frq[r,3]<-(dim(cnv_gene[which(cnv_gene$CopyNum >0),])[1])/total_counts ## 在整个拷贝数数据集中
    cnv_frq[r,4]<-(dim(cnv_gene[which(cnv_gene$CopyNum <0),])[1])/total_counts ## CNV_loss
    cnv_frq[r,5]<-(total_counts/total_counts)-cnv_frq[r,3]-cnv_frq[r,4] ## none_CNV
    r<-r+1
  }
  cnv_frq<-cnv_frq[order(cnv_frq$CNV_loss),]
}

cnv.all[1:4,1:5]
table(substr(colnames(cnv.all[,-1]),14,15))

cnv_freq=get_CNV_Freq(df_cnv=cnv.all,genes_custom = APOBEC.genes$Gene_symbol)

library(tidyr)
cnv_freq_1 <- cnv_freq %>%
  pivot_longer(c(CNV_gain, CNV_loss,none_CNV), names_to = "CNV", values_to = "value")

cnv_freq_1$gene = factor(cnv_freq_1$gene, levels=cnv_freq$gene)
fig2d=ggplot(cnv_freq_1,aes(x=gene,y=value,fill=CNV))+
  geom_bar(position = "fill",stat="identity")+
  theme_bw()+scale_fill_manual(values = pal_lancet('lanonc')(9)[c(2,1,3)])+
  labs(x="APOBEC gene family",y="Frequency",fill="CNV")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="top")
fig2d
savePDF('PDFs/Fig2D.pdf',fig2d,width = 4,height = 4)
##################################################################
################################# 不同样本类型 之间 基因的比较
##################################################################
get_Exp_Compare=function(df_exp,genes_custom=APOBEC.genes$Gene_symbol){
  df_exp=data.frame(geneSymbol=rownames(df_exp),df_exp,stringsAsFactors = F,check.names = F)
  df_exp=reshape2::melt(df_exp,id.vars='geneSymbol',measure.vars=colnames(df_exp)[-1]
                        ,variable.name='Sample',value.name='Expr')
  head(df_exp)
  df_exp$sampleType=substr(df_exp$Sample,14,15)
  mode(df_exp$sampleType)='integer'
  
  
  head(sampleType.codes)
  
  df_exp$cancer_type=NA
  for(code in unique(df_exp$sampleType)){
    inds=which(df_exp$sampleType==code)
    df_exp$cancer_type[inds]=sampleType.codes[which(sampleType.codes$Code==code),3]
  }
  head(df_exp)
  df_exp=df_exp[,c(1,2,3,5)]
  if(!is.null(genes_custom)){
    df_exp=df_exp[df_exp$geneSymbol %in% genes_custom,]
    rownames(df_exp)=NULL
  }
  return(df_exp)
}

dim(luad.tcga.exp)
range(luad.tcga.exp)

table(substr(colnames(luad.tcga.exp),14,15))

tcga.exp.cmp=get_Exp_Compare(luad.tcga.exp,genes_custom = APOBEC.genes$Gene_symbol)

table(tcga.exp.cmp$cancer_type)
head(tcga.exp.cmp)

####### 原发 跟 正常样本之间的比较
tcga.exp.cmp=tcga.exp.cmp[tcga.exp.cmp$cancer_type %in% c('TP','NT'),]

tcga.exp.cmp$cancer_type[which(tcga.exp.cmp$cancer_type=='TP')]='Primary Solid Tumor'
# tcga.exp.cmp$cancer_type[which(tcga.exp.cmp$cancer_type=='TM')]='Metastatic'
tcga.exp.cmp$cancer_type[which(tcga.exp.cmp$cancer_type=='NT')]='Solid Tissue Normal'

library(ggpubr)
p <- ggboxplot(tcga.exp.cmp, x="geneSymbol", y="Expr", color = "cancer_type", 
               palette = pal_lancet('lanonc')(9)[c(2,1,3)]
               # , facet.by = "fa"
               , short.panel.labs = F
)
p=p+stat_compare_means(aes(group=cancer_type), label = "p.signif", method = "t.test")
# p=p+facet_wrap(~ fa, scales = "free",ncol=1) 
fig2f=p+theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

################################################
################ 不同 拷贝数 之间的基因表达的比较
################################################
get_CNV_Custom=function(df_cnv,df_exp,genes_custom,abbr){
  df_cnv=reshape2::melt(df_cnv,id.vars='Gene Symbol',measure.vars=colnames(df_cnv)[-1]
                        ,variable.name='Sample',value.name='CopyNum')
  head(df_cnv)
  
  df_cnv$Sample=as.character(df_cnv$Sample)
  df_cnv=df_cnv[df_cnv$`Gene Symbol` %in% genes_custom,]
  dim(df_cnv)
  
  exp.samples=colnames(df_exp)
  print(table(substr(exp.samples,14,15)))
  #############
  sampl.num=length(unique(df_cnv$Sample))##297个样本
  ######## 拷贝数 变异频率
  all.mt=rbind()
  dat.al=rbind()
  
  for(g in genes_custom){
    # g=genes[1]
    inds=which(df_cnv$`Gene Symbol`==g)
    sub_matrix<-df_cnv[inds,]
    dim(sub_matrix)
    
    all.mt=rbind(all.mt,
                 c(sum(sub_matrix$CopyNum>0)/sampl.num,
                   sum(sub_matrix$CopyNum<0)/sampl.num))
    ##### 分别统计每个样本的突变类型和基因的表达情况
    amplification.samp=sub_matrix[which(sub_matrix$CopyNum>0),2]
    deletion.samp=sub_matrix[which(sub_matrix$CopyNum<0),2]
    diploid.samp=setdiff(unique(df_cnv$Sample),c(amplification.samp,deletion.samp))
    
    length(c(amplification.samp,deletion.samp,diploid.samp))
    
    ######## CNV+RNA-seq均有数据的样本名称
    amplification.samp=intersect(amplification.samp,exp.samples)
    deletion.samp=intersect(deletion.samp,exp.samples)
    diploid.samp=intersect(diploid.samp,exp.samples)
    
    ###############
    ind1=match(deletion.samp,exp.samples)## del
    ind2=match(diploid.samp,exp.samples) ## amp
    ind3=match(amplification.samp,exp.samples) ## diploid
    ind4=grep('-11$',exp.samples) ## 原发组织样本
    
    sp=c(rep('deletion',length(ind1)),rep('diploid',length(ind2))
         ,rep('amplification',length(ind3)),rep('Normal',length(ind4)))
    # print(length(sp))
    samples=c(deletion.samp,diploid.samp,amplification.samp,exp.samples[grep('-11$',exp.samples)])
    ###############
    dat=data.frame(Expression=as.numeric(df_exp[which(row.names(df_exp)==g), c(ind1,ind2,ind3,ind4)])
                   ,category=sp,type=rep(g,length(sp)),sample=samples)
    dim(dat)
    dat.al=rbind(dat.al,dat)  ## 一行代表一个样本
  }
  dat.al=cbind(dat.al,fa=rep(abbr,nrow(dat.al)))
  rownames(all.mt)=genes_custom
  colnames(all.mt)=c('Amp','Del')
  return(dat.al)
}

table(substr(colnames(cnv.all)[-1],14,15))

tcga.exp.cmp.for_CNV=get_CNV_Custom(df_cnv=cnv.all,df_exp=luad.tcga.exp
                                    ,genes_custom = APOBEC.genes$Gene_symbol
                                    ,abbr = 'TCGA-LUAD')

library(ggpubr)
table(tcga.exp.cmp.for_CNV$category)

table(tcga.exp.cmp.for_CNV$category)

table(substr(tcga.exp.cmp.for_CNV$sample,14,15))

table(tcga.exp.cmp.for_CNV$category)
tcga.exp.cmp.for_CNV$category=factor(tcga.exp.cmp.for_CNV$category,levels = c('Normal','deletion','diploid','amplification'))
head(tcga.exp.cmp.for_CNV)
unique(tcga.exp.cmp.for_CNV$fa)

table(substr(tcga.exp.cmp.for_CNV$sample,14,15))

unique(names(which(table(substr(tcga.exp.cmp.for_CNV$sample,1,12))>1)))

tcga.exp.cmp.for_CNV$patient=substr(tcga.exp.cmp.for_CNV$sample,1,12)

########### 原发肿瘤样本中
p <- ggboxplot(tcga.exp.cmp.for_CNV, x="type", y="Expression", color = "category", 
               # palette = "jco",
               scales = "free_x",
               facet.by = "fa",ncol=3, short.panel.labs = T)
p=p+scale_color_manual(values = pal_lancet('lanonc')(9)[c(4,1,3,2)])+stat_compare_means(aes(group=category), label = "p.signif", method = "anova",
                                                                                      label.y = 10)
fig2e=p+theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
# fig2ef=mg_merge_plot(fig2e,fig2f,nrow = 2,ncol = 1,labels =c('E','F'),heights = c(1.2,1))
savePDF('PDFs/Fig2E.pdf',fig2e,width = 12,height = 6)

#############################################
############### Signature analysis
#############################################
library(maftools)
library('NMF')

######### 临床数据
luad.tcga.cli.forMut=luad.tcga.cli[,-1]
luad.tcga.cli.forMut$os=luad.tcga.cli.forMut$os/365

###########
tcga.project='LUAD'
maf=read.maf(paste0(MG_Grobal_DBPath,'/TCGA/Matrix/mutect2_mafs/TCGA.',tcga.project,'.mutect.somatic.maf')
             # ,clinicalData =luad.tcga.cli.forMut
             ,isTCGA = T)

dim(maf@variants.per.sample)## 565例原发肿瘤样本

oncoplot(maf = maf,top=20)

########## Estimate Tumor Mutation Burden
plotPoint=function(dat,variable=NULL,main=NULL,xlab=NULL,ylab=NULL){
  dat=dat[order(dat[,variable],decreasing = F),]
  medval=median(dat[,variable],na.rm = T)
  par(mar = c(2, 4, 2, 0))
  pointcol = "#009688"
  medlinecol = "#FF5722"
  yat = pretty(range(dat[variable != 0][, variable],finite=T,na.rm = T))
  plot(NA, xlim = c(0, nrow(dat)), ylim = range(yat), 
       axes = FALSE, xlab = NA, ylab = NA)
  abline(h = yat, lty = 2, lwd = 1, col = "gray90")
  abline(h = medval, lty = 2, lwd = 1, col = medlinecol)
  points(x = 1:nrow(dat), y = dat[,variable], 
         pch = 19, col = pointcol)
  title(main = main, sub = paste0("Median: ",sprintf("%.2f",medval),xlab), line = 0, adj = 0)
  axis(side = 2, at = yat, las = 2)
  mtext(text = ylab, side = 2, line = 2.5)
}
############
tcga.mutload=tmb(maf=maf)
tcga.mutload=data.frame(tcga.mutload)
tcga.mutload$total_perMB
median(tcga.mutload$total_perMB)
log10(median(tcga.mutload$total_perMB))

pdf('PDFs/Fig3A.pdf',height = 5,width = 5)
plotPoint(dat = tcga.mutload,variable='total_perMB',main='Mutation Burden',xlab ="/M",ylab = 'TMB/MB' )
dev.off()

########## Tumor heterogeneity and MATH scores
#计算mutant-allele tumor heterogeneity
tcga.MATH <- data.frame()
for (smp in unique(maf@data$Tumor_Sample_Barcode)){
  out.math = inferHeterogeneity(maf = maf, tsb = smp)
  Tumor_Sample_Barcode=unique(out.math$clusterData$Tumor_Sample_Barcode)
  MATH = unique(out.math$clusterData$MATH)
  out = data.frame(Tumor_Sample_Barcode, MATH)
  tcga.MATH = rbind(tcga.MATH, out)
}
head(tcga.MATH)

tcga.MATH$MATH_log=log2(tcga.MATH$MATH+1)
median(tcga.MATH$MATH)

# pdf('PDFs/Fig3B.pdf',height = 6,width = 6)
plotPoint(dat = tcga.MATH,variable='MATH',main='ITH Score',xlab =NULL,ylab = 'ITH Score' )
# dev.off()
########## 
maf.tnm = trinucleotideMatrix(maf = maf, prefix = NULL, add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
# plotApobecDiff(tnm = maf.tnm, maf = maf, pVal = 0.2)

table(maf.tnm$APOBEC_scores$APOBEC_Enriched)
dim(maf.tnm$APOBEC_scores)

### 567个样本
length(unique(c(as.character(maf@variants.per.sample$Tumor_Sample_Barcode),as.character(maf@maf.silent$Tumor_Sample_Barcode))))
setdiff(maf.tnm$APOBEC_scores$Tumor_Sample_Barcode,c(as.character(maf@variants.per.sample$Tumor_Sample_Barcode),as.character(maf@maf.silent$Tumor_Sample_Barcode)))

########
selected.column=c('Tumor_Sample_Barcode'
                  ,'tCw_to_G+tCw_to_T','n_C>G_and_C>T','n_bg_tcw','n_bg_wga','n_bg_C','n_bg_G'
                  ,'APOBEC_Enrichment','APOBEC_Enriched')
setdiff(selected.column,colnames(maf.tnm$APOBEC_scores))
head(maf.tnm$APOBEC_scores)

tcga.APOBEC.enrichment=data.frame(maf.tnm$APOBEC_scores,stringsAsFactors = F,check.names = F)[,selected.column]
rownames(tcga.APOBEC.enrichment)=tcga.APOBEC.enrichment$Tumor_Sample_Barcode

table(tcga.APOBEC.enrichment$APOBEC_Enrichment<=1)
table(tcga.APOBEC.enrichment$APOBEC_Enrichment<=2)
table(tcga.APOBEC.enrichment$APOBEC_Enrichment<=3)
table(tcga.APOBEC.enrichment$APOBEC_Enrichment<=4)

######## 重新划分
tcga.APOBEC.enrichment$AMES=ifelse(!is.na(tcga.APOBEC.enrichment$APOBEC_Enrichment),ifelse(tcga.APOBEC.enrichment$APOBEC_Enrichment<=1,"Low",ifelse(tcga.APOBEC.enrichment$APOBEC_Enrichment<=2,"Moderate","High")),NA)
tcga.APOBEC.enrichment$AMES=factor(tcga.APOBEC.enrichment$AMES,levels = c('Low','Moderate','High'))

table(tcga.APOBEC.enrichment$AMES)
table(tcga.APOBEC.enrichment$APOBEC_Enriched)

pdf('PDFs/Fig3B.pdf',height = 5,width = 5)
plotPoint(dat = tcga.APOBEC.enrichment,variable='tCw_to_G+tCw_to_T',main='TCW mutation',xlab =NULL,ylab = 'TCW mutation' )
dev.off()

pdf('PDFs/Fig3C.pdf',height = 5,width = 5)
plotPoint(dat = tcga.APOBEC.enrichment,variable='APOBEC_Enrichment',main='APOBEC enrichment scores',xlab =NULL,ylab = 'AMES' )
dev.off()

############ 预后生存差异
tcga.APOBEC_scores=merge(x=tcga.APOBEC.enrichment
                         , y = luad.tcga.cli.forMut
                         , by = 'Tumor_Sample_Barcode', sort = F, all.x = T)
tcga.APOBEC_scores=crbind2DataFrame(tcga.APOBEC_scores)

table(tcga.APOBEC_scores$APOBEC_Enrichment<=1)
table(tcga.APOBEC_scores$APOBEC_Enrichment>2)

tcga.APOBEC_scores$AMES=ifelse(!is.na(tcga.APOBEC_scores$APOBEC_Enrichment),ifelse(tcga.APOBEC_scores$APOBEC_Enrichment<=1,"Low",ifelse(tcga.APOBEC_scores$APOBEC_Enrichment<=2,"Moderate","High")),NA)

median(tcga.APOBEC_scores$APOBEC_Enrichment)

table(tcga.APOBEC_scores$AMES)
table(tcga.APOBEC_scores$APOBEC_Enriched)
library(survival)
ggplotKMCox(data.frame(time = tcga.APOBEC_scores$os
                             , event = tcga.APOBEC_scores$os_status
                             , groups=tcga.APOBEC_scores$APOBEC_Enriched)[!is.na(tcga.APOBEC_scores$os_status),]
                  , add_text = '')


ggplotKMCox(data.frame(time = tcga.APOBEC_scores$os
                       , event = tcga.APOBEC_scores$os_status
                       , groups=tcga.APOBEC_scores$AMES)[!is.na(tcga.APOBEC_scores$os_status),]
            , add_text = '')

##################################
dim(tcga.mutload)
dim(tcga.MATH)

tcga.mutload_use=merge(tcga.mutload,tcga.APOBEC.enrichment,by='Tumor_Sample_Barcode',sort=F)
rownames(tcga.mutload_use)=tcga.mutload_use$Tumor_Sample_Barcode
tcga.mutload_use$AMES=factor(tcga.mutload_use$AMES,levels = c('Low','Moderate','High'))

tcga.MATH_use=merge(tcga.MATH,tcga.APOBEC.enrichment,by='Tumor_Sample_Barcode',sort=F)
rownames(tcga.MATH_use)=tcga.MATH_use$Tumor_Sample_Barcode
tcga.MATH_use$AMES=factor(tcga.MATH_use$AMES,levels = c('Low','Moderate','High'))

fig3d1=mg_PlotMutiBoxplot(data=data.frame(AMES=tcga.mutload_use[,'total_perMB'])
                   , group = tcga.mutload_use$AMES
                   , legend.pos = 'tr'
                   , xangle = 0
                   , add = 'boxplot'
                   , xlab=''
                   , ylab = 'TMB'
                   , group_cols = ggsci::pal_aaas()(9)[c(3, 1, 2)]
                   , test_method = 'kruskal.test'
                   , fill = T)
   
fig3d2=mg_PlotMutiBoxplot(data=data.frame(AMES=tcga.MATH_use$MATH)
                   , group = tcga.MATH_use$AMES
                   , xangle = 0
                   , legend.pos = 'tr'
                   , add = 'boxplot'
                   , xlab=''
                   , ylab = 'ITH score'
                   , group_cols = ggsci::pal_aaas()(9)[c(3, 1, 2)]
                   , test_method = 'kruskal.test'
                   , fill = T
                   )

############# TNB 分析
tcga.TNB=read.table('origin_datas/TCIA/TCIA-NeoantigensData.LUAD.tsv',sep = '\t',header = T)
tcga.TNB=data.frame(table(tcga.TNB$patientBarcode))
colnames(tcga.TNB)=c('Tumor_Sample_Barcode','TNB')

tcga.TNB_use=merge(tcga.TNB,tcga.APOBEC.enrichment,by='Tumor_Sample_Barcode',sort=F)
rownames(tcga.TNB_use)=tcga.TNB_use$Tumor_Sample_Barcode
tcga.TNB_use$AMES=factor(tcga.TNB_use$AMES,levels = c('Low','Moderate','High'))

fig3d3=mg_PlotMutiBoxplot(data=data.frame(AMES=tcga.TNB_use$TNB)
                   , group = tcga.TNB_use$AMES
                   , legend.pos = 'tr'
                   , add = 'boxplot'
                   , xlab=''
                   , ylab = 'TNB'
                   , xangle=0
                   , group_cols = ggsci::pal_aaas()(9)[c(3, 1, 2)]
                   , test_method = 'kruskal.test'
                   , fill = T)

ggplot(tcga.mutload_use, aes(x=APOBEC_Enrichment, y=total_perMB)) + 
  xlab('APOBEC_Enrichment')+ylab('TMB')+
  geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =APOBEC_Enrichment, y =total_perMB))

ggplot(tcga.MATH_use, aes(x=APOBEC_Enrichment, y=MATH)) + 
  xlab('APOBEC_Enrichment')+ylab('ITH Score')+
  geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =APOBEC_Enrichment, y =MATH))

ggplot(tcga.TNB_use, aes(x=APOBEC_Enrichment, y=TNB)) + 
  xlab('APOBEC_Enrichment')+ylab('TNB')+
  geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =APOBEC_Enrichment, y =TNB))

############
library(ggplot2)
library(ggpubr)
library(ggExtra)

######### 同义突变
tcga.silent.per.sample=data.frame(table(maf@maf.silent$Tumor_Sample_Barcode))
colnames(tcga.silent.per.sample)=c('Tumor_Sample_Barcode','Variants')
tcga.silent.per.sample=merge(tcga.silent.per.sample,tcga.APOBEC.enrichment,by='Tumor_Sample_Barcode',sort=F)
tcga.silent.per.sample$Variants=log2(tcga.silent.per.sample$Variants)
######## 非同义突变
tcga.variants.per.sample_use=merge(maf@variants.per.sample,tcga.APOBEC.enrichment,by='Tumor_Sample_Barcode',sort=F)
tcga.variants.per.sample_use$Variants=log2(tcga.variants.per.sample_use$Variants)


tcga.all.mut.per.sample=rbind(data.frame(tcga.variants.per.sample_use,mut.type='Non-synonymous')
                              ,data.frame(tcga.silent.per.sample,mut.type='Synonymous'))

# corT=cor.test(x,y,method="spearman")
# cor=corT$estimate
# pValue=corT$p.value

fig3e=ggplot(tcga.all.mut.per.sample, aes(x=APOBEC_Enrichment, y=Variants,colour =mut.type)) + 
  xlab('AMES')+ylab('Mutation counts (log2)')+
  geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  scale_colour_manual(values = pal_lancet()(9)[c(7,8)])+
  theme(legend.position = c(0.8, 0.2))+
  stat_cor(method = 'spearman', aes(x =APOBEC_Enrichment, y =Variants))
fig3e

fig3d=mg_merge_plot(fig3d1,fig3d2,fig3d3,nrow = 1,ncol = 3,labels = c(''),common.legend = T)
fig3de=mg_merge_plot(fig3d,fig3e,nrow = 1,ncol = 2,labels = c('D','E'),widths = c(2,1))
savePDF('PDFs/Fig3DE.pdf',fig3de,width = 15,height = 5)
####################### 
setdiff(APOBEC.genes$Gene_symbol,rownames(luad.tcga.t.exp))

luad.tcga.t.exp.sub=luad.tcga.t.exp[APOBEC.genes$Gene_symbol,]
luad.tcga.t.exp.sub=t(luad.tcga.t.exp.sub)
library(ggcor)
p1=quickcor(luad.tcga.t.exp.sub, cor.test = TRUE) +
  geom_square(data = get_data(type = "lower", show.diag = FALSE)) +
  scale_fill_gradient2n(colours = rev(red_blue())) +
  geom_mark(data = get_data(type = "upper", show.diag = FALSE),size=2,sep = "\n") +
  geom_abline(slope = -1, intercept = 12) + remove_x_axis()

###############
luad.tcga.t.exp.sub_use=data.frame(Tumor_Sample_Barcode=substr(rownames(luad.tcga.t.exp.sub),1,12),luad.tcga.t.exp.sub,stringsAsFactors = F,check.names = F)
rownames(luad.tcga.t.exp.sub_use)=luad.tcga.t.exp.sub_use$Tumor_Sample_Barcode

tcga.mutload_use=crbind2DataFrame(tcga.mutload_use)
tcga.TNB_use=crbind2DataFrame(tcga.TNB_use)
tcga.MATH_use=crbind2DataFrame(tcga.MATH_use)
tcga.APOBEC.enrichment=crbind2DataFrame(tcga.APOBEC.enrichment)

tcga.smp.selected=Reduce(intersect,list(tcga.mutload_use$Tumor_Sample_Barcode
                      ,tcga.TNB_use$Tumor_Sample_Barcode
                      ,tcga.MATH_use$Tumor_Sample_Barcode
                      ,tcga.APOBEC.enrichment$Tumor_Sample_Barcode
                      ,luad.tcga.t.exp.sub_use$Tumor_Sample_Barcode))

setdiff(tcga.smp.selected,luad.tcga.t.exp.sub_use$Tumor_Sample_Barcode)

tcga.all.indicator=data.frame(TMB=tcga.mutload_use[tcga.smp.selected,'total_perMB']
                              ,TNB=tcga.TNB_use[tcga.smp.selected,'TNB']
                              ,ITH=tcga.MATH_use[tcga.smp.selected,'MATH']
                              ,AMES=tcga.APOBEC.enrichment[tcga.smp.selected,'APOBEC_Enrichment'])
head(tcga.all.indicator)

p2=quickcor(tcga.all.indicator,luad.tcga.t.exp.sub_use[tcga.smp.selected,2:ncol(luad.tcga.t.exp.sub_use)], cor.test = TRUE) +
  geom_square() +
  scale_fill_gradient2n(colours = rev(red_blue()))

fig3f=cowplot::plot_grid(p1, p2, ncol = 1
                         , align = "v"
                         , rel_heights = c(1.7,1))
fig3f

savePDF('PDFs/Fig3F.pdf',fig3f,height = 5,width = 5)
#########################################################################################
##################### Detecting cancer driver genes based on positional clustering
#########################################################################################
AMES.smp.low=tcga.APOBEC.enrichment$Tumor_Sample_Barcode[which(tcga.APOBEC.enrichment$AMES=='Low')]
AMES.smp.moderate=tcga.APOBEC.enrichment$Tumor_Sample_Barcode[which(tcga.APOBEC.enrichment$AMES=='Moderate')]
AMES.smp.high=tcga.APOBEC.enrichment$Tumor_Sample_Barcode[which(tcga.APOBEC.enrichment$AMES=='High')]

maf.AMES.low=subsetMaf(maf=maf,tsb = AMES.smp.low)
maf.AMES.moderate=subsetMaf(maf=maf,tsb = AMES.smp.moderate)
maf.AMES.high=subsetMaf(maf=maf,tsb = AMES.smp.high)

# maf.sig = oncodrive(maf = maf, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore')
# plotOncodrive(res = maf.sig, fdrCutOff = 0.05, useFraction = TRUE, labelSize = 0.5)
############
maf.sig.low = oncodrive(maf = maf.AMES.low, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore')
maf.sig.moderate = oncodrive(maf = maf.AMES.moderate, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore')
maf.sig.high = oncodrive(maf = maf.AMES.high, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore')

pdf('PDFs/Fig3G-1.pdf',height = 5,width = 4)
plotOncodrive(res = maf.sig.low, fdrCutOff = 0.05, useFraction = TRUE
              , labelSize = 0.5,bubbleSize=0.3)
dev.off()

pdf('PDFs/Fig3G-2.pdf',height = 5,width = 4)
plotOncodrive(res = maf.sig.moderate, fdrCutOff = 0.05, useFraction = TRUE
              , labelSize = 0.5,bubbleSize=0.3)
dev.off()

pdf('PDFs/Fig3G-3.pdf',height = 5,width = 4)
plotOncodrive(res = maf.sig.high, fdrCutOff = 0.05, useFraction = TRUE
              , labelSize = 0.5,bubbleSize=0.3)
dev.off()

######## KRAS 突变
maf.sub.drive=subsetMaf(maf=maf,genes='KRAS') 
tcga.oncodrive.dat=merge(maf.sub.drive@variants.per.sample,tcga.APOBEC.enrichment,by='Tumor_Sample_Barcode',all.x=TRUE)
tcga.oncodrive.stat=data.frame(table(tcga.oncodrive.dat$AMES))
colnames(tcga.oncodrive.stat)=c('AMES','Mutation frequency')
tcga.oncodrive.stat$`Mutation frequency`=tcga.oncodrive.stat$`Mutation frequency`/sum(tcga.oncodrive.stat$`Mutation frequency`)
tcga.oncodrive.stat
tcga.oncodrive.stat$AMES=factor(tcga.oncodrive.stat$AMES,levels = c('Low','Moderate','High'))

fig3h=ggbarplot(tcga.oncodrive.stat, x="AMES", y="Mutation frequency", fill = "AMES", color = "white", 
          # orientation = "horiz",   #横向显示
          # palette = "aaas",    #配色方案
          legend = "top",    #图例位置
          sort.val = "asc",    #上升排序，区别于desc
          sort.by.groups=TRUE)+    #按组排序
  scale_fill_manual(values = ggsci::pal_aaas()(9)[c(3, 1, 2)])+
  scale_y_continuous(expand=c(0, 0)) +
  scale_x_discrete(expand=c(0,0))
fig3h
savePDF('PDFs/Fig3H.pdf',fig3h,width = 5,height = 5)
################ fisher exact test
maf.AMES.Low.Moderate=merge_mafs(list(maf.AMES.low,maf.AMES.moderate))
# maf.AMES.High.Moderate=merge_mafs(list(maf.AMES.high,maf.AMES.moderate))

tcga.high.vs.low.moderate <- mafCompare(m1 = maf.AMES.Low.Moderate, m2 = maf.AMES.high, m1Name = 'L+M', m2Name = 'H', minMut = 5)
print(tcga.high.vs.low.moderate)

table(tcga.high.vs.low.moderate$results$pval<0.001)

pdf('PDFs/Fig3I.pdf',height = 5,width = 5)
forestPlot(mafCompareRes = tcga.high.vs.low.moderate, pVal = 0.001,color='black')
dev.off()

#############
maf.data=maf@data
maf.data=merge(maf.data,tcga.APOBEC.enrichment,by='Tumor_Sample_Barcode',all.x=TRUE)
genes.APOBEC.stat=table(maf.data$Hugo_Symbol,maf.data$AMES)
head(genes.APOBEC.stat,n=10)
genes.APOBEC.stat=as.matrix(genes.APOBEC.stat)
genes.APOBEC.stat=cbind(genes.APOBEC.stat,ave_mut=rowMeans(genes.APOBEC.stat[,1:3]))
genes.APOBEC.stat=data.frame(genes.APOBEC.stat)
genes.APOBEC.stat$group="all non-silent muations"
genes.APOBEC.stat=genes.APOBEC.stat[order(genes.APOBEC.stat$ave_mut,decreasing = T),]
head(genes.APOBEC.stat,n=10)
genes.APOBEC.stat$group[1:10]='Top 10 non-silent muations'

table(genes.APOBEC.stat$group)

writeMatrix(genes.APOBEC.stat,'files/genes.APOBEC.stat.txt',header = T)

# genes.APOBEC.stat=read.table('../files/genes.APOBEC.stat.txt',header = T,sep='\t',row.names = 1)
# genes.APOBEC.stat=genes.APOBEC.stat[order(genes.APOBEC.stat$ave_mut,decreasing = F),]
# range(genes.APOBEC.stat$ave_mut)
# colnames(genes.APOBEC.stat)
# library("ggtern")
# p1<-ggtern(data=genes.APOBEC.stat,aes(x=Low,y=Moderate,z=High))+geom_point(aes(size=0.01,color=group),alpha=0.8)
# p1=p1+scale_colour_manual(values = c("grey","red"))
# p1=p1+theme_legend_position(x="tl")+theme_bw()
# p1
# ggsave('../PDFs/Fig3J.pdf',p1,width = 4,height = 4)

head(genes.APOBEC.stat,n=10)
############ top10
somaticInteractions(maf = maf.AMES.low, top = 10, pvalue = c(0.05, 0.01))
somaticInteractions(maf = maf.AMES.moderate, top = 10, pvalue = c(0.05, 0.01))
somaticInteractions(maf = maf.AMES.high, top = 10, pvalue = c(0.05, 0.01))

############# Oncogenic Signaling Pathways
OncogenicPathways(maf = maf.AMES.low)
PlotOncogenicPathways(maf = maf.AMES.low, pathways = "RTK-RAS")

############# DDR pathways
DDR.pathways[,1]=gsub(" \\(.*\\)","",DDR.pathways[,1])
DDR.pathways=apply(DDR.pathways, 1,function(x){data.frame(Pathway=x[1],Gene=unlist(strsplit(x[2],split=", ")),stringsAsFactors = F)})
DDR.pathways=do.call(rbind,DDR.pathways)
DDR.pathways=DDR.pathways[,c('Gene','Pathway')]
head(DDR.pathways)

setdiff(DDR.pathways$Gene,maf.AMES.low@gene.summary$Hugo_Symbol)
setdiff(DDR.pathways$Gene,maf.AMES.moderate@gene.summary$Hugo_Symbol)
setdiff(DDR.pathways$Gene,maf.AMES.high@gene.summary$Hugo_Symbol)

setdiff(DDR.pathways$Gene,maf@gene.summary$Hugo_Symbol)

maf.onco.pathways.low=OncogenicPathways(maf = maf.AMES.low,pathways =DDR.pathways)
maf.onco.pathways.moderate=OncogenicPathways(maf = maf.AMES.moderate,pathways =DDR.pathways)
maf.onco.pathways.high=OncogenicPathways(maf = maf.AMES.high,pathways =DDR.pathways)

maf.onco.path.all.res=rbind(data.frame(maf.onco.pathways.low,type="Low")
                            ,data.frame(maf.onco.pathways.moderate,type="Moderate")
                            ,data.frame(maf.onco.pathways.high,type="High"))
maf.onco.path.all.res=maf.onco.path.all.res[,c(1,7,6)]
maf.onco.path.all.res$Fraction_mutated_samples=maf.onco.path.all.res$Fraction_mutated_samples*100
head(maf.onco.path.all.res)
maf.onco.path.all.res$Fraction_mutated_samples
sprintf("%.1f", maf.onco.path.all.res$Fraction_mutated_samples)
signif(maf.onco.path.all.res$Fraction_mutated_samples, 3)   

get_circularBarplot=function(data){
  ###############################
  ## 同时对分组和得分排序，可获得排序的圆环条形图
  colnames(data)=c('ID','group','value')
  data <- data %>% arrange(group,value)
  table(data$group)
  data$value=signif(data$value, 3)  
  ###############################
  ## 在每个分组数据的后面插入几行缺失值
  empty_bar <- 2
  to_add <- data.frame(matrix(NA, empty_bar*nlevels(data$group), ncol(data)))
  colnames(to_add) <- colnames(data) # 设置数据表的名称
  to_add$group <- rep(levels(data$group), each=empty_bar)#为数据表添加分组变量
  data <- rbind(data, to_add) # 合并两个数据
  data <- data %>% arrange(group) # 将数据根据分组进行排序
  data$id <- seq(1, nrow(data))
  # 获取每个样本的名称在y轴的位置和倾斜角度
  label_data <- data
  number_of_bar <- nrow(label_data) # 计算条的数量
  ## 每个条上标签的轴坐标的倾斜角度
  angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar 
  label_data$hjust <- ifelse( angle < -90, 1, 0) # 调整标签的对其方式
  label_data$angle <- ifelse(angle < -90, angle+180, angle) ## 标签倾斜角度
  ## 为数据准备基础弧线的数据
  base_data <- data %>% group_by(group) %>% 
    summarize(start=min(id), end=max(id) - empty_bar) %>% 
    rowwise() %>% mutate(title=mean(c(start, end)))
  # 为网格标尺准备数据
  grid_data <- base_data
  grid_data$end <- grid_data$end[c(nrow(grid_data), 1:nrow(grid_data)-1)] + 1
  grid_data$start <- grid_data$start - 1
  grid_data <- grid_data[-1,]
  
  ###############################
  ## 可视化分组圆环条形图
  p1 <- ggplot(data)+
    ## 添加条形图
    geom_bar(aes(x=as.factor(id), y=value, fill=group),stat="identity",
             alpha=0.8) +
    ##为条形图添加一些划分等级的线(20/40/60/80)(按比例添加是因为知道满分100)
    # geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80),
    #              colour = "blue", alpha=0.5, size=0.5 ,inherit.aes = FALSE)+
    # geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60),
    #              colour = "blue", alpha=0.5, size=0.5 ,inherit.aes = FALSE )+
    # geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40),
    #              colour = "blue", alpha=0.5, size=0.5 , inherit.aes = FALSE )+
    # geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20),
    #              colour = "blue", alpha=0.5, size=0.5 , inherit.aes = FALSE )+
    # # 添加文本表示(20/40/60/80)表示每条线的大小
    # annotate("text", x = rep(max(data$id),4), y = c(20, 40, 60, 80), 
    #          label = c("20", "40", "60", "80") , color="blue", size=3,
    #          angle=0, fontface="bold", hjust=1) +
    ylim(-50,110) + ## 设置y轴坐标表的取值范围,可流出更大的圆心空白
    ## 设置使用的主题并使用极坐标系可视化条形图
    theme_minimal() +
    theme(legend.position = "none", # 不要图例
          axis.text = element_blank(),# 不要x轴的标签
          axis.title = element_blank(), # 不要坐标系的名称
          panel.grid = element_blank(), # 不要网格线
          plot.margin = unit(rep(-3,4), "cm")
          )+ ## 整个图与周围的边距
    coord_polar() + ## 极坐标系
    ## 为条形图添加文本
    geom_text(data=label_data, 
              aes(x=id, y=value+5, label=ID,hjust=hjust),
              color="black",fontface="bold",alpha=0.8, size=2.5, 
              angle= label_data$angle, inherit.aes = FALSE) +
    geom_text(data=label_data, 
              aes(x=id, y=value-10, label=value,hjust=hjust),
              color="black",fontface="bold",alpha=0.8, size=2.5, 
              angle= label_data$angle, inherit.aes = FALSE) +
    
    # 为图像添加基础线的信息
    geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5),
                 colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )+
    ## 添加分组文本信息
    geom_text(data=base_data, aes(x = title, y = -18, label=group),alpha=0.8,
              colour = "black", size=4,fontface="bold", inherit.aes = FALSE)
  # p1=p1+scale_fill_brewer(palette = 'Set2')
  p1=p1+scale_fill_manual(values = ggsci::pal_aaas()(9)[c(3, 1, 2)])
  return(p1)
}

fig3j=get_circularBarplot(maf.onco.path.all.res)
##############
maf.DDR.pws=subsetMaf(maf = maf,genes=DDR.pathways$Gene)
maf.DDR.pws.smp=getSampleSummary(maf.DDR.pws)
dim(maf.DDR.pws.smp)
head(maf.DDR.pws.smp)
maf.DDR.pws.smp_use=merge(maf.DDR.pws.smp,tcga.APOBEC.enrichment,by='Tumor_Sample_Barcode',all.y=TRUE)
head(maf.DDR.pws.smp_use)
maf.DDR.pws.smp_use$DDR=ifelse(!is.na(maf.DDR.pws.smp_use$total),'Mut','WT')

fig3k=mg_PlotMutiBoxplot(data = data.frame(DDR=maf.DDR.pws.smp_use$APOBEC_Enrichment)
                   , group = maf.DDR.pws.smp_use$DDR
                   , legend.pos = 'tr'
                   , add = 'boxplot'
                   # , xlab='DDR'
                   , ylab = 'AMES'
                   , xangle=0
                   , group_cols = pal_nejm()(8)[c(1,4)]
                   , test_method = 'kruskal.test')

fig3jk=mg_merge_plot(fig3j,fig3k,nrow = 1,ncol = 2,labels = c('J','K'),widths = c(5,3))
savePDF('PDFs/Fig3JK.pdf',fig3jk,height = 5,width = 8)

DDR.genes=maf@gene.summary[maf@gene.summary$Hugo_Symbol %in% unique(DDR.pathways$Gene),]
DDR.genes=DDR.genes[order(DDR.genes$total,decreasing = T),]
DDR.genes.top10=DDR.genes$Hugo_Symbol[1:10]

pdf('PDFs/Fig3L.pdf',height = 5,width = 7)
oncoplot(maf,genes = DDR.genes.top10)
dev.off()

########################################################
############ Mutational Signatures
########################################################
# maf.sign = estimateSignatures(mat = maf.tnm, nTry = 6,nrun = 100
#                               # ,pConstant = 1e-12
#                               )
# save(maf.sign,file = 'origin_datas/TCGA/maf.sign.RData')
load('origin_datas/TCGA/maf.sign.RData')
pdf('PDFs/FigS1.pdf',height = 5,width = 5)
plotCophenetic(res = maf.sign)
dev.off()

maf.signatures = extractSignatures(mat = maf.tnm, n = 4)

maf.v3.cosm = compareSignatures(nmfRes = maf.signatures, sig_db = "SBS")
library('pheatmap')
pheatmap::pheatmap(mat = maf.v3.cosm$cosine_similarities
                   , cluster_rows = FALSE
                   , main = "cosine similarity against validated signatures")

pdf('PDFs/Fig4A.pdf',height = 12,width = 6)
plotSignatures(nmfRes = maf.signatures, title_size = 1.2, sig_db = "SBS")
dev.off()

rowSums(maf.signatures$contributions)/sum(maf.signatures$contributions)

########################
# corT=cor.test(x,y,method="spearman")
# cor=corT$estimate
# pValue=corT$p.value


maf.signatures.contri=maf.signatures$contributions

maf.signatures.contri_use=data.frame(Tumor_Sample_Barcode=colnames(maf.signatures.contri),t(maf.signatures.contri))
head(maf.signatures.contri_use)
maf.signatures.contri_use=merge(maf.signatures.contri_use,tcga.APOBEC.enrichment,by='Tumor_Sample_Barcode',all.y=TRUE)
colnames(maf.signatures.contri_use)

# p.all=list()
# for(var in c('Signature_1','Signature_2','Signature_3','Signature_4')){
#   # var='Signature_3'
#   p=ggplot(maf.signatures.contri_use, aes(x=APOBEC_Enrichment, y=get(var))) + 
#     xlab('AMES')+ylab(var)+
#     geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
#     stat_cor(method = 'spearman', aes(x =APOBEC_Enrichment, y =get(var)))
#   p.all=c(p.all,list(p))
# }
# length(p.all)

p1=ggplot(maf.signatures.contri_use, aes(x=APOBEC_Enrichment, y=Signature_1)) +
    xlab('AMES')+ylab('Signature_1')+
    geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
    stat_cor(method = 'spearman', aes(x =APOBEC_Enrichment, y =Signature_1))

p2=ggplot(maf.signatures.contri_use, aes(x=APOBEC_Enrichment, y=Signature_2)) +
  xlab('AMES')+ylab('Signature_2')+
  geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =APOBEC_Enrichment, y =Signature_2))

p3=ggplot(maf.signatures.contri_use, aes(x=APOBEC_Enrichment, y=Signature_3)) +
  xlab('AMES')+ylab('Signature_3')+
  geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =APOBEC_Enrichment, y =Signature_3))

p4=ggplot(maf.signatures.contri_use, aes(x=APOBEC_Enrichment, y=Signature_4)) +
  xlab('AMES')+ylab('Signature_4')+
  geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =APOBEC_Enrichment, y =Signature_4))

fig4b=mg_merge_plot(list(p1,p2,p3,p4),nrow = 4,ncol = 1,labels = c('B'))
fig4b
savePDF('PDFs/Fig4B.pdf',fig4b,height = 12,width = 4)

p1=ggplot(maf.signatures.contri_use, aes(x=APOBEC_Enrichment, y=Signature_2)) + 
  xlab('AMES')+ylab('Signature_2')+
  geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =APOBEC_Enrichment, y =Signature_2))
p1

p1=ggplot(maf.signatures.contri_use, aes(x=APOBEC_Enrichment, y=Signature_3)) + 
  xlab('AMES')+ylab('Signature_3')+
  geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =APOBEC_Enrichment, y =Signature_3))
p1

p1=ggplot(maf.signatures.contri_use, aes(x=APOBEC_Enrichment, y=Signature_4)) + 
  xlab('AMES')+ylab('Signature_4')+
  geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =APOBEC_Enrichment, y =Signature_4))
p1

##################### 
######### 差异表达分析
mg_limma_DEG_use=function(exp,group,ulab,dlab=NULL,paired=FALSE){
  library(limma)
  ind1=which(group==ulab)
  if(is.null(dlab)){
    ind2=which(group!=ulab)
  }else{
    ind2=which(group==dlab)
  }
  sml <- c(rep('G1',length(ind1)),rep('G0',length(ind2)))    # set group names
  eset=exp[,c(ind1,ind2)]
  fl <- as.factor(sml)
  if(paired){
    pairinfo=factor(rep(1:length(ind1),each=2))
    design <- model.matrix(~fl+0+pairinfo)
    colnames(design)[1:length(levels(fl))] <- levels(fl)
  }else{
    design <- model.matrix(~fl+0)
    colnames(design)[1:length(levels(fl))] <- levels(fl)
  }
  cont.matrix<-limma::makeContrasts(contrasts='G1-G0',levels=design)
  #print(head(eset))
  fit<-lmFit (eset,design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  #print(sml)
  tT <- topTable(fit2, adjust="fdr", coef=1,sort.by="B", number=nrow(eset))
  regulated=ifelse(tT$logFC>0,'Up','Down')
  lfcs=c(log2(1.2),log2(1.3),log2(1.5),1)
  all.deg.cnt=cbind()
  for(lfc in lfcs){
    deg1=regulated[which(abs(tT$logFC)>lfc&tT$P.Value<0.05)]
    deg2=regulated[which(abs(tT$logFC)>lfc&tT$P.Value<0.01)]
    deg3=regulated[which(abs(tT$logFC)>lfc&tT$adj.P.Val<0.05)]
    deg4=regulated[which(abs(tT$logFC)>lfc&tT$adj.P.Val<0.01)]
    all.deg.cnt=cbind(all.deg.cnt,c(paste0(sum(deg1=='Up'),'|',sum(deg1=='Down'))
                                    ,paste0(sum(deg2=='Up'),'|',sum(deg2=='Down'))
                                    ,paste0(sum(deg3=='Up'),'|',sum(deg3=='Down'))
                                    ,paste0(sum(deg4=='Up'),'|',sum(deg4=='Down'))))
  }
  row.names(all.deg.cnt)=c('p<0.05','p<0.01','FDR<0.05','FDR<0.01')
  colnames(all.deg.cnt)=paste0(c('1.2','1.3','1.5','2'),'-fold')
  return(list(Exp=eset,Group=group[c(ind1,ind2)],DEG=tT,Summary=all.deg.cnt))
}

get_DEG=function(df_deg,p.cutoff=0.05,logfc.cutoff=1){
  df.deg.res=df_deg$DEG
  df.deg.sig=df.deg.res[which(df.deg.res$adj.P.Val<p.cutoff & abs(df.deg.res$logFC)>logfc.cutoff),]
}

######### 通路分析
getGeneFC=function(gene.exp,group,ulab=ulab,dlab=dlab){
  degs_C1_C3=mg_limma_DEG(gene.exp, 
                          group,
                          ulab=ulab,
                          dlab=dlab)
  table(degs_C1_C3$DEG$adj.P.Val<0.05)
  table(degs_C1_C3$DEG$adj.P.Val<0.01)
  
  ## 未过滤任何基因
  degs_C1_C3_sig<-degs_C1_C3$DEG[which(degs_C1_C3$DEG$adj.P.Val <= 1),]
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  
  degs_C1_C3_sig_gene<-rownames(degs_C1_C3_sig)
  
  degs_gene_entz=bitr(degs_C1_C3_sig_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  degs_gene_entz <- dplyr::distinct(degs_gene_entz,SYMBOL,.keep_all=TRUE)
  
  gene_df <- data.frame(logFC=degs_C1_C3_sig$logFC,
                        SYMBOL = rownames(degs_C1_C3_sig))
  gene_df <- merge(gene_df,degs_gene_entz,by="SYMBOL")
  head(gene_df)
  
  geneList<-gene_df$logFC
  names(geneList)=gene_df$ENTREZID 
  head(geneList)
  
  geneList=sort(geneList,decreasing = T)
  head(geneList) ## 降序排列
  return(geneList)
}

getGeneFC_use=function(gene.exp,group,ulab=ulab,dlab=NULL){
  degs_C1_C3=mg_limma_DEG_use(gene.exp, 
                              group,
                              ulab=ulab,
                              dlab=dlab)
  table(degs_C1_C3$DEG$adj.P.Val<0.05)
  table(degs_C1_C3$DEG$adj.P.Val<0.01)
  
  ## 未过滤任何基因
  degs_C1_C3_sig<-degs_C1_C3$DEG[which(degs_C1_C3$DEG$adj.P.Val <= 1),]
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  
  degs_C1_C3_sig_gene<-rownames(degs_C1_C3_sig)
  
  degs_gene_entz=bitr(degs_C1_C3_sig_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  degs_gene_entz <- dplyr::distinct(degs_gene_entz,SYMBOL,.keep_all=TRUE)
  
  gene_df <- data.frame(logFC=degs_C1_C3_sig$logFC,
                        SYMBOL = rownames(degs_C1_C3_sig))
  gene_df <- merge(gene_df,degs_gene_entz,by="SYMBOL")
  head(gene_df)
  
  geneList<-gene_df$logFC
  names(geneList)=gene_df$ENTREZID 
  head(geneList)
  
  geneList=sort(geneList,decreasing = T)
  head(geneList) ## 降序排列
  return(geneList)
}

###########
dim(tcga.APOBEC.enrichment)
table(tcga.APOBEC.enrichment$AMES)

tcga.group=tcga.APOBEC.enrichment[,c('Tumor_Sample_Barcode','AMES','APOBEC_Enrichment')]
rownames(tcga.group)=paste0(rownames(tcga.group),"-01")
tcga.group$AMES=factor(tcga.group$AMES,levels = c('Low','Moderate','High'))
#############
tcga.t.exp=luad.tcga.t.exp
setdiff(rownames(tcga.group),colnames(luad.tcga.t.exp))
setdiff(rownames(tcga.group),colnames(luad.tcga.t.exp))

tcga.group=tcga.group[intersect(rownames(tcga.group),colnames(luad.tcga.t.exp)),]
dim(tcga.group)
tcga.t.exp.sub=tcga.t.exp[,rownames(tcga.group)]
dim(tcga.t.exp.sub)

table(tcga.group$AMES)

tcga.geneList=getGeneFC_use(gene.exp=tcga.t.exp.sub,group=tcga.group$AMES
                            ,ulab='High',dlab = 'Low')

## KEGG 通路gmt文件
gmt2list <- function(annofile){
  
  if (!file.exists(annofile)) {
    stop("There is no such a gmt file!")
  }
  
  if (tools::file_ext(annofile) == "xz") {
    annofile <- xzfile(annofile)
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
    close(annofile)
  } else if (tools::file_ext(annofile) == "gmt") {
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
  } else {
    stop ("Only gmt and gmt.xz are accepted for gmt2list")
  }
  
  y <- strsplit(x, "\t")
  names(y) <- sapply(y, `[[`, 1)
  
  annoList <- lapply(y, `[`, c(-1,-2))
  return(annoList)
}
kegmt=gmt2list('origin_datas/c2.cp.kegg.v7.5.1.entrez.gmt')

##########
library(fgsea)
# tcga.geneList_use=data.frame(entrezgene_id=names(tcga.geneList),log2FoldChange=tcga.geneList)
# head(tcga.geneList_use)
# tcga.geneList_use <- tcga.geneList_use %>% 
#   mutate(rank = rank(log2FoldChange,  ties.method = "random")) %>%
#   arrange(desc(rank))
# head(tcga.geneList_use)
# gene_list=tcga.geneList_use$log2FoldChange
# names(gene_list)=tcga.geneList_use$entrezgene_id
# head(gene_list)

fgseaRes <- fgsea(pathways = kegmt, 
                  stats =tcga.geneList ,
                  minSize=15,
                  maxSize=500,
                  nperm=1000)
head(fgseaRes)
sum(fgseaRes[, pval < 0.01])
sum(fgseaRes[, padj < 0.05])

topPathwaysUp <- fgseaRes[ES > 0][head(order(padj), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(padj), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

pdf('PDFs/Fig5A.pdf',height = 6,width = 12)
plotGseaTable(kegmt[topPathwaysUp], tcga.geneList, fgseaRes, 
              gseaParam = 0.5)
dev.off()

###################### 
### 免疫浸润分析
######### 免疫特征
get.IOBR.immu.format=function(luad.tcga.t.exp.cibersort){
  luad.tcga.t.exp.cibersort = data.frame(luad.tcga.t.exp.cibersort)
  rownames(luad.tcga.t.exp.cibersort) = luad.tcga.t.exp.cibersort$ID
  luad.tcga.t.exp.cibersort = luad.tcga.t.exp.cibersort[, -1]
  colnames(luad.tcga.t.exp.cibersort) = gsub('(.*)_.*', "\\1", colnames(luad.tcga.t.exp.cibersort))
  return(luad.tcga.t.exp.cibersort)
}

get_PlotMutiBoxplot=function(exp.immune.score,group,group.val=NULL,group_cols=NULL,xangle=45){
  if(is.null(group_cols)){
    p=mg_PlotMutiBoxplot(exp.immune.score[rownames(group),]
                         , group = group[,group.val]
                         , legend.pos = 'tr'
                         , add = 'boxplot'
                         , xangle=xangle
                         , ylab = 'Estimated Proportion'
                         , group_cols = pal_nejm()(8)[c(1,3,2,4)]
                         , test_method = 'kruskal.test'
    )
  }else{
    p=mg_PlotMutiBoxplot(exp.immune.score[rownames(group),]
                         , group = group[,group.val]
                         , legend.pos = 'tr'
                         , add = 'boxplot'
                         , xangle=xangle
                         , binwidth=0.02
                         , ylab = 'Estimated Proportion'
                         , group_cols = group_cols
                         , test_method = 'kruskal.test'
    )
  }
  
  return(p)
}

load('origin_datas/immune/tcga.exp.cibersort.RData')
load('origin_datas/immune/tcga.exp.estimate.RData')

tcga.exp.cibersort=get.IOBR.immu.format(tcga.exp.cibersort)
tcga.exp.estimate=get.IOBR.immu.format(tcga.exp.estimate)

tcga.group$group=tcga.group$APOBEC_Enrichment
tcga.group$group[which(tcga.group$group<1)]='<1'
tcga.group$group[which(tcga.group$group>=1 & tcga.group$group<2)]='1~2'
tcga.group$group[which(tcga.group$group>=2 & tcga.group$group<3)]='2~3'
tcga.group$group[which(tcga.group$group>=3 & tcga.group$group<4)]='3~4'
# tcga.group$group[which(tcga.group$group>=2 & tcga.group$group<4)]='2~4'

table(tcga.group$group)

get_PlotMutiBoxplot(tcga.exp.cibersort[rownames(tcga.group),1:22],group = tcga.group,group.val='group')
get_PlotMutiBoxplot(tcga.exp.estimate[rownames(tcga.group),1:3],group=tcga.group,group.val='group',xangle=0)

############
library(igraph)

df_M=cbind(AMES=tcga.group$APOBEC_Enrichment,tcga.exp.cibersort[rownames(tcga.group),c(1:4,6:22)])

df_cor = cor(df_M,method = 'spearman')##### 相关性
df_pval = corrplot::cor.mtest(df_M, conf.level = 0.95)$p ###pvalue

library(tidyverse)
### Pivot data from wide to long
g = pivot_longer(data=rownames_to_column(as.data.frame(df_cor),var = "from"),
                 cols = 2:(ncol(df_cor)+1), ## Columns to pivot into longer format
                 names_to = "to",
                 values_to = "cor")
gp = pivot_longer(data=rownames_to_column(as.data.frame(df_pval)),
                  cols = 2:(ncol(df_pval)+1),
                  names_to = "gene",
                  values_to = "p")
g$p = gp$p

############# 去掉相关性为1的节点
g = g[g$from!=g$to,] 
############# 
g$group = case_when(g$cor>0.1 & g$p<0.05 ~ "positive",
                    g$cor< (-0.1) & g$p<0.05 ~ "negative",
                    TRUE~"not" )
head(g)

####### 
head(g[,c(1,2,3,5)])
########## 节点
cibersort.celltype=readMatrix('origin_datas/cibersort.cell.type.txt',row = F)
cibersort.celltype_use=rbind(c('AMES','AMES','AMES')
                             ,cibersort.celltype)
cibersort.celltype_use$class=ifelse(cibersort.celltype_use$Group=='AMES','yes','no')
head(cibersort.celltype_use)

##########
network =  graph_from_data_frame(d=g[g$group!="not",c(1,2,3,5)]
                                 , vertices=cibersort.celltype_use
                                 , directed=F) 

############ 网络设置
########## 节点的设置
head(cibersort.celltype_use)
unique(cibersort.celltype_use$Group)

V(network)$color=c(pal_npg()(9)[-8],pal_lancet()(9)[5])[as.numeric(as.factor(V(network)$Group))]
V(network)$label=V(network)$abbr
V(network)$label.size=0.8
V(network)$label.color='black'

########## 边的设置
E(network)$arrow.mode <- 0 
E(network)$color=c("#2874C5","#f87669")[as.numeric(as.factor(E(network)$group))]
E(network)$width <- abs(E(network)$cor)*5

pdf('PDFs/Fig5B.pdf',height = 6,width = 6)
par(bg="white", mar=c(0,0,0,0))
plot(network,
     vertex.size=10,
     vertex.label.degree=0,
     layout=layout.circle,
     vertex.frame.color="transparent",
     edge.curved = 0.2)
legend("bottomright"
       , legend=levels(as.factor(V(network)$Group))  
       , col = c(pal_npg()(9)[-8],pal_lancet()(9)[5]) , bty = "o"
       , pch=20 
       , pt.cex = 2
       , cex=0.7
       , ncol = 2
       , text.col=c(pal_npg()(9)[-8],pal_lancet()(9)[5]) , horiz = F
       )
legend("bottomleft",
       c("positive", 
         "negative"),
       col = c("#FB9A99", "#8DA0CB"), bty="n", 
       # inset=c(0,-0.2),
       cex = 1, lty = 1, lwd = 2)
dev.off()

#############
fig5c=mg_PlotMutiBoxplot(data = data.frame(AMES=tcga.exp.cibersort[rownames(tcga.group),'Macrophages_M1'])
                   , group = tcga.group$AMES
                   , legend.pos = 'tr'
                   , add = 'boxplot'
                   , ylab = 'Macrophages M1 cell abundance(CIBERSORT)'
                   , xangle=0
                   , group_cols = pal_aaas()(9)[c(3, 1, 2)]
                   , test_method = 'kruskal.test')

fig5c
savePDF('PDFs/Fig5C.pdf',fig5c,height = 6,width = 4)

################## 分子亚型之间的通路分析
######### 差异表达分析
mg_limma_DEG_use=function(exp,group,ulab,dlab=NULL,paired=FALSE){
  library(limma)
  ind1=which(group==ulab)
  if(is.null(dlab)){
    ind2=which(group!=ulab)
  }else{
    ind2=which(group==dlab)
  }
  sml <- c(rep('G1',length(ind1)),rep('G0',length(ind2)))    # set group names
  eset=exp[,c(ind1,ind2)]
  fl <- as.factor(sml)
  if(paired){
    pairinfo=factor(rep(1:length(ind1),each=2))
    design <- model.matrix(~fl+0+pairinfo)
    colnames(design)[1:length(levels(fl))] <- levels(fl)
  }else{
    design <- model.matrix(~fl+0)
    colnames(design)[1:length(levels(fl))] <- levels(fl)
  }
  cont.matrix<-limma::makeContrasts(contrasts='G1-G0',levels=design)
  #print(head(eset))
  fit<-lmFit (eset,design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  #print(sml)
  tT <- topTable(fit2, adjust="fdr", coef=1,sort.by="B", number=nrow(eset))
  regulated=ifelse(tT$logFC>0,'Up','Down')
  lfcs=c(log2(1.2),log2(1.3),log2(1.5),1)
  all.deg.cnt=cbind()
  for(lfc in lfcs){
    deg1=regulated[which(abs(tT$logFC)>lfc&tT$P.Value<0.05)]
    deg2=regulated[which(abs(tT$logFC)>lfc&tT$P.Value<0.01)]
    deg3=regulated[which(abs(tT$logFC)>lfc&tT$adj.P.Val<0.05)]
    deg4=regulated[which(abs(tT$logFC)>lfc&tT$adj.P.Val<0.01)]
    all.deg.cnt=cbind(all.deg.cnt,c(paste0(sum(deg1=='Up'),'|',sum(deg1=='Down'))
                                    ,paste0(sum(deg2=='Up'),'|',sum(deg2=='Down'))
                                    ,paste0(sum(deg3=='Up'),'|',sum(deg3=='Down'))
                                    ,paste0(sum(deg4=='Up'),'|',sum(deg4=='Down'))))
  }
  row.names(all.deg.cnt)=c('p<0.05','p<0.01','FDR<0.05','FDR<0.01')
  colnames(all.deg.cnt)=paste0(c('1.2','1.3','1.5','2'),'-fold')
  return(list(Exp=eset,Group=group[c(ind1,ind2)],DEG=tT,Summary=all.deg.cnt))
}

get_DEG=function(df_deg,p.cutoff=0.05,logfc.cutoff=1){
  df.deg.res=df_deg$DEG
  df.deg.sig=df.deg.res[which(df.deg.res$adj.P.Val<p.cutoff & abs(df.deg.res$logFC)>logfc.cutoff),]
}

######### 通路分析
getGeneFC=function(gene.exp,group,ulab=ulab,dlab=dlab){
  degs_C1_C3=mg_limma_DEG(gene.exp, 
                          group$Cluster,
                          ulab=ulab,
                          dlab=dlab)
  table(degs_C1_C3$DEG$adj.P.Val<0.05)
  table(degs_C1_C3$DEG$adj.P.Val<0.01)
  
  ## 未过滤任何基因
  degs_C1_C3_sig<-degs_C1_C3$DEG[which(degs_C1_C3$DEG$adj.P.Val <= 1),]
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  
  degs_C1_C3_sig_gene<-rownames(degs_C1_C3_sig)
  
  degs_gene_entz=bitr(degs_C1_C3_sig_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  degs_gene_entz <- dplyr::distinct(degs_gene_entz,SYMBOL,.keep_all=TRUE)
  
  gene_df <- data.frame(logFC=degs_C1_C3_sig$logFC,
                        SYMBOL = rownames(degs_C1_C3_sig))
  gene_df <- merge(gene_df,degs_gene_entz,by="SYMBOL")
  head(gene_df)
  
  geneList<-gene_df$logFC
  names(geneList)=gene_df$ENTREZID 
  head(geneList)
  
  geneList=sort(geneList,decreasing = T)
  head(geneList) ## 降序排列
  return(geneList)
}

getGeneFC_use=function(gene.exp,group,ulab=ulab,dlab=NULL){
  degs_C1_C3=mg_limma_DEG_use(gene.exp, 
                          group$Cluster,
                          ulab=ulab,
                          dlab=dlab)
  table(degs_C1_C3$DEG$adj.P.Val<0.05)
  table(degs_C1_C3$DEG$adj.P.Val<0.01)
  
  ## 未过滤任何基因
  degs_C1_C3_sig<-degs_C1_C3$DEG[which(degs_C1_C3$DEG$adj.P.Val <= 1),]
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  
  degs_C1_C3_sig_gene<-rownames(degs_C1_C3_sig)
  
  degs_gene_entz=bitr(degs_C1_C3_sig_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  degs_gene_entz <- dplyr::distinct(degs_gene_entz,SYMBOL,.keep_all=TRUE)
  
  gene_df <- data.frame(logFC=degs_C1_C3_sig$logFC,
                        SYMBOL = rownames(degs_C1_C3_sig))
  gene_df <- merge(gene_df,degs_gene_entz,by="SYMBOL")
  head(gene_df)
  
  geneList<-gene_df$logFC
  names(geneList)=gene_df$ENTREZID 
  head(geneList)
  
  geneList=sort(geneList,decreasing = T)
  head(geneList) ## 降序排列
  return(geneList)
}

#################### 两两之间比较
luad.tcga.t.exp_use=tcga.t.exp.sub
luad.tcga.t.exp_use=luad.tcga.t.exp_use[rownames(luad.tcga.t.exp_use) %in% ann.pcg$gene_name,]
dim(luad.tcga.t.exp_use)
table(tcga.group$AMES)

tcga.deg.c1c2=mg_limma_DEG(luad.tcga.t.exp_use,group=tcga.group$AMES,ulab = 'High',dlab = 'Low')
tcga.deg.c1c2$Summary

tcga.deg.c1c2.sig=get_DEG(tcga.deg.c1c2,logfc.cutoff=log2(1.5),p.cutoff = 0.05)
nrow(tcga.deg.c1c2.sig)

fig6a=mg_volcano(logfc = tcga.deg.c1c2$DEG$logFC
                 ,tcga.deg.c1c2$DEG$adj.P.Val
                 ,cutFC=log2(1.5),cutPvalue=0.05)
dim(tcga.deg.c1c2$DEG)
writeMatrix(tcga.deg.c1c2$DEG,'files/文件/tcga.deg.high2low.txt')

###################### 功能富集分析
dim(tcga.deg.c1c2.sig)

tcga.deg.c1c2.sig.up=tcga.deg.c1c2.sig[which(tcga.deg.c1c2.sig$logFC>0),]
tcga.deg.c1c2.sig.down=tcga.deg.c1c2.sig[which(tcga.deg.c1c2.sig$logFC<0),]
nrow(tcga.deg.c1c2.sig.up)
nrow(tcga.deg.c1c2.sig.down)

################ 差异表达基因KEGG富集分析
tcga.c1c2.up.kegg.res=mg_clusterProfiler(rownames(tcga.deg.c1c2.sig.up))
table(tcga.c1c2.up.kegg.res$Enrich_tab$DB)

tcga.c1c2.down.kegg.res=mg_clusterProfiler(rownames(tcga.deg.c1c2.sig.down))

table(tcga.c1c2.up.kegg.res$Enrich_tab$DB)
table(tcga.c1c2.down.kegg.res$Enrich_tab$DB)

tcga.up.kegg=tcga.c1c2.up.kegg.res$Enrich_tab
tcga.down.kegg=tcga.c1c2.down.kegg.res$Enrich_tab

writeMatrix(dat = tcga.up.kegg,outpath = 'files/文件/tcga.up.kegg.txt')
writeMatrix(dat = tcga.down.kegg,outpath = 'files/文件/tcga.down.kegg.txt')

tcga.up.kegg_use=tcga.up.kegg
tcga.down.kegg_use=tcga.down.kegg

tcga.up.kegg_use$DB[which(tcga.up.kegg_use$DB=='pathway_KEGG')]='KEGG'
tcga.up.kegg_use$DB[which(tcga.up.kegg_use$DB=='geneontology_Biological_Process')]='GO_BP'
tcga.up.kegg_use$DB[which(tcga.up.kegg_use$DB=='geneontology_Cellular_Component')]='GO_CC'
tcga.up.kegg_use$DB[which(tcga.up.kegg_use$DB=='geneontology_Molecular_Function')]='GO_MF'

tcga.down.kegg_use$DB[which(tcga.down.kegg_use$DB=='pathway_KEGG')]='KEGG'
tcga.down.kegg_use$DB[which(tcga.down.kegg_use$DB=='geneontology_Biological_Process')]='GO_BP'
tcga.down.kegg_use$DB[which(tcga.down.kegg_use$DB=='geneontology_Cellular_Component')]='GO_CC'
tcga.down.kegg_use$DB[which(tcga.down.kegg_use$DB=='geneontology_Molecular_Function')]='GO_MF'

db.color=pal_aaas(alpha = 0.8)(9)[1:4]
names(db.color)=c('KEGG','GO_BP','GO_CC','GO_MF')

##########
tcga.up.kegg_forVis=tcga.up.kegg_use %>% 
  group_by(DB) %>%
  top_n(n=10)
tcga.up.kegg_forVis=data.frame(tcga.up.kegg_forVis,check.names = F,stringsAsFactors = F)
tcga.up.kegg_forVis$DB=factor(tcga.up.kegg_forVis$DB,levels = c('GO_BP','GO_CC','GO_MF','KEGG'))
tcga.up.kegg_forVis$FDR=-log10(tcga.up.kegg_forVis$FDR)

##########
tcga.down.kegg_forVis=tcga.down.kegg_use %>% 
  group_by(DB) %>%
  top_n(n=10)
tcga.down.kegg_forVis=data.frame(tcga.down.kegg_forVis,check.names = F,stringsAsFactors = F)
tcga.down.kegg_forVis$DB=factor(tcga.down.kegg_forVis$DB,levels = c('GO_BP','GO_CC','GO_MF','KEGG'))
tcga.down.kegg_forVis$FDR=-log10(tcga.down.kegg_forVis$FDR)

p1=ggbarplot(tcga.up.kegg_forVis, x="description", y="FDR", fill = "DB", color = "white", 
          orientation = "horiz",   #横向显示
          palette = db.color,    #配色方案
          legend = "right",    #图例位置
          sort.val = "asc",    #上升排序，区别于desc
          sort.by.groups=TRUE)+    #按组排序
  scale_y_continuous(expand=c(0, 0)) + scale_x_discrete(expand=c(0,0)) +
  labs(title = "Up-regulated")

p2=ggbarplot(tcga.down.kegg_forVis, x="description", y="FDR", fill = "DB", color = "white", 
                orientation = "horiz",   #横向显示
                palette = db.color,    #配色方案
                legend = "right",    #图例位置
                sort.val = "asc",    #上升排序，区别于desc
                sort.by.groups=TRUE)+    #按组排序
  scale_y_continuous(expand=c(0, 0)) + scale_x_discrete(expand=c(0,0)) +
  # theme(legend.position = c(1,0))+
  labs(title = "Down-regulated")
fig6b=mg_merge_plot(p1,p2,nrow = 1,ncol = 2,common.legend = T)

fig6ab=mg_merge_plot(fig6a,fig6b,nrow = 1,ncol = 2,labels = c('A','B'),widths = c(1,3))
savePDF('PDFs/Fig6AB.pdf',fig6ab,height = 5,width = 15)

########### 蛋白互作网络
writeMatrix(rownames(tcga.deg.c1c2.sig),outpath = 'files/tcga.deg.forPPI.txt',header=F)

tcga.ppi.hub.genes=read.csv('files/PPI/string_interactions_Cluster1_node.csv')
tcga.ppi.hub.genes=tcga.ppi.hub.genes$name
length(tcga.ppi.hub.genes)

tcga.ppi.hub.genes.kegg=mg_clusterProfiler(tcga.ppi.hub.genes)

tcga.ppi.hub.genes.kegg_use=tcga.ppi.hub.genes.kegg$Enrich_tab

tcga.ppi.hub.genes.kegg_use$DB[which(tcga.ppi.hub.genes.kegg_use$DB=='pathway_KEGG')]='KEGG'
tcga.ppi.hub.genes.kegg_use$DB[which(tcga.ppi.hub.genes.kegg_use$DB=='geneontology_Biological_Process')]='GO_BP'
tcga.ppi.hub.genes.kegg_use$DB[which(tcga.ppi.hub.genes.kegg_use$DB=='geneontology_Cellular_Component')]='GO_CC'
tcga.ppi.hub.genes.kegg_use$DB[which(tcga.ppi.hub.genes.kegg_use$DB=='geneontology_Molecular_Function')]='GO_MF'

tcga.ppi.hub.genes.kegg_forVis=tcga.ppi.hub.genes.kegg_use %>% 
  group_by(DB) %>%
  top_n(n=10)
tcga.ppi.hub.genes.kegg_forVis=data.frame(tcga.ppi.hub.genes.kegg_forVis,check.names = F,stringsAsFactors = F)
tcga.ppi.hub.genes.kegg_forVis$DB=factor(tcga.ppi.hub.genes.kegg_forVis$DB,levels = c('GO_BP','GO_CC','GO_MF','KEGG'))
tcga.ppi.hub.genes.kegg_forVis$FDR=-log10(tcga.ppi.hub.genes.kegg_forVis$FDR)

fig6d=ggbarplot(tcga.ppi.hub.genes.kegg_forVis, x="description", y="FDR", fill = "DB", color = "white", 
             orientation = "horiz",   #横向显示
             palette = db.color,    #配色方案
             legend = "right",    #图例位置
             sort.val = "asc",    #上升排序，区别于desc
             sort.by.groups=TRUE)+    #按组排序
  scale_y_continuous(expand=c(0, 0)) + scale_x_discrete(expand=c(0,0)) +
  # theme(legend.position = c(1,0))+
  labs(title = "Cluster 1")
savePDF('PDFs/Fig6D.pdf',fig6d,height = 6,width = 8)

##########################################
################
luad.tcga.t.exp.os=luad.tcga.t.exp.os[rownames(tcga.group),]
luad.tcga.t.exp.cli.tnm=luad.tcga.t.exp.cli.tnm[rownames(tcga.group),]
###########################
all.deg.sig=rownames(tcga.deg.c1c2.sig)
length(all.deg.sig)
writeMatrix(all.deg.sig,'files/文件/all.deg.sig.txt.txt')
########## 单因素Cox分析
tcga.cox=cox_batch(t(scale(t(luad.tcga.t.exp_use[(all.deg.sig),])))
                   ,time = luad.tcga.t.exp.os$OS.time
                   ,event = luad.tcga.t.exp.os$OS)
table(tcga.cox$p.value<0.05)
table(tcga.cox$p.value<0.01)
table(tcga.cox$p.value<0.001)

############
p.cutoff=0.05
tcga.cox_use=tcga.cox
tcga.cox_use$coef=log(tcga.cox_use$HR)
tcga.cox_use$Gene=rownames(tcga.cox_use)
tcga.cox_use$type=rep('None',nrow(tcga.cox_use))
tcga.cox_use$type[which(tcga.cox_use$p.value<p.cutoff & tcga.cox_use$coef>0)]='Risk'
tcga.cox_use$type[which(tcga.cox_use$p.value<p.cutoff & tcga.cox_use$coef<0)]='Protective'

table(tcga.cox_use$type)

p1 <- ggplot(data = tcga.cox_use,
             aes(x = coef,
                 y = -log10(p.value)))
p1=p1+geom_point(alpha=0.4, size=3.5, aes(color=type))
p1=p1+scale_color_manual(values=c(mg_colors[2],'grey',mg_colors[1]),limits = c("Protective",'None', "Risk"),name='State')
p1=p1+geom_hline(yintercept = -log10(p.cutoff),lty=4,col="black",lwd=0.8)
p1=p1+ylab('-log10(pvalue)')+xlab('Cox coefficient')
# p1=p1+ggrepel::geom_text_repel(data=module.genes.cox_use[which(module.genes.cox_use$p.value<0.05),],aes(label=Gene))
p1=p1+theme_bw()
p1=p1+theme(
  axis.text.y=element_text(family="Times",face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
  axis.title.y=element_text(family="Times",face="plain"), #设置y轴标题的字体属性
  legend.text=element_text(face="plain", family="Times", colour="black"  #设置图例的子标题的字体属性
  ),
  legend.title=element_text(face="plain", family="Times", colour="black" #设置图例的总标题的字体属性
  ),
  legend.justification=c(1,1), legend.position='top'
  ,legend.background = element_rect(fill = NA, colour = NA)
)

table(tcga.cox_use$type)

fig7a=p1
fig7a
dim(tcga.cox_use)
head(tcga.cox_use)

tcga.cox_forVis=tcga.cox_use
tcga.cox_forVis=tcga.cox_forVis[which(tcga.cox_forVis$type %in% c('Risk','Protective')),]
tcga.cox_forVis$p.value=-log10(tcga.cox_forVis$p.value)
range(tcga.cox_forVis$p.value)

library(ggpubr)
fig7a=ggdotchart(tcga.cox_forVis
           , x="Gene", y="coef",color = "type",
           palette = "aaas",     #配色方案
           legend = "right",     #图例位置
           sorting = "descending",   #上升排序，区别于desc
           add = "segments",    #增加线段
           rotate = TRUE,       #横向显示
           dot.size = tcga.cox_forVis$p.value/2,        #圆圈大小
           # label = round(tcga.cox_use$coef),   #圆圈内数值
           # font.label = list(color="white",size=10, vjust=0.5),   #圆圈内数值字体设置
           ggtheme = theme_pubr())+font("y.text", size = 10, color = "black")+
  theme(legend.position = 'top')+ylab('Cox coefficient')+xlab('')
fig7a
savePDF('PDFs/Fig7A.pdf',fig7a,height = 8,width = 6)

######### lasso
tcga.gene.sig=rownames(tcga.cox)[which(tcga.cox$p.value<p.cutoff)]
length(tcga.gene.sig)

tcga.exp.sig=luad.tcga.t.exp_use[tcga.gene.sig,]
tcga.exp.sig=t(tcga.exp.sig)
dim(tcga.exp.sig)

mg_lasso_cox_use=function(dat,time,event,nfolds=3,lambda.min=T,show_text=T,figLabels=c('A','B')){
  library("glmnet") 
  library('survival')
  t.inds=which(!is.na(time)&!is.na(event)&time>0)
  dat=dat[t.inds,]
  time=as.numeric(time[t.inds])
  event=as.numeric(event[t.inds])
  y=Surv(time,event)
  set.seed(888826666)
  # set.seed(123456789)
  # set.seed(1111111111.33333)
  fit1_cv = cv.glmnet(as.matrix(dat), y, family = "cox", nfolds=nfolds)
  fit<-glmnet(dat, y, family = "cox")
  if(lambda.min){
    lambda=fit1_cv$lambda.min
  }else{
    lambda=fit1_cv$lambda.1se
  }
  coefficients<-coef(fit,s=lambda)
  Active.Index<-which(coefficients[,1]!=0)
  genes=row.names(coefficients)[Active.Index]
  Active.coefficients<-coefficients[Active.Index]  
  g=mg_plot_lasso(fit,fit1_cv,lambda = lambda,show_text=show_text,figLabels=figLabels)
  return(list(Mode1=fit,Model2=fit1_cv,Genes=genes,Coef=Active.coefficients,lambda=lambda,plot=g))
}
options(ggrepel.max.overlaps = Inf)
tcga.lasso.res=mg_lasso_cox_use(t((t(tcga.exp.sig)))
                                , time = luad.tcga.t.exp.os$OS.time
                                , event = luad.tcga.t.exp.os$OS
                                , nfolds = 10
                                , lambda.min = T
                                , figLabels=c('B','C'))
tcga.lasso.res$Genes

tcga.lasso.res$lambda

tcga.lasso.res$plot

tcga.exp.for.cox=luad.tcga.t.exp_use[match(tcga.lasso.res$Genes,row.names(luad.tcga.t.exp_use)),]
dim(tcga.exp.for.cox)

lst.modl=createCoxModel((t(tcga.exp.for.cox))
                        ,time = luad.tcga.t.exp.os$OS.time
                        ,event = luad.tcga.t.exp.os$OS
                        ,isStep = T)
lst.modl$Genes
lst.modl$fmla

tcga.risk.score=lst.modl$Score
# tcga.risk.score=mosaic::zscore(tcga.risk.score)

lst.modl$Coef

gene.coef=data.frame(Gene=lst.modl$Genes,Coef=lst.modl$Coef)
gene.coef$Type=ifelse(lst.modl$Coef>0,'Risk','Protective')
gene.coef$Type=factor(gene.coef$Type,levels=c('Risk','Protective'))
table(gene.coef$Type)
library(dplyr)
fig7d=gene.coef %>% 
  ggplot(aes(reorder(Gene, Coef), Coef)) +
  geom_col(aes(fill = Type)) +
  coord_flip() +
  labs(x = "") +
  labs(y = "LASSO Cox coefficient") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1),legend.position = 'top')
fig7d

#################
cox_as_data_frame=function(coxphsummary){
  coxphsummary <- summary(coxphsummary)
  allHRCIs <- matrix(ncol = 4, nrow = 0)
  colnames(allHRCIs) <- c("HR", "HR.95L", "HR.95H", "pvalue")
  for (i in 1:dim(coxphsummary$coefficients)[1])
  {
    hrcis     <- matrix(
      c( exp(coxphsummary$coefficients[i, "coef"] + c(0, qnorm(0.025), qnorm(0.975)) * coxphsummary$coefficients[i, "se(coef)"]),
         # exp(-coxphsummary$coefficients[i, "coef"] + c(qnorm(0.025), 0, qnorm(0.975)) * coxphsummary$coefficients[i, "se(coef)"]),
         coxphsummary$coefficients[i, "Pr(>|z|)"]),
      nrow = 1, byrow = TRUE)
    allHRCIs <- rbind(allHRCIs, hrcis)
  }
  coefficient_labels <- rownames(coxphsummary$coefficients)
  rownames(allHRCIs)=coefficient_labels
  allHRCIs=data.frame(allHRCIs,check.names = F,stringsAsFactors = F)
  return(allHRCIs)
}

mg_Forestplot=function(df_m,outFile,width=6,height=3){
  gene=rownames(df_m)
  hr=sprintf("%.3f",df_m$"HR")
  hrLow=sprintf("%.3f",df_m$"HR.95L")
  hrHigh=sprintf("%.3f",df_m$"HR.95H")
  Hazard.ratio=paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal=ifelse(df_m$pvalue<0.001, "<0.001", sprintf("%.3f", df_m$pvalue))
  #########
  pdf(file=outFile, width = width, height =height,onefile = FALSE)
  n=nrow(df_m)
  nRow=n+1
  ylim=c(1,nRow)
  layout.show(layout(matrix(c(1,2),nc=2),width=c(3,2)))
  #森林图左边的基因信息
  xlim = c(0,3)
  par(mar=c(4,2,1.5,1.5))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
  
  #绘制森林图
  par(mar=c(4,1,1.5,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.03,col="black",lwd=2.5,lty=5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'blue')
  points(as.numeric(hr), n:1, pch = 20, col = boxcolor, cex=1.3)
  axis(1)
  dev.off()
}

lst.modl.cox.res=cox_as_data_frame(lst.modl$Cox)
rownames(lst.modl.cox.res)=lst.modl$Genes

mg_Forestplot(lst.modl.cox.res,outFile='PDFs/Fig7D.pdf',width = 6,height = 4)

fig7=mg_merge_plot(mg_merge_plot(fig7a,fig7d,nrow = 2,ncol = 1,labels = c('A','D'))
                   ,tcga.lasso.res$plot,nrow = 1,ncol = 2,heights = c(1,1))

savePDF('PDFs/Fig7.pdf',fig7,height = 12,width = 12)

####
tcga.cutoff <- survminer::surv_cutpoint(data.frame(time=luad.tcga.t.exp.os$OS.time/365,
                                                   event=luad.tcga.t.exp.os$OS,
                                                   risk=tcga.risk.score), time = "time", event = "event",
                                        variables = c("risk"))
tcga.cutoff=tcga.cutoff$cutpoint$cutpoint
# tcga.cutoff=0
# tcga.cutoff=median(tcga.risk.score)

fig8A=plotCoxModel_Batch_use(riskScore = tcga.risk.score
                             ,dat = t(tcga.exp.for.cox[match(lst.modl$Genes, row.names(tcga.exp.for.cox)), ])
                             , time = luad.tcga.t.exp.os$OS.time/365
                             , event = luad.tcga.t.exp.os$OS
                             , cutoff = tcga.cutoff
                             , title='AMrs'
                             , mks = c(1, 3, 5))
fig8a=fig8A[[1]]
fig8a

plotRiskScoreModel_use(riskScore = tcga.risk.score
                       ,dat = t(tcga.exp.for.cox[match(lst.modl$Genes, row.names(tcga.exp.for.cox)), ])
                       , time = luad.tcga.t.exp.os$OS.time/365
                       , event = luad.tcga.t.exp.os$OS
                       , cutoff = tcga.cutoff)

tcga.subtype.AMrs=ifelse(tcga.risk.score>tcga.cutoff,'High','Low')
tcga.subtype.AMrs=data.frame(tcga.subtype.AMrs)
colnames(tcga.subtype.AMrs)='Cluster'
table(tcga.subtype.AMrs$Cluster)
tcga.subtype.AMrs$Cluster=factor(tcga.subtype.AMrs$Cluster,levels = c('Low','High'))

######### GSE72094
match(lst.modl$Genes,row.names(gse72094.t.exp))
lst.modl$Genes

rownames(gse72094.t.exp)[grep("C7orf68",rownames(gse72094.t.exp))]
rownames(gse72094.t.exp)[which(rownames(gse72094.t.exp)=='C7orf68')]='HILPDA'

gse72094.model.dat=gse72094.t.exp[intersect(lst.modl$Genes,row.names(gse72094.t.exp)),]

lst.vd.mod1=createCoxModel((t(gse72094.model.dat)),time=gse72094.t.cli.os$OS.time/365,event = gse72094.t.cli.os$OS)
# gse72094.risk.score=lst.vd.mod1$Score
gse72094.risk.score=predictScoreByCoxModel(coxModel = lst.modl,(t(gse72094.model.dat)))
# gse72094.risk.score=mosaic::zscore(gse72094.risk.score)

lst.modl$fmla
lst.vd.mod1$fmla

gse72094.cutoff <- survminer::surv_cutpoint(data.frame(time=gse72094.t.cli.os$OS.time/365,
                                                       event=gse72094.t.cli.os$OS,
                                                       risk=gse72094.risk.score), time = "time", event = "event",
                                            variables = c("risk"))
gse72094.cutoff=gse72094.cutoff$cutpoint$cutpoint
# gse72094.cutoff=0

fig8B=plotCoxModel_Batch_use(riskScore = gse72094.risk.score
                             , dat = t(gse72094.t.exp[intersect(lst.modl$Genes, row.names(gse72094.t.exp)),])
                             , time = gse72094.t.cli.os$OS.time/365
                             , event = gse72094.t.cli.os$OS
                             , cutoff = gse72094.cutoff
                             , title='AMrs'
                             , mks = c(1,2,3,5))
fig8b=fig8B[[2]]
fig8b

gse72094.subtype.AMrs=ifelse(gse72094.risk.score>gse72094.cutoff,'High','Low')
gse72094.subtype.AMrs=data.frame(gse72094.subtype.AMrs,stringsAsFactors = F)
colnames(gse72094.subtype.AMrs)='Cluster'
table(gse72094.subtype.AMrs)
gse72094.subtype.AMrs$Cluster=factor(gse72094.subtype.AMrs$Cluster,levels = c('Low','High'))

######### GSE31210
match(lst.modl$Genes,row.names(gse31210.t.exp))
lst.modl$Genes

rownames(gse31210.t.exp)[grep("C7orf68",rownames(gse31210.t.exp))]
rownames(gse31210.t.exp)[which(rownames(gse31210.t.exp)=='C7orf68')]='HILPDA'

gse31210.model.dat=gse31210.t.exp[intersect(lst.modl$Genes,row.names(gse31210.t.exp)),]

lst.vd.mod1=createCoxModel((t(gse31210.model.dat)),time=gse31210.t.cli.os$OS.time/365,event = gse31210.t.cli.os$OS)
# gse31210.risk.score=lst.vd.mod1$Score
gse31210.risk.score=predictScoreByCoxModel(coxModel = lst.modl,t(gse31210.model.dat))
# gse31210.risk.score=mosaic::zscore(gse31210.risk.score)

lst.modl$fmla
lst.vd.mod1$fmla

gse31210.cutoff <- survminer::surv_cutpoint(data.frame(time=gse31210.t.cli.os$OS.time/365,
                                                   event=gse31210.t.cli.os$OS,
                                                   risk=gse31210.risk.score), time = "time", event = "event",
                                        variables = c("risk"))
gse31210.cutoff=gse31210.cutoff$cutpoint$cutpoint
# gse31210.cutoff=0

fig8C=plotCoxModel_Batch_use(riskScore = gse31210.risk.score
                             , dat = t(gse31210.t.exp[intersect(lst.modl$Genes, row.names(gse31210.t.exp)),])
                             , time = gse31210.t.cli.os$OS.time/365
                             , event = gse31210.t.cli.os$OS
                             , cutoff = gse31210.cutoff
                             , title='AMrs'
                             , mks = c(1,3,5))
fig8c=fig8C[[2]]
fig8c

gse31210.subtype.AMrs=ifelse(gse31210.risk.score>gse31210.cutoff,'High','Low')
gse31210.subtype.AMrs=data.frame(gse31210.subtype.AMrs,stringsAsFactors = F)
colnames(gse31210.subtype.AMrs)='Cluster'
table(gse31210.subtype.AMrs)
gse31210.subtype.AMrs$Cluster=factor(gse31210.subtype.AMrs$Cluster,levels = c('Low','High'))

######### GSE50081
match(lst.modl$Genes,row.names(gse50081.t.exp))
lst.modl$Genes

rownames(gse50081.t.exp)[grep("DOWN16",rownames(gse50081.t.exp))]

gse50081.model.dat=gse50081.t.exp[intersect(lst.modl$Genes,row.names(gse50081.t.exp)),]

lst.vd.mod2=createCoxModel((t(gse50081.model.dat))
                           ,time=gse50081.t.cli.os$OS.time/365
                           ,event = gse50081.t.cli.os$OS)
# gse50081.risk.score=lst.vd.mod2$Score
gse50081.risk.score=predictScoreByCoxModel(coxModel = lst.modl,t(gse50081.model.dat))
# gse50081.risk.score=mosaic::zscore(gse50081.risk.score)
lst.modl$fmla
lst.vd.mod2$fmla

gse50081.cutoff <- survminer::surv_cutpoint(data.frame(time=gse50081.t.cli.os$OS.time/365,
                                                       event=gse50081.t.cli.os$OS,
                                                       risk=gse50081.risk.score), time = "time", event = "event",
                                            variables = c("risk"))
gse50081.cutoff=gse50081.cutoff$cutpoint$cutpoint
# gse50081.cutoff=0

fig8D=plotCoxModel_Batch_use(riskScore = gse50081.risk.score
                             , dat = t(gse50081.t.exp[intersect(lst.modl$Genes, row.names(gse50081.t.exp)),])
                             , time = gse50081.t.cli.os$OS.time/365
                             , event = gse50081.t.cli.os$OS
                             , cutoff = gse50081.cutoff
                             , title='AMrs'
                             , mks = c(1, 3, 5))
fig8d=fig8D[[2]]
fig8d

gse50081.subtype.AMrs=ifelse(gse50081.risk.score>gse50081.cutoff,'High','Low')
gse50081.subtype.AMrs=data.frame(gse50081.subtype.AMrs)
colnames(gse50081.subtype.AMrs)='Cluster'
table(gse50081.subtype.AMrs)

####
fig8bc=mg_merge_plot(fig8b,fig8c,fig8d,nrow = 1,ncol = 3,labels = c('D','E','F'))
fig8=mg_merge_plot(fig8a,fig8bc,nrow = 2,ncol = 1,heights = c(1.8,2))
savePDF('PDFs/Fig8.pdf',fig8,width = 16,height = 20)

#### AMrs分组与临床信息的比较分析
tcga.t.cli_use=cbind(luad.tcga.t.exp.os,luad.tcga.t.exp.cli.tnm)
colnames(tcga.t.cli_use)
median(tcga.t.cli_use$Age,na.rm = T)

tcga.t.cli_use$Age1=ifelse(tcga.t.cli_use$Age>60,'>60','<=60')
tcga.t.cli_use$Status=ifelse(tcga.t.cli_use$OS==0,'Alive','Dead')

tcga.t.cli_use=dplyr::rename(tcga.t.cli_use,c('T Stage'='Clinical_T'
                                              ,'N Stage'='Clinical_N'
                                              ,'M Stage'='Clinical_M'
                                              ,'Stage'='Clinical_Stage'))
tcga.t.cli_use$Cluster=tcga.subtype.AMrs[rownames(tcga.t.cli_use),1]
tcga.t.cli_use$AMES=tcga.group[rownames(tcga.t.cli_use),'AMES']
tcga.t.cli_use$AMES=factor(tcga.t.cli_use$AMES,levels = c('Low','Moderate','High'))

head(tcga.t.cli_use)
tcga.t.cli_use$Stage=gsub("Stage ","",tcga.t.cli_use$Stage)
table(tcga.t.cli_use$Stage)
###########
colnames(tcga.t.cli_use)
colnames(tcga.t.cli_use)[c(3:6,9,8,12)]

tcga.risk.cli.cmp=list()
for(i in c(3:6,9,8,12)){
  p=mg_violin_use(data.frame(tcga.t.cli_use[names(tcga.risk.score),i],tcga.risk.score)
                  ,test_method = 'wilcox.test'
                  ,cmp_test_method = 'wilcox.test'
                  # ,group.col = pal_aaas()(9)[c(2,1)]
                  , xlab=colnames(tcga.t.cli_use)[i]
                  , ylab = 'AMrs', melt = T, legend.pos = 't')
  tcga.risk.cli.cmp=c(tcga.risk.cli.cmp,list(p))
}
length(tcga.risk.cli.cmp)

################# 桑基图
colnames(tcga.t.cli_use)
plot_sankey=function(df_m){
  library(ggalluvial)
  library(ggplot2)
  library(dplyr)
  corLodes=to_lodes_form(df_m, axes = 1:ncol(df_m), id = "Cohort")
  mycol <- rep(c("#029149","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),15)
  ggplot(corLodes, aes(x = x, stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) +
    scale_x_discrete(expand = c(0, 0)) +  
    #用aes.flow控制线调颜色，forward说明颜色和前面一致，backward说明与后面一致
    geom_flow(width = 2/10,aes.flow = "forward") + 
    geom_stratum(alpha = .9,width = 2/10) +
    scale_fill_manual(values = mycol) +
    #size = 2代表基因名字大小
    geom_text(stat = "stratum", size = 2,color="black") +
    xlab("") + ylab("") + theme_bw() + 
    theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank()) + #ȥ????????
    theme(panel.grid =element_blank()) + 
    theme(panel.border = element_blank()) + 
    ggtitle("") + guides(fill = FALSE)   
}

fig9b=plot_sankey(tcga.t.cli_use[,c('Cluster','AMES')])

fig9a=mg_merge_plot(c(tcga.risk.cli.cmp,list(fig9b)),nrow = 2,ncol = 4
                    ,labels = LETTERS[1])

################# TCGA
tcga.t.cli_forKM=cbind(tcga.t.cli_use,AMrs=tcga.subtype.AMrs$Cluster)
tcga.t.cli_forKM=data.frame(tcga.t.cli_forKM)

tcga.t.cli_forKM$Age1=factor(tcga.t.cli_forKM$Age1,levels = c('<=60','>60'))

tcga.t.cli_forKM$Gender=as.factor(tcga.t.cli_forKM$Gender)
tcga.t.cli_forKM$Status=ifelse(tcga.t.cli_forKM$OS==1,'Dead','Alive')

tcga.t.cli_forKM$Stage=ifelse(!is.na(tcga.t.cli_forKM$Stage),ifelse(tcga.t.cli_forKM$Stage %in% c('I','II'),'I+II','III+IV'),NA)
# tcga.t.cli_forKM$Grade=ifelse(!is.na(tcga.t.cli_forKM$Grade),ifelse(tcga.t.cli_forKM$Grade %in% c('G1','G2'),'G1+G2','G3+G4'),NA)

tcga.t.cli_forKM$AMrs

colnames(tcga.t.cli_forKM)
colnames(tcga.t.cli_forKM)[c(3:6,9,8,10,11)]
colnames(tcga.t.cli_forKM)[c(6,9,8)]

tcga.subtype.AMrs.km.cmp=list()
for(i in c(6,9,8)){
  df_cli=tcga.t.cli_forKM[,c(colnames(tcga.t.cli_forKM)[i],'OS.time','OS','AMrs')]
  print(names(table(df_cli[,1])))
  if(!(colnames(tcga.t.cli_forKM)[i] %in% c('Viral.etiology','Fibrosis'))){
    for(variate in na.omit(names(table(df_cli[,1])))){
      df_dat=df_cli[which(df_cli[,1]==variate),]
      p=ggplotKMCox(data.frame(time = df_dat$OS.time/365
                               , event = df_dat$OS
                               , df_dat$AMrs)
                    , pal = NULL
                    , title = paste(colnames(tcga.t.cli_forKM)[i],variate,sep = "="))
      tcga.subtype.AMrs.km.cmp=c(tcga.subtype.AMrs.km.cmp,list(p))
      
    }
    
  }else{
    for(variate in na.omit(names(table(df_cli[,1])))){
      if(variate=='Negative'){
        next
      }
      df_dat=df_cli[which(df_cli[,1]==variate),]
      p=ggplotKMCox(data.frame(time = df_dat$OS.time/365
                               , event = df_dat$OS
                               , df_dat$NPRS.Score)
                    ,pal = NULL
                    ,title = variate)
      tcga.subtype.AMrs.km.cmp=c(tcga.subtype.AMrs.km.cmp,list(p))
    }
  }
}

length(tcga.subtype.AMrs.km.cmp)

fig9c=mg_merge_plot(tcga.subtype.AMrs.km.cmp,nrow = 2,ncol = 3,labels = 'B')

fig9=mg_merge_plot(fig9a,fig9c,nrow = 2,ncol = 1,heights = c(2,2))

savePDF(filename = 'PDFs/Fig9.pdf',fig9,height = 20,width = 16)

##### NPRS分组中的突变特征
tcga.subtype.AMrs_use=tcga.subtype.AMrs
tcga.subtype.AMrs_use$Tumor_Sample_Barcode=substr(rownames(tcga.subtype.AMrs_use),1,12)
rownames(tcga.subtype.AMrs_use)=tcga.subtype.AMrs_use$Tumor_Sample_Barcode
tcga.subtype.AMrs_use=tcga.subtype.AMrs_use[,c(2,1)]
table(tcga.subtype.AMrs_use$Cluster)

## 在线绘制
####### 在线绘制
tcga.immu.lands.p1=readMatrix(paste0(MG_Grobal_baseFolder,'/source/PMC5982584_supplement_2.txt'))
tcga.immu.lands.p1<-tcga.immu.lands.p1[tcga.immu.lands.p1$`TCGA Study`=='LUAD',]
rownames(tcga.immu.lands.p1)=paste0(rownames(tcga.immu.lands.p1),'-01')

tcga.subtype.AMrs.forMut=tcga.subtype.AMrs
rownames(tcga.subtype.AMrs.forMut)=substr(rownames(tcga.subtype.AMrs.forMut),1,12)
writeMatrix(tcga.subtype.AMrs.forMut,'files/tcga.subtype.AMrs.forMut.txt')

####
tcga.subtype.AMrs_forAlt=tcga.subtype.AMrs
rownames(tcga.subtype.AMrs_forAlt)=substr(rownames(tcga.subtype.AMrs_forAlt),1,15)

table(is.na(match(row.names(tcga.subtype.AMrs_forAlt),row.names(tcga.immu.lands.p1))))

tcga.immu.lands.p1=tcga.immu.lands.p1[intersect(row.names(tcga.subtype.AMrs_forAlt),row.names(tcga.immu.lands.p1)),]
dim(tcga.immu.lands.p1)

table(tcga.immu.lands.p1$`TCGA Subtype`,tcga.subtype.AMrs_forAlt[rownames(tcga.immu.lands.p1),1])
plotMutiBar_use(table(tcga.immu.lands.p1$`TCGA Subtype`
                      ,tcga.subtype.AMrs_forAlt[rownames(tcga.immu.lands.p1),1])[,])

table(tcga.immu.lands.p1$`Immune Subtype`,tcga.subtype.AMrs_forAlt[rownames(tcga.immu.lands.p1),1])
plotMutiBar_use(table(tcga.immu.lands.p1$`Immune Subtype`
                      ,tcga.subtype.AMrs_forAlt[rownames(tcga.immu.lands.p1),1]))

colnames(tcga.immu.lands.p1)

col.selected=c('Aneuploidy Score','Homologous Recombination Defects','Fraction Altered','Number of Segments','Nonsilent Mutation Rate')
length(col.selected)

pal_aaas()(10)[c(2,3,1,4:9)]

tcga.AMrs.geneAlt.p.all=list()
tcga.AMrs.geneAlt.cor.p.all=list()

for(val in col.selected){
  if(val=='Nonsilent Mutation Rate'){
    p1=mg_violin_use(data.frame(tcga.subtype.AMrs_forAlt[rownames(tcga.immu.lands.p1),1]
                                ,tcga.immu.lands.p1[,val])
                     ,melt = T
                     ,ylab = 'Tumor mutation burden'
                     ,test_method = 'wilcox.test'
                     ,cmp_test_method = 'wilcox.test'
                     ,legend.pos = 'tr'
                     ,jitter=T
                     ,group.col = pal_aaas()(9)[c(2,1)]
                     ,show_compare = T)
  }else{
    p1=mg_violin_use(data.frame(tcga.subtype.AMrs_forAlt[rownames(tcga.immu.lands.p1),1]
                                ,tcga.immu.lands.p1[,val])
                     ,melt = T
                     ,ylab = val
                     ,jitter=T
                     ,group.col = pal_aaas()(9)[c(2,1)]
                     ,test_method = 'wilcox.test'
                     ,cmp_test_method = 'wilcox.test'
                     ,legend.pos = 'tr'
                     ,show_compare = T)
  }
  
  p2=mg_cor_point(tcga.risk.score,tcga.immu.lands.p1[substr(names(tcga.risk.score),1,12),val]
                  , xlab = 'NPRS'
                  , ylab = val
                  # , top_col = 'red'
                  # , right_col = 'brown'
                  , marginal.type = 'density')
  tcga.AMrs.geneAlt.p.all=c(tcga.AMrs.geneAlt.p.all,list(p1))
  tcga.AMrs.geneAlt.cor.p.all=c(tcga.AMrs.geneAlt.cor.p.all,list(p2))
}


# fig8a=mg_merge_plot(tcga.AMrs.geneAlt.p.all,nrow = 1,ncol = 5,common.legend = T)
# fig8b=mg_merge_plot(tcga.AMrs.geneAlt.cor.p.all,nrow = 1,ncol = 5)
# 
# fig8ab=mg_merge_plot(fig8a,fig8b,nrow = 2,ncol = 1,labels = c('A','B'))
# savePDF('PDFs/Fig8AB.pdf',fig8ab,height = 8,width = 16)

####################### KEGG通路评分与风险评分的相关性分析
# library(GSVA)
# library(GSVAdata)
# ########### hallmark
# geneSets <- getGmt("./origin_datas/h.all.v7.5.1.symbols.gmt")
# 
# tcga.exp.h.all <- gsva(as.matrix(luad.tcga.exp), geneSets,min.sz=10,method='ssgsea')
# tcga.exp.h.all <- data.frame(tcga.exp.h.all,check.names = F)
# save(tcga.exp.h.all,file = 'origin_datas/GSVA/tcga.exp.h.all.RData')
# ######### KEGG
# gmtFile='origin_datas/c2.cp.kegg.v7.5.1.symbols.gmt'
# c2KEGG <- getGmt(gmtFile,
#                  collectionType=BroadCollection(category="c2"),
#                  geneIdType=SymbolIdentifier())
# 
# tcga.kegg.ssgsea <- gsva(as.matrix(luad.tcga.exp), 
#                          c2KEGG,
#                          method = 'ssgsea',
#                          min.sz = 10,
#                          max.sz = 500,
#                          verbose = TRUE)
# save(tcga.kegg.ssgsea,file = 'origin_datas/GSVA/tcga.kegg.ssgsea.RData')
load('origin_datas/GSVA/tcga.exp.h.all.RData')
tcga.exp.h.all=tcga.exp.h.all[,names(tcga.risk.score)]
dim(tcga.exp.h.all)

sort(-abs(cor(t(tcga.exp.h.all),(as.numeric(tcga.risk.score)))))[1:32]

write.table(tcga.exp.h.all,'files/文件/tcga.exp.h.all.txt',sep = '\t',quote = F)
write.table(cor(t(tcga.exp.h.all),tcga.risk.score),'files/文件/tcga.h.all.ssgsea.AMrs.cor.txt',sep = '\t',quote = F)
sum(abs(cor(t(tcga.exp.h.all),(as.numeric(tcga.risk.score))))>=0.4)

tcga.AMrs.pathway <- rownames(tcga.exp.h.all[order(abs(cor(t(tcga.exp.h.all),(tcga.risk.score))),
                                                     decreasing=TRUE)[1:12],order(tcga.risk.score)])
length(tcga.AMrs.pathway)

tcga.AMrs.pathway_datas <- t(tcga.exp.h.all[tcga.AMrs.pathway, ])
tcga.AMrs.pathway_datas <- crbind2DataFrame(tcga.AMrs.pathway_datas)
tcga.AMrs.pathway_datas$AMrs <- tcga.risk.score

# pdf('PDFs/Fig10D.pdf',width = 9,height = 9,onefile = F)
# mg_quick_cor_plot(dat1 = tcga.AMrs.pathway_datas,
#                   dat2 = NULL,
#                   cluster_col= T,
#                   R_range=c(-1,0,1),
#                   # P_range=c(5,6)
# )
# dev.off()

pdf('PDFs/Fig10D.pdf',width = 10,height = 10,onefile = F)
quickcor(tcga.AMrs.pathway_datas, cor.test = TRUE) +
  geom_square(data = get_data(type = "lower", show.diag = FALSE)) +
  scale_fill_gradient2n(colours = rev(red_blue())) +
  geom_mark(data = get_data(type = "upper", show.diag = FALSE),size=3,sep = "\n") +
  geom_abline(slope = -1, intercept = 14) + remove_x_axis()
dev.off()

############# 差异分析
all(rownames(tcga.subtype.AMrs)==colnames(tcga.exp.h.all))

tcga.exp.h.all_use=t(tcga.exp.h.all)
tcga.exp.h.all_use=data.frame(tcga.exp.h.all_use,check.names = F)
tcga.exp.h.all_use[1:4,1:5]
tcga.exp.h.all_use$sample=rownames(tcga.exp.h.all_use)
head(tcga.exp.h.all_use)

tcga.exp.h.all_long=reshape2::melt(tcga.exp.h.all_use,id.vars='sample')
head(tcga.exp.h.all_long)
colnames(tcga.exp.h.all_long)=c('sample','pathway','ssgsea')
head(tcga.exp.h.all_long)

##############
tcga.subtype.AMrs_use=tcga.subtype.AMrs
tcga.subtype.AMrs_use$sample=rownames(tcga.subtype.AMrs_use)
###########
tcga.exp.h.all_long=merge(tcga.exp.h.all_long,tcga.subtype.AMrs_use,by='sample')
head(tcga.exp.h.all_long)

########### 差异比较
tcga.h.all.AMrs.cmp=compare_means(ssgsea ~ Cluster,group.by='pathway',data=tcga.exp.h.all_long)
tcga.h.all.AMrs.cmp=data.frame(tcga.h.all.AMrs.cmp)
rownames(tcga.h.all.AMrs.cmp)=tcga.h.all.AMrs.cmp$pathway
tcga.h.all.AMrs.cmp=tcga.h.all.AMrs.cmp[order(tcga.h.all.AMrs.cmp$p.adj,decreasing = F),]
head(tcga.h.all.AMrs.cmp)

table(tcga.h.all.AMrs.cmp$p.adj<0.0001)
match(tcga.AMrs.pathway,tcga.h.all.AMrs.cmp$pathway)

################################### 山峦图
tcga.h.all.selected=tcga.h.all.AMrs.cmp$pathway[1:10]

head(tcga.subtype.AMrs)
table(tcga.subtype.AMrs$Cluster)

all(colnames(tcga.exp.h.all)==rownames(tcga.subtype.AMrs))

#############
head(tcga.h.all.AMrs.cmp)
tcga.exp.h.all_forCmp=tcga.exp.h.all
rownames(tcga.exp.h.all_forCmp)=paste(tcga.h.all.AMrs.cmp[rownames(tcga.exp.h.all_forCmp),'p.signif']
                                      ,rownames(tcga.exp.h.all_forCmp)
                                      ,sep=" ")
rownames(tcga.exp.h.all_forCmp)

p1=mg_ridges_plot(data=t(tcga.exp.h.all_forCmp[tcga.h.all.selected,which(tcga.subtype.AMrs$Cluster=='Low')]),xlab = "ssGSEA score",col=mg_colors[1:12])+
  theme_bw()+theme(legend.position = 'none')+
  labs(title = "AMrs-low")

p2=mg_ridges_plot(data=t(tcga.exp.h.all_forCmp[tcga.h.all.selected,which(tcga.subtype.AMrs$Cluster=='High')]),xlab = "ssGSEA score",col=mg_colors[1:12])+
  theme_bw()+theme(legend.position = 'none')+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank() )+
  labs(title = "AMrs-high")


fig10=cowplot::plot_grid(p1,p2,align = 'h',
                         nrow = 1,ncol = 2,rel_widths = c(2.5,1))
fig10
savePDF('PDFs/Fig10C.pdf',fig10,height = 6,width = 8)

get_PlotMutiBoxplot(t(tcga.exp.h.all_forCmp)[rownames(tcga.subtype.AMrs),tcga.h.all.selected]
                    ,tcga.subtype.AMrs
                    ,group.val='Cluster'
                    ,xangle=30
                    ,group_cols=ggsci::pal_lancet()(9)[c(1,2)])

######### 免疫特征
fig10a=get_PlotMutiBoxplot(tcga.exp.cibersort[rownames(tcga.subtype.AMrs),1:22]
                           ,tcga.subtype.AMrs
                           ,group.val='Cluster'
                           ,xangle=30
                           ,group_cols=ggsci::pal_lancet()(9)[c(1,2)])
# fig10b=get_PlotMutiBoxplot(gse72094.t.exp.cibersort[rownames(gse72094.subtype.AMrs),1:22]
#                            ,gse72094.subtype.AMrs
#                            ,xangle=30
#                            ,group_cols=ggsci::pal_lancet()(9)[c(1,2)])

fig10c=get_PlotMutiBoxplot(tcga.exp.estimate[rownames(tcga.subtype.AMrs),c(1:3)]
                           ,tcga.subtype.AMrs
                           ,group.val='Cluster'
                           ,xangle=0
                           ,group_cols=ggsci::pal_lancet()(9)[c(1,2)])
# fig10d=get_PlotMutiBoxplot(gse72094.t.exp.estimate[rownames(gse72094.subtype.AMrs),1:3]
#                            ,gse72094.subtype.AMrs
#                            ,xangle=0
#                            ,group_cols=ggsci::pal_lancet()(9)[c(1,2)])

dim(tcga.exp.cibersort)

(apply(tcga.exp.cibersort[,1:22], 2, function(x){sd(x)}))

x_matrix=cbind(tcga.exp.cibersort[names(tcga.risk.score),c(1:4,6:22)],AMrs=tcga.risk.score)
fig10e=mg_quick_cor_plot(dat1 = x_matrix
                         ,R_range=c(-1,0,1))

# pdf('PDFs/Fig10C.pdf',width = 10,height = 10)
# quickcor(x_matrix, cor.test = TRUE) +
#   geom_square(data = get_data(type = "lower", show.diag = FALSE)) +
#   scale_fill_gradient2n(colours = rev(red_blue())) +
#   geom_mark(data = get_data(type = "upper", show.diag = FALSE),size=3,sep = "\n",sig.thres=0.05) +
#   geom_abline(slope = -1, intercept = 23) + remove_x_axis()
# dev.off()
# fig10ad=mg_merge_plot(fig10a,fig10c,fig10b,fig10d
#                       ,nrow = 2,ncol = 2
#                       ,widths = c(0.6,0.35),labels = LETTERS[1:4]
#                       ,common.legend = T)
# fig10ad

savePDF('PDFs/Fig10A.pdf',fig10a,width = 8,height = 6)
pdf('PDFs/Fig10B.pdf',width = 8,height = 8)
fig10e$plot
dev.off()

######################### 免疫治疗响应预测
################# 评分
##########################
inflam.genes=read.table('origin_datas/T_cell-inflamed.txt',header = T)
inflam.genes=inflam.genes$Gene
setdiff(inflam.genes,rownames(tcga.t.exp))

dim(tcga.t.exp.sub)
dim(tcga.group)
all(rownames(tcga.group)==colnames(tcga.t.exp.sub))
range(scale(t(tcga.t.exp.sub)),na.rm = T)

table(tcga.group$group)
table(tcga.group$AMES)

tcga.smp.low.high=rownames(tcga.group)[tcga.group$group %in% c('<1','2~3')]

col_an=HeatmapAnnotation(type = tcga.group[,'AMES'], 
                         show_annotation_name = F, 
                         col = list(type = c("Low" = "#008B45FF","Moderate"="#3B4992FF","High"="#EE0000FF")), 
                         show_legend = T,  
                         annotation_legend_param = list(title = "Type",ncol = 1), 
                         which = "col" )

Heatmap(as.matrix(t(scale(t(tcga.t.exp.sub[inflam.genes,]))))
        , name = "Expr"
        # , column_km =3
        , cluster_rows = T
        , row_title_gp = gpar(fill = pal_lancet('lanonc',alpha =0.6)(9)[c(7,1)])
        , show_row_dend = F
        , column_split = tcga.subtype.AMrs[,'Cluster']
        , cluster_columns = T
        , cluster_column_slices=T
        , show_column_dend = T
        , show_column_names = F
        , col = circlize::colorRamp2(c(-4, 0,4), c('#3B4992FF','white','#EE0000FF'))
        , column_title_gp = gpar(fill = (pal_aaas()(9)[c(3,2)]))
        # , top_annotation = col_an
        , border = TRUE)
#######
tcga.inflam.ssgsea=ssGSEAScore_by_genes(gene.exp = tcga.t.exp.sub,genes = as.character(inflam.genes))
tcga.inflam.ssgsea=data.frame(t(tcga.inflam.ssgsea),check.names = F)
head(tcga.inflam.ssgsea)

fig11a=mg_PlotMutiBoxplot(data = data.frame(AMrs=tcga.inflam.ssgsea[rownames(tcga.subtype.AMrs),1])
                   , group = tcga.subtype.AMrs$Cluster
                   , legend.pos = 'tr'
                   , add = 'boxplot'
                   , ylab = 'T‐cell–inflamed GEP score'
                   , xangle=0
                   , group_cols = pal_lancet()(9)[c(1, 2)]
                   # , test_method = 't.test'
                   )

cor.test(x=tcga.inflam.ssgsea$GeneSet,tcga.risk.score,method = 'spearman')
#################
tcga.IFNγScore=immu_Th1_IFNγScore(tcga.t.exp.sub)
tcga.IFNγScore=data.frame(tcga.IFNγScore,check.names = F)
head(tcga.IFNγScore)

fig11b=mg_PlotMutiBoxplot(data = data.frame(AMrs=(tcga.IFNγScore[,1]))
                   , group = tcga.subtype.AMrs[,'Cluster']
                   , legend.pos = 'tr'
                   , add = 'boxplot'
                   , ylab = 'Th1/IFNγ gene signature ssGSEA score'
                   , xangle=0
                   , group_cols = pal_lancet()(9)[c(1, 2)]
                   , test_method = 'kruskal.test')

####### 免疫检查点抑制
mg_PlotMutiBoxplot(tcga.t.exp.icg[rownames(tcga.subtype.AMrs),]
                   , group = tcga.subtype.AMrs$Cluster
                   , legend.pos = 'right'
                   , add = 'boxplot'
                   , ylab = 'Normalized gene expression'
                   , group_cols = ggsci::pal_lancet()(9)[c(1, 2)]
                   , test_method = 'kruskal.test')

####################
setdiff(c('PDCD1','CD274','CTLA4'),rownames(tcga.t.exp.sub))
p.all=list()
for(gene in c('PDCD1','CD274','CTLA4')){
  p=mg_PlotMutiBoxplot(data = data.frame(AMrs=as.numeric(tcga.t.exp.sub[gene,]))
                       , group = tcga.subtype.AMrs[,'Cluster']
                       , legend.pos = 'tr'
                       , add = 'boxplot'
                       , ylab = 'log2(TPM+1)'
                       , xangle=0
                       , group_cols = pal_lancet()(9)[c(1, 2)]
                       # , test_method = 'kruskal.test'
                       )+labs(title=gene)
    
  p.all=c(p.all,list(p))
}
length(p.all)
fig11c=mg_merge_plot(p.all,nrow = 1,ncol = 3,labels = c('A'),common.legend = T)

##################
icg.genes=c('PDCD1','CD274','CTLA4')
tcga.t.exp.sub_use=data.frame(t(tcga.t.exp.sub))
tcga.t.exp.sub_use$sample=rownames(tcga.t.exp.sub_use)
tcga.t.exp.sub_use[1:4,1:3]

tcga.exp.icg_long=reshape2::melt(tcga.t.exp.sub_use[,c('sample',icg.genes)],id.vars='sample')
head(tcga.exp.icg_long)
colnames(tcga.exp.icg_long)=c('sample','gene','expr')
head(tcga.exp.icg_long)

##############
tcga.exp.icg_long=merge(tcga.exp.icg_long,tcga.subtype.AMrs_use,by='sample')
head(tcga.exp.icg_long)

tcga.icg.cmp=compare_means(expr ~ Cluster,group.by='gene',data=tcga.exp.icg_long)
tcga.icg.cmp=data.frame(tcga.icg.cmp)

rownames(tcga.icg.cmp)=tcga.icg.cmp$gene
tcga.icg.cmp=tcga.icg.cmp[order(tcga.icg.cmp$p.adj,decreasing = F),]
head(tcga.icg.cmp)

tcga.exp.icg_m=reshape2::dcast(tcga.exp.icg_long,sample ~ gene,value.var='expr')
rownames(tcga.exp.icg_m)=tcga.exp.icg_m$sample
tcga.exp.icg_m=tcga.exp.icg_m[,-1]
tcga.exp.icg_m=tcga.exp.icg_m[rownames(tcga.subtype.AMrs),]

dim(tcga.exp.icg_m)
colnames(tcga.exp.icg_m)

colnames(tcga.exp.icg_m)=paste(tcga.icg.cmp[colnames(tcga.exp.icg_m),'p.signif']
                                      ,colnames(tcga.exp.icg_m)
                                      ,sep=" ")

p1=mg_ridges_plot(data=tcga.exp.icg_m[which(tcga.subtype.AMrs$Cluster=='Low'),],xlab = "log2(TPM+1)")+
  theme_bw()+theme(legend.position = 'none')+
  labs(title = "AMrs-low")

p2=mg_ridges_plot(data=tcga.exp.icg_m[which(tcga.subtype.AMrs$Cluster=='High'),],xlab = "log2(TPM+1)")+
  theme_bw()+theme(legend.position = 'none')+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank() )+
  labs(title = "AMrs-high")

fig11c=cowplot::plot_grid(p1,p2,align = 'h',
                         nrow = 1,ncol = 2,rel_widths = c(1,1))
fig11c

####################
tcga.cytscore=immu_cytScore(tcga.t.exp.sub)
tcga.cytscore=data.frame(tcga.cytscore,check.names = F)
head(tcga.cytscore)

fig11d=mg_PlotMutiBoxplot(data = data.frame(AMrs=(tcga.cytscore[,1]))
                   , group = tcga.subtype.AMrs$Cluster
                   , legend.pos = 'tr'
                   , add = 'boxplot'
                   , ylab = 'Cytolytic activity'
                   , xangle=0
                   , group_cols = pal_lancet()(9)[c(1, 2)]
                   , test_method = 'kruskal.test')

################### TIDE 软件 预测 免疫治疗响应
# writeMatrix(t(scale(t(luad.tcga.t.exp),scale = F)),outpath = 'origin_datas/TCGA/luad.tcga.t.exp.zscore.txt')

tcga.tide<-read.csv('origin_datas/TCGA-LUAD_TIDE.csv',row.names = 1,stringsAsFactors = F)
tcga.tide=tcga.tide[rownames(tcga.subtype.AMrs),]
dim(tcga.tide)

########
tide.selected=c('MDSC','CAF','TAM.M2','Exclusion','Dysfunction','TIDE')
####### GSE72094
tcga.AMrs.tide.p.all=c()
for(i in tide.selected){
  p=mg_violin_use(data.frame(tcga.subtype.AMrs$Cluster,tcga.tide[rownames(tcga.subtype.AMrs),i])
                  , test_method = 'wilcox.test'
                  , cmp_test_method='wilcox.test'
                  , group.col = ggsci::pal_lancet()(9)[c(1, 2)]
                  , ylab = i, melt = T, legend.pos = 'tr')
  tcga.AMrs.tide.p.all=c(tcga.AMrs.tide.p.all,list(p))
}

fig11e=mg_merge_plot(tcga.AMrs.tide.p.all,nrow = 1,ncol = 6,common.legend=T)

fig11e=mg_PlotMutiBoxplot(data = data.frame(AMrs=(tcga.tide[rownames(tcga.subtype.AMrs),'TIDE']))
                          , group = tcga.subtype.AMrs$Cluster
                          , legend.pos = 'tr'
                          , add = 'boxplot'
                          , ylab = 'TIDE'
                          , xangle=0
                          , group_cols = pal_lancet()(9)[c(1, 2)]
                          , test_method = 'kruskal.test')

head(tcga.mutload_use)

fig11e=mg_PlotMutiBoxplot(data = data.frame(AMrs=(tcga.mutload_use[substr(rownames(tcga.subtype.AMrs),1,12),'total_perMB']))
                   , group = tcga.subtype.AMrs$Cluster
                   , legend.pos = 'tr'
                   , add = 'boxplot'
                   , ylab = 'TMB'
                   , xangle=0
                   , group_cols = pal_lancet()(9)[c(1, 2)]
                   , test_method = 'kruskal.test')

##################### 药物IC50
############# GDSC数据库细胞系表达谱
GDSC_exp<-read.table('origin_datas/GDSC/Cell_line_RMA_proc_basalExp/Cell_line_RMA_proc_basalExp.txt',header = T,check.names = F,sep = '\t')
GDSC_exp[1:4,1:5]
setdiff(lst.modl$Genes,GDSC_exp$GENE_SYMBOLS)
table(duplicated(GDSC_exp$GENE_SYMBOLS))
######## 去掉重复的GENE_SYMBOLS
GDSC_exp<-GDSC_exp[!duplicated(GDSC_exp$GENE_SYMBOLS),]
rownames(GDSC_exp)<-GDSC_exp$GENE_SYMBOLS
GDSC_exp$GENE_title<-NULL
GDSC_exp$GENE_SYMBOLS<-NULL
GDSC_exp[1:4,1:5]
range(GDSC_exp)
dim(GDSC_exp)

setdiff(lst.modl$Genes,row.names(GDSC_exp))
grep("TLSP",row.names(GDSC_exp))

GDSC.model.dat=GDSC_exp[intersect(lst.modl$Genes,row.names(GDSC_exp)),]
dim(GDSC.model.dat)

############## 风险模型
lst.modl.coef=lst.modl$Coef
names(lst.modl.coef)=lst.modl$Genes

rownames(GDSC.model.dat)

GDSC.risk.score=as.matrix(t(GDSC.model.dat)) %*% as.numeric(lst.modl.coef[rownames(GDSC.model.dat)])
GDSC.risk.score=mosaic::zscore(GDSC.risk.score)
GDSC.risk.score=data.frame(AMrs=GDSC.risk.score,sample=row.names(GDSC.risk.score))
head(GDSC.risk.score)
rownames(GDSC.risk.score)<-gsub('DATA.','',rownames(GDSC.risk.score))
GDSC.risk.score$sample<-rownames(GDSC.risk.score)
head(GDSC.risk.score)

############# 药物IC50
gdsc_drug_ic<-read.csv('origin_datas/GDSC/PANCANCER_IC_Fri Jun 10 02_20_40 2022.csv',header = T)
dim(gdsc_drug_ic)
head(gdsc_drug_ic)
colnames(gdsc_drug_ic)

gdsc_drug_ic[,10:13]<-NULL
gdsc_drug_ic[,5:8]<-NULL
# colnames(gdsc_drug_ic)
# gdsc_drug_ic=gdsc_drug_ic[,c('Drug.name','Drug.Id','Cell.line.name','Cosmic.sample.Id','AUC')]
colnames(gdsc_drug_ic)[4]<-"sample"
head(gdsc_drug_ic)

#################################
GDSC_AMrs_drug<-merge(gdsc_drug_ic,GDSC.risk.score,by='sample')
head(GDSC_AMrs_drug)

drug_name<-as.data.frame(table(GDSC_AMrs_drug$Drug.name))
colnames(drug_name)<-c("drug","counts")

library(psych)
GDSC_AMrs_drug.cor=GDSC_AMrs_drug %>%
  group_by(Drug.name) %>%
  summarise(estimate.rho=corr.test(AUC,AMrs,method = "spearman",adjust = "fdr")$r
            ,pvalue=corr.test(AUC,AMrs,method = "spearman",adjust = "fdr")$p
            ,FDR=corr.test(AUC,AMrs,method = "spearman",adjust = "fdr")$p.adj)
GDSC_AMrs_drug.cor=as.data.frame(GDSC_AMrs_drug.cor)
head(GDSC_AMrs_drug.cor)

########
estimate.rho.cutoff=0.3
GDSC_AMrs_drug.cor$Class=case_when(GDSC_AMrs_drug.cor$estimate.rho>0.3 & GDSC_AMrs_drug.cor$pvalue<0.05 ~ "positive",
                                   GDSC_AMrs_drug.cor$estimate.rho < (-0.3) & GDSC_AMrs_drug.cor$pvalue < 0.05 ~ "negative",
                                   TRUE ~ "ns" )
head(GDSC_AMrs_drug.cor)

drug_AMrs<-GDSC_AMrs_drug.cor[which(GDSC_AMrs_drug.cor$Class != 'ns'),]
drug_AMrs<-drug_AMrs[order(drug_AMrs$estimate.rho),]
drug_AMrs$`-log10(pvalue)`<-(-log10(drug_AMrs$pvalue))

write.table(GDSC_AMrs_drug.cor,'files/GDSC_AMrs.corr.txt',row.names=T,col.names=T,quote=F,sep="\t")

range(drug_AMrs$`-log10(pvalue)`)
range(drug_AMrs$estimate.rho)
library(tidyverse)
fig11f=drug_AMrs %>% 
  ggplot(aes(reorder(Drug.name, estimate.rho), estimate.rho)) + 
  geom_col(aes(fill = `-log10(pvalue)`)) + 
  scale_fill_gradient2(low = "blue", 
                       mid = 'white',
                       high = "red", 
                       midpoint = 1.3) + 
  #coord_flip() + 
  labs(x = "")+
  labs(y = "Rs of drug sensitivity and AMrs")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

########################### 药物靶向的通路
drug_pathway<-read.csv('origin_datas/GDSC/PANCANCER_ANOVA_Fri Jun 10 02_18_35 2022.csv',header = T,stringsAsFactors =F)
dim(drug_pathway)
colnames(drug_pathway)
drug_pathway[,c(5:22)]<-NULL
drug_pathway[,c(2:3)]<-NULL
head(drug_pathway)
drug_pathway$target_pathway

drug_pathway<-drug_pathway[which(drug_pathway$drug_name %in% drug_AMrs$Drug.name),]
head(drug_pathway)
drug_pathway_use<-as.data.frame(table(drug_pathway))
head(drug_pathway_use)
drug_pathway_use<-drug_pathway_use[which(drug_pathway_use$Freq != 0),]
colnames(drug_pathway_use)<-c("Drugs","pathway","freq")
head(drug_pathway_use)
dim(drug_pathway_use) ### 8个通路

drug_pathway_use_final<-merge(drug_pathway_use,drug_AMrs,by.x='Drugs',by.y='Drug.name')
head(drug_pathway_use_final)

drug_pathway_use_final[,3]<-NULL
drug_pathway_use_final[,6:8]<-NULL
head(drug_pathway_use_final)
length(drug_pathway_use_final$pathway)

drug_pathway_use_final=drug_pathway_use_final[which(drug_pathway_use_final$estimate.rho<0),]

pathway<-as.data.frame(table(drug_pathway_use_final$pathway))
pathway_use<-as.character(pathway$Var1)
pathway_use
length(pathway_use)

########## 相关性
head(drug_pathway_use_final)
down=reshape2::dcast(drug_pathway_use_final,Drugs ~ pathway,value.var='estimate.rho',fill=0)
rownames(down)=down$Drugs
down=down[,-1]
down=t(down)
down=as.matrix(down)
######### 相关性显著性p值
up=reshape2::dcast(drug_pathway_use_final,Drugs ~ pathway,value.var='pvalue',fill=0)
rownames(up)=up$Drugs
up=up[,-1]
up=t(up)
up=(-log10(up))
up[which(is.infinite(up),arr.ind = )]=0
up=as.matrix(up)

############ 可视化
library("ComplexHeatmap")
library("circlize")
plotMutiHeatmap=function(up,down,up_break,up_colors,down_break,down_colors,title){
  UpColor <- colorRamp2(breaks = up_break, colors = up_colors)
  DnColor <- colorRamp2(breaks = down_break, colors = down_colors)
  
  DiagFunc <- function(up, down){
    function(j, i, x, y, width, height, fill){
      grid.polygon(unit.c(x - 0.5*width, x - 0.5*width, x + 0.5*width), 
                   unit.c(y - 0.5*height, y + 0.5*height, y - 0.5*height),
                   gp = gpar(fill = DnColor(down[i, j]), col = "grey")) 
      
      grid.polygon(unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width), 
                   unit.c(y + 0.5*height, y - 0.5*height, y + 0.5*height),
                   gp = gpar(fill = UpColor(up[i, j]), col = "grey"))
      ##### 显著性
      if(up[i, j]>=1.3){
        txt="***"
        if(up[i, j]>=1.3&up[i, j]<2){
          txt='*'
        }else if(up[i, j]>=2&up[i, j]<3){
          txt='**'
        }else if(up[i, j]>=3&up[i, j]<4){
          txt='***'
        }
        grid.text(label=txt,x=(x + 0.5*width),
                  y=(y+ 0.5*height),just = c('right','top'))
      }
      if(down[i, j]>0.1){
      }
    }
  }
  ######## 热图
  p1 <- Heatmap(up, column_title = title
                , rect_gp = gpar(type = "none")
                , show_heatmap_legend = F
                , cluster_rows = T
                , cluster_columns = T
                , cell_fun = DiagFunc(up = up, down = down) 
  ) 
  ############## legend correlation
  col_fun = colorRamp2(down_break, down_colors) 
  lgd <- Legend(title = "Correlation", 
                col_fun = col_fun, 
                # at = c(-0.3,0,0.3), 
                # labels = c("-0.3","0","0.3"),  
                direction = "horizontal" 
  )
  ############## legend correlation p value
  col_fun2 = colorRamp2(up_break, up_colors) 
  lgd2 <- Legend(title = "-log10(p value)", 
                 col_fun = col_fun2, 
                 at = c(0,1,2,3,4,5), 
                 labels = c('0',"1","2","3","4",">5"),  
                 direction = "horizontal"
  )
  
  draw(p1, annotation_legend_list = list(lgd,lgd2), annotation_legend_side = "bottom"
       ,heatmap_legend_side = "bottom", merge_legend = TRUE)
}

up_break=c(0, 5)
down_break=c(-0.5,0,0.5)
up_colors=c("#FFFFFF","#6f9a8d")
down_colors=c("blue",'white',"red")
range(down)
range(up)

pdf('PDFs/Fig11G.pdf',width = 6,height = 6)
plotMutiHeatmap(up=up,down=down
                ,up_break=up_break
                ,up_colors=up_colors
                ,down_break=down_break
                ,down_colors=down_colors,title='')
dev.off()

fig11ab=mg_merge_plot(fig11e,fig11a,fig11b,nrow = 1,ncol = 3,labels = c('A','B','C'),common.legend = T)
fig11df=mg_merge_plot(fig11c,fig11d,fig11f,nrow = 1,ncol = 3,labels = c('D','E','F'),common.legend = T,widths = c(1.2,1,1.5))

fig11=mg_merge_plot(fig11ab,fig11def
                    ,nrow = 2,ncol = 1
                    ,common.legend = T)

savePDF('PDFs/Fig11.pdf',fig11,height =8,width = 12)

########
mg_compare_uni_muti_cox_use=function(dat,event,os){
  # dat=crbind2DataFrame(dat)
  sig.clini=colnames(dat)
  dat$time=os
  dat$status=event
  dat=dat[which(!is.na(os)&!is.na(event)),]
  all.cox=rbind()
  rnames=c()
  for(s in sig.clini){
    fmla <- as.formula(paste0("Surv(time, status) ~",s))
    cox <- coxph(fmla, data = dat)
    #summary(cox)[[7]]
    #print(summary(cox))
    re=cbind(summary(cox)[[7]][,5],summary(cox)[[7]][,2],summary(cox)[[8]][,3],summary(cox)[[8]][,4])
    if(nrow(re)==1){
      rnames=c(rnames,s)
    }else{
      rnames=c(rnames,row.names(summary(cox)[[7]]))
    }
    all.cox=rbind(all.cox,re)
  }
  row.names(all.cox)=rnames
  colnames(all.cox)=c('p.value','HR','Low 95%CI','High 95%CI')
  
  fmla <- as.formula(paste0("Surv(time, status) ~",paste0(sig.clini,collapse = '+')))
  cox <- coxph(fmla, data = dat)
  muti.re=cbind(summary(cox)[[7]][,5],summary(cox)[[7]][,2],summary(cox)[[8]][,3],summary(cox)[[8]][,4])
  row.names(muti.re)=row.names(summary(cox)[[7]])
  colnames(muti.re)=c('p.value','HR','Low 95%CI','High 95%CI')
  return(list(muti=crbind2DataFrame(muti.re),uni=crbind2DataFrame(all.cox)))
}
getForestplotData=function(res_sig){
  res_sig<-signif(res_sig,digits=2)
  res_sig$CI_for_HR=paste0(" (",res_sig$`Low 95%CI`, "-", res_sig$`High 95%CI`, ")")
  colnames(res_sig)=c("p.value","HR","CI_lower","CI_upper", "(95%_CI_for_HR)")
  
  forest_table<-data.frame(Features=rownames(res_sig),HR=res_sig$HR,`(95%CI)`=res_sig$`(95%_CI_for_HR)`, pvalue=res_sig$p.value,check.names = F,stringsAsFactors = F)
  forest_table$sig<-mg_format_p_values(forest_table$pvalue)
  
  forest_table2<-data.frame(Features="Features",HR="HR",`(95%CI)`="(95%CI)", pvalue="p-value",sig='Significant',check.names = F,stringsAsFactors=F)
  tabletext<-rbind(forest_table2,forest_table)
  
  forest_stastic<-data.frame(mean=as.numeric(as.character(res_sig$HR)),lower=as.numeric(as.character(res_sig$CI_lower)),upper=as.numeric(as.character(res_sig$CI_upper)))
  forest_stastic1<-data.frame(mean=NA,lower=NA,upper=NA)
  cochrane_from_rmeta<-rbind(forest_stastic1,forest_stastic)
  return(list(tabletext,cochrane_from_rmeta))
}

library(forestplot)
all.t.cli.forModel=cbind(tcga.t.cli_use,Group=tcga.subtype.AMrs$Cluster,AMrs=tcga.risk.score)
all.t.cli.forModel=data.frame(all.t.cli.forModel)

all.t.cli.forModel$Stage=ifelse(is.na(all.t.cli.forModel$Stage),NA,ifelse(all.t.cli.forModel$Stage %in% c('I','II'),'I+II','III+IV'))
all.t.cli.forModel$Stage=factor(all.t.cli.forModel$Stage,levels = c('I+II','III+IV'))

all.t.cli.forModel$T.Stage=ifelse(is.na(all.t.cli.forModel$T.Stage),NA,ifelse(all.t.cli.forModel$T.Stage %in% c('T1','T2'),'T1+T2','T3+T4'))
all.t.cli.forModel$N.Stage=ifelse(is.na(all.t.cli.forModel$N.Stage),NA,ifelse(all.t.cli.forModel$N.Stage %in% c('N0'),'N0','NX'))
all.t.cli.forModel$M.Stage=ifelse(is.na(all.t.cli.forModel$M.Stage),NA,ifelse(all.t.cli.forModel$M.Stage %in% c('M0'),'M0','MX'))

table(all.t.cli.forModel$T.Stage)
table(all.t.cli.forModel$N.Stage)
table(all.t.cli.forModel$M.Stage)

all.t.cli.forModel$T.Stage=factor(all.t.cli.forModel$T.Stage,levels = c('T1+T2','T3+T4'))
all.t.cli.forModel$N.Stage=factor(all.t.cli.forModel$N.Stage,levels = c('N0','NX'))
all.t.cli.forModel$M.Stage=factor(all.t.cli.forModel$M.Stage,levels = c('M0','MX'))

all.t.cli.forModel$Age1=factor(all.t.cli.forModel$Age1,levels = c('<=60','>60'))
all.t.cli.forModel$Gender=factor(all.t.cli.forModel$Gender,levels = c('FEMALE','MALE'))
all.t.cli.forModel$Status=ifelse(all.t.cli.forModel$OS==1,'Dead','Alive')

str(all.t.cli.forModel)
colnames(all.t.cli.forModel)

table(all.t.cli.forModel$T.Stage)
table(all.t.cli.forModel$N.Stage)
table(all.t.cli.forModel$M.Stage)
table(all.t.cli.forModel$Stage)

all.t.cli.forCox=all.t.cli.forModel
all.t.cli.forCox=data.frame(all.t.cli.forCox)

colnames(all.t.cli.forCox)
colnames(all.t.cli.forCox)[c(3:6,7,8,14)]

cm.cox=mg_compare_uni_muti_cox_use(all.t.cli.forCox[,c(3:6,7,8,14)]
                                   ,os=all.t.cli.forCox$OS.time
                                   ,event = all.t.cli.forCox$OS)

cm.cox.uni=signif(cm.cox$uni,digits=3)
cm.cox.uni$`Hazard Ratio(95%CI)`=paste0(cm.cox.uni$HR,"(",cm.cox.uni$`Low 95%CI`, "-", cm.cox.uni$`High 95%CI`, ")")

table(cm.cox.uni$p.value<0.05)
cm.cox.uni[cm.cox.uni$p.value<0.05,]

colnames(all.t.cli.forCox)[c(3:6,14)]

cm.cox2=mg_compare_uni_muti_cox_use(all.t.cli.forCox[,c(3:6,14)],os=all.t.cli.forCox$OS.time,event = all.t.cli.forCox$OS)
cm.cox.muti=signif(cm.cox2$muti,digits=3)
cm.cox.muti$`Hazard Ratio(95%CI)`=paste0(cm.cox.muti$HR,"(",cm.cox.muti$`Low 95%CI`, "-", cm.cox.muti$`High 95%CI`, ")")

cm.cox.uni
cm.cox.muti

cm.cox.muti[which(cm.cox.muti$p.value<0.05),]

writeMatrix(cm.cox.uni,outpath = 'files/tcga.cli.single.cox.txt')
writeMatrix(cm.cox.muti,outpath = 'files/tcga.cli.multi.cox.txt')

cm.cox.uni=getForestplotData(cm.cox$uni)
cm.cox.muti=getForestplotData(cm.cox2$muti)

pdf('PDFs/Fig12E.pdf',height = 5,width = 6,onefile = F)
tabletext=cm.cox.uni[[1]]
cochrane_from_rmeta=cm.cox.uni[[2]]
forestplot(tabletext, 
           mean=cochrane_from_rmeta$mean,
           lower=cochrane_from_rmeta$lower,
           upper=cochrane_from_rmeta$upper,
           graph.pos=2,
           hrzl_lines=list('2'=gpar(lty=1,col="black"),
                           '6'=gpar(lty=1,col="black")),
           zero = 1, 
           xlog=FALSE,
           fn.ci_norm = fpDrawDiamondCI,
           boxsize = 0.3, 
           col=fpColors(line = "darkred", box="darkred",zero = 'black',summary='black'), 
           lty.ci = 7,  
           lwd.ci = 3,
           ci.vertices.height = 0.2, 
           lineheight = "auto", 
           xlab="Hazard ratio" )
dev.off()

pdf('PDFs/Fig12F.pdf',height = 5,width = 6,onefile = F)
tabletext=cm.cox.muti[[1]]
cochrane_from_rmeta=cm.cox.muti[[2]]
forestplot(tabletext, 
           mean=cochrane_from_rmeta$mean,
           lower=cochrane_from_rmeta$lower,
           upper=cochrane_from_rmeta$upper,
           graph.pos=2,
           hrzl_lines=list('2'=gpar(lty=1,col="black"),
                           '4'=gpar(lty=1,col="black")),
           zero = 1, 
           xlog=FALSE,
           fn.ci_norm = fpDrawDiamondCI,
           boxsize = 0.3, 
           col=fpColors(line = "darkred", box="darkred",zero = 'black',summary='black'), 
           lty.ci = 7,  
           lwd.ci = 3,
           ci.vertices.height = 0.2, 
           lineheight = "auto", 
           xlab="Hazard ratio" )

dev.off()

mg_plotDCA=function(status,fmlas,modelNames,data){
  set.seed(123456789)
  all.mod=list()
  for(i in 1:length(fmlas)){
    fmla <- as.formula(paste0("status~",fmlas[i]))
    model<-rmda::decision_curve(fmla,
                                data=data,
                                bootstraps=500)
    all.mod=c(all.mod,list(model))
  }
  rmda::plot_decision_curve(all.mod,
                            curve.names=modelNames,
                            # col=mg_colors[c(1,10:12,4,5,7:8)],
                            col=mg_colors[c(1,11,12,8)],
                            xlim=c(0,1),legend.position="topright",
                            lwd=1,
                            confidence.intervals=FALSE)
}

colnames(all.t.cli.forCox)
colnames(all.t.cli.forCox)[c(3:4,14)]

nom.plot=mg_nomogram(all.t.cli.forCox[,c(3:4,14)]
                     ,os = all.t.cli.forCox$OS.time
                     ,status = all.t.cli.forCox$OS)

mg_nomogram_buti(nom.plot$Mod,cut.time = c(365,365*3,365*5))

save.image(file='20220530_LUAD_APOBEC_final.RData')


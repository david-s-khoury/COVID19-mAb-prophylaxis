
# computing the geometric mean concentration ------------------------------

fun_geomean_conc <- function(days,concs,time_min,time_max){
  # inputs:
  # days      days concentration was measured
  # concs     measured concentration
  # time_min  lower time limit for geometric mean
  # time_max  upper time limit for geometric mean
  exp(pracma::trapz(c(time_min,days[days>time_min & days<time_max],time_max),
            c(approx(days,log(concs),xout=time_min)$y,log(concs[days>time_min & days<time_max]),approx(days,log(concs),xout=time_max)$y))/
        (time_max-time_min))
  }


# visualization of concentration and efficacy -----------------------------

# Plot for each study
conc.eff.plot <- function(trial_name,conc_log10=TRUE,same_scale=TRUE,
                          point_size=1,errorbar_w=5,errorbar_h=3){
  if(trial_name%in%c("Levin","Isa")){
    conc.tmp <- data_conc[data_conc$trial==trial_name,]
    if(conc_log10){
      conc.tmp$conc <- log10(conc.tmp$conc)
    }
    if(same_scale){
      conc.min.tmp <- ifelse(conc_log10,log10(min(data_conc$conc[data_conc$conc>0])),min(data_conc$conc[data_conc$conc>0]))
      conc.max.tmp <- ifelse(conc_log10,log10(max(data_conc$conc)),min(data_conc$conc))
    }else{
      conc.min.tmp <- min(conc.tmp$conc[is.finite(conc.tmp$conc)])
      conc.max.tmp <- max(conc.tmp$conc)
    }
    coeff <- 100/(conc.max.tmp-conc.min.tmp) # coefficient for second y-axis
    data.tmp <- data_all[data_all$trial==trial_name & !is.na(data_all$efficacy),]
    x.max <- max(data.tmp$day.max)
    # breaks and labels:
    if(same_scale){
        breaks.tmp <- log10(c(0.1,0.2,0.5,1,2,5,10,20,50,100,200))
        labels.tmp <- c(0.1,0.2,0.5,1,2,5,10,20,50,100,200)
    }else{
      if(trial_name=="Isa"){
        breaks.tmp <- log10(c(0.1,0.2,0.5,1,2,5,10,20,50,100,200))
        labels.tmp <- c(0.1,0.2,0.5,1,2,5,10,20,50,100,200)
      }else if(trial_name=="Levin"){
        breaks.tmp <- log10(seq(6,24,by=2))
        labels.tmp <- seq(6,24,by=2)
      }
    }
    # title:
    if(trial_name=="Isa"){
      title.tmp <- "Casirivimab + imdevimab (repeated administration every 4 weeks)"
    }else if(trial_name=="Levin"){
      title.tmp <- "Cilgavimab + tixagevimab"
    }
    
    ggplot(data.tmp,aes(x=day,y=efficacy/coeff+conc.min.tmp)) +
      # concentration data over time:
      geom_line(data=conc.tmp,inherit.aes=FALSE,aes(x=day,y=conc)) + 
      # efficacy estimates and 95% CIs:
      geom_errorbar(aes(ymin=efficacy.lower/coeff+conc.min.tmp, ymax=efficacy.upper/coeff+conc.min.tmp),width=errorbar_w,size=0.25,color="dodgerblue") +
      geom_errorbarh(aes(xmin=day.min, xmax=day.max),height=errorbar_h/coeff,size=0.25,color="dodgerblue") +
      geom_point(size=point_size,color="dodgerblue") +
      # axis, theme, etc.:
      scale_y_continuous(name = "Concentration [mg/L]", breaks = breaks.tmp, labels = labels.tmp,
                         sec.axis = sec_axis(~.*1, breaks=seq(0,100,by=25)/coeff+conc.min.tmp, labels=seq(0,100,by=25), name="Efficacy [%]")) +
      coord_cartesian(xlim = c(0,x.max), ylim=c(conc.min.tmp,conc.max.tmp)) +
      labs(x="Time [days]", title = title.tmp) +
      theme_bw() +
      theme(legend.key.size = unit(0.4, 'cm'),panel.grid = element_line(colour="gray95",size = 0.2),
            axis.title.y.right = element_text(color = "dodgerblue"),axis.text.y.right = element_text(color = "dodgerblue"))
    
  }else if(trial_name=="O'Brien & Herman"){ # plot the O'Brien & Herman data differently:
    conc.tmp <- data_conc[data_conc$trial=="O'Brien & Herman",]
    if(conc_log10){
      conc.tmp$conc <- log10(conc.tmp$conc)
    }
    conc.tmp <- dplyr::mutate(conc.tmp, obrien.data=ifelse(conc.tmp$day<=168,1,0)) # O'Brien concentration data to day XX, then Herman data
    if(same_scale){
      conc.min.tmp <- ifelse(conc_log10,log10(min(data_conc$conc[data_conc$conc>0])),min(data_conc$conc[data_conc$conc>0]))
      conc.max.tmp <- ifelse(conc_log10,log10(max(data_conc$conc)),min(data_conc$conc))
    }else{
      conc.min.tmp <- min(conc.tmp$conc[is.finite(conc.tmp$conc)])
      conc.max.tmp <- max(conc.tmp$conc)
    }
    coeff <- 100/(conc.max.tmp-conc.min.tmp) # coefficient for second y-axis
    data.tmp <- data_all[data_all$trial%in%c("O'Brien","Herman"),]
    x.max <- max(data.tmp$day.max)
    # breaks and labels:
    if(same_scale){
      breaks.tmp <- log10(c(0.1,0.2,0.5,1,2,5,10,20,50,100,200))
      labels.tmp <- c(0.1,0.2,0.5,1,2,5,10,20,50,100,200)
    }else{
      breaks.tmp <- log10(c(0.3,0.5,1,3,5,10,30,50,100))
      labels.tmp <- c(0.3,0.5,1,3,5,10,30,50,100)
    }
    # title:
    title.tmp <- "Casirivimab + imdevimab"
    
    my.shapes <- c(19,5)
    
    ggplot(data.tmp,aes(x=day,y=efficacy/coeff+conc.min.tmp,shape=trial)) +
      # concentration data over time:
      geom_line(data=conc.tmp[conc.tmp$obrien.data==1,],inherit.aes=FALSE,aes(x=day,y=conc)) + 
      geom_line(data=conc.tmp[conc.tmp$day>=168,],inherit.aes=FALSE,aes(x=day,y=conc),linetype=2) + 
      # efficacy estimates and 95% CIs:
      geom_errorbar(aes(ymin=efficacy.lower/coeff+conc.min.tmp, ymax=efficacy.upper/coeff+conc.min.tmp),width=errorbar_w,size=0.25,color="dodgerblue") +
      geom_errorbarh(aes(xmin=day.min, xmax=day.max),height=errorbar_h/coeff,size=0.25,color="dodgerblue") +
      geom_point(size=point_size,color="dodgerblue") +
      # axis, theme, etc.:
      scale_shape_manual(values=my.shapes,name="Study:") + 
      scale_linetype_manual(values=c(1,2),name="Study:") +
      scale_y_continuous(name = "Concentration [mg/L]", breaks = breaks.tmp, labels = labels.tmp,
                         sec.axis = sec_axis(~.*1, breaks=seq(0,100,by=25)/coeff+conc.min.tmp, labels=seq(0,100,by=25), name="Efficacy [%]")) +
      coord_cartesian(xlim = c(0,x.max), ylim=c(conc.min.tmp,conc.max.tmp)) +
      labs(x="Time [days]", title = title.tmp) +
      theme_bw() +
      theme(legend.key.size = unit(0.4, 'cm'),panel.grid = element_line(colour="gray95",size = 0.2),legend.position = "bottom",
            axis.title.y.right = element_text(color = "dodgerblue"),axis.text.y.right = element_text(color = "dodgerblue"))
    
  }else{# Schmidt trial
    conc.tmp <- data_conc[data_conc$trial=="Schmidt",]
    if(conc_log10){
      conc.tmp$conc <- log10(conc.tmp$conc)
    }
    if(same_scale){
      conc.min.tmp <- ifelse(conc_log10,log10(min(data_conc$conc[data_conc$conc>0])),min(data_conc$conc[data_conc$conc>0]))
      conc.max.tmp <- ifelse(conc_log10,log10(max(data_conc$conc)),min(data_conc$conc))
    }else{
      conc.min.tmp <- min(conc.tmp$conc[is.finite(conc.tmp$conc)])
      conc.max.tmp <- max(conc.tmp$conc)
    }
    coeff <- 100/(conc.max.tmp-conc.min.tmp) # coefficient for second y-axis
    data.tmp <- data_all[data_all$trial%in%c("Schmidt (delta)","Schmidt (omicron)"),]
    data.tmp$trial <- ifelse(data.tmp$trial=="Schmidt (delta)","delta","omicron")
    x.max <- max(data.tmp$day.max)
    # breaks and labels:
    if(same_scale){
      breaks.tmp <- log10(c(0.1,0.2,0.5,1,2,5,10,20,50,100,200))
      labels.tmp <- c(0.1,0.2,0.5,1,2,5,10,20,50,100,200)
    }else{
      breaks.tmp <- log10(c(6,10,15,20,25,30,35))
      labels.tmp <- c(6,10,15,20,25,30,35)
    }
    # title:
    title.tmp <- "Adintrevimab"
    
    my.shapes <- c(19,2)
    
    ggplot(data.tmp,aes(x=day,y=efficacy/coeff+conc.min.tmp,shape=trial)) +
      # concentration data over time:
      geom_line(data=conc.tmp,inherit.aes=FALSE,aes(x=day,y=conc),linetype=2) + 
      # efficacy estimates and 95% CIs:
      geom_errorbar(aes(ymin=efficacy.lower/coeff+conc.min.tmp, ymax=efficacy.upper/coeff+conc.min.tmp),width=errorbar_w,size=0.25,color="dodgerblue") +
      geom_errorbarh(aes(xmin=day.min, xmax=day.max),height=errorbar_h/coeff,size=0.25,color="dodgerblue") +
      geom_point(size=point_size,color="dodgerblue") +
      # axis, theme, etc.:
      scale_shape_manual(values=my.shapes,name="Variant:") + 
      scale_y_continuous(name = "Concentration [mg/L]", breaks = breaks.tmp, labels = labels.tmp,
                         sec.axis = sec_axis(~.*1, breaks=seq(0,100,by=25)/coeff+conc.min.tmp, labels=seq(0,100,by=25), name="Efficacy [%]")) +
      coord_cartesian(xlim = c(0,x.max), ylim=c(conc.min.tmp,conc.max.tmp)) +
      labs(x="Time [days]", title = title.tmp) +
      theme_bw() +
      theme(legend.key.size = unit(0.4, 'cm'),panel.grid = element_line(colour="gray95",size = 0.2),legend.position = "bottom",
            axis.title.y.right = element_text(color = "dodgerblue"),axis.text.y.right = element_text(color = "dodgerblue"))
    
  }
}


# efficacy functions ------------------------------------------------------

# other efficacy functions:
efficacy <- function(d,par,model=use_efficacy_model){
  if(model=="logistic"){
    # 3 parameters: maximum (transform specified in max_transform), slope (log-transformed), IC50 (log-transformed)
    m <- ifelse(max_transform=="exp",exp(par[1]),ifelse(max_transform=="1-exp",1-exp(par[1]),ifelse(max_transform=="1-abs",1-abs(par[1]),par[1])))
    m/(1+(2*m-1)*exp(-exp(par[2])*(log10(d)-log10(exp(par[3])))))
  }else if(model=="single hit"){
    # 1 parameter for slope
    1-exp(-par*log10(d))
  }else if(model=="powerlaw"){
    # parameters: slope (as for single hit model), exponent
    1-exp(-par[1]*log10(d)^par[2])
  }else if(model=="threshold"){
    # parameters: threshold for jump, efficacy before jump, efficacy after jump
    ifelse(d<par[1],par[2],par[3])
  }else if(model=="slope threshold"){
    # parameters: threshold for jump, efficacy before jump, efficacy after jump, slope before jump
    ifelse(d<par[1],par[2]+par[4]*log10(d),par[3])
  }else if(model=="double logistic"){
    # parameters: threshold, maximum efficacy, slope before threshold, slope after threshold
    # same transformations as for the logistic model, i.e. log-transform on slopes and specified transform on maximum
    m <- ifelse(max_transform=="exp",exp(par[2]),ifelse(max_transform=="1-exp",1-exp(par[2]),ifelse(max_transform=="1-abs",1-abs(par[2]),par[2])))
    ifelse(d<par[1],m/(m+exp(-exp(par[3])*log10(d))),m/(m+exp(-exp(par[3])*log10(par[1])-exp(par[4])*(log10(d)-log10(par[1])))))
  }else if(model=="logistic with slope 1"){
    # parameters: maximum efficacy, IC50
    # same transformations as for the logistic model, i.e. log-transformed IC50 and specified transform on maximum
    m <- ifelse(max_transform=="exp",exp(par[1]),ifelse(max_transform=="1-exp",1-exp(par[1]),ifelse(max_transform=="1-abs",1-abs(par[1]),par[1])))
    m/(1+(2*m-1)*exp(-(log10(d)-log10(exp(par[2])))))
  }else if(model=="logistic with max 1"){
    # parameters: slope, IC50
    # same transformations as for the logistic model, i.e. log-transformed slope and IC50
    1/(1+exp(-exp(par[1])*(log10(d)-log10(exp(par[2])))))
  }
}
efficacy <- Vectorize(efficacy, vectorize.args = "d")

# inverse logistic efficacy function: computes concentration for a given efficacy
inv_eff_logist <- function(eff,par){
  # input parameters:
  # eff   efficacy (between 0 and 1, not in percent)
  # par   parameters of the efficacy function (same as for the logistic efficacy function)
  m <- ifelse(max_transform=="exp",exp(par[1]),ifelse(max_transform=="1-exp",1-exp(par[1]),ifelse(max_transform=="1-abs",1-abs(par[1]),par[1])))
  ifelse(eff>=m | eff<=0,NA,exp(par[3])*10^(-(log((m-eff)/(eff*(2*m-1))))/(exp(par[2]))))
}

# inverse logistic efficacy function with maximum 1: computes concentration for a given efficacy
inv_eff_logist_max_1 <- function(eff,par){
  # input parameters:
  # eff   efficacy (between 0 and 1, not in percent)
  # par   parameters of the efficacy function (same as for the logistic efficacy function with maximum 1)
  ifelse(eff>=1 | eff<=0,NA,exp(par[2])*10^(-(log((1-eff)/eff))/(exp(par[1]))))
}

# negative log-likelihoods function for various efficacy functions:
# parameters: model dependent
nllh <- function(par,data_fit,model){-sum(log(dbinom(data_fit$e.control,data_fit$n.control,par[c(1:nrow(data_fit))])))-
    sum(log(dbinom(data_fit$e.treatment,data_fit$n.treatment,par[c(1:nrow(data_fit))]*(1-efficacy(data_fit$geomean.conc.ic50,par[(nrow(data_fit)+1):length(par)],model)))))}


# for multiple imputation -------------------------------------------------

# sample a day from an interval:
sample.int.vec <- function(day_min,day_max){sapply(c(1:length(day_min)),
                                                   function(x){day_min[x]-1+sample.int(day_max[x]-day_min[x]+1,size=1)})}

# find the Ab concentration on that day:
fun.conc <- function(data_imputed,data_conc_daily){sapply(c(1:nrow(data_imputed)),
                                                          function(x){data_conc_daily$conc[data_conc_daily$day==data_imputed$day[x] & 
                                                                                             data_conc_daily$trial==ifelse(data_imputed$trial[x]%in%c("O'Brien","Herman"),"O'Brien & Herman",
                                                                                                                           ifelse(data_imputed$trial[x]%in%c("Schmidt (delta)","Schmidt (omicron)"),"Schmidt",
                                                                                                                                  data_imputed$trial[x]))]})}

# sample IC50s:
sample.ic50s <- function(Means_table,data_imputed){sapply(c(1:nrow(data_imputed)),
                                                          function(x){exp(rnorm(1,mean=Means_table$log_IC50[Means_table$Antibodies==data_imputed$drug[x] & Means_table$VariantCategory==ifelse(data_imputed$trial[x]=="Schmidt (delta)","Delta",ifelse(data_imputed$trial[x]=="Schmidt (omicron)","Omicron/BA.1","Wild Type"))],
                                                                                sd=Means_table$log_IC50_SE[Means_table$Antibodies==data_imputed$drug[x] & Means_table$VariantCategory==ifelse(data_imputed$trial[x]=="Schmidt (delta)","Delta",ifelse(data_imputed$trial[x]=="Schmidt (omicron)","Omicron/BA.1","Wild Type"))]))})}

# variance within imputations:
fun.var.within <- function(hessians){Reduce('+',lapply(hessians,function(x){solve(x)[ncol(x)-c((length(LowerB_all_models[[which(names(LowerB_all_models)==use_efficacy_model)]])-1):0),ncol(x)-c((length(LowerB_all_models[[which(names(LowerB_all_models)==use_efficacy_model)]])-1):0)]}))/(length(hessians))}

# variance between imputations:
fun.var.between <- function(param_best){Reduce('+',lapply(as.list(c(1:nrow(param_best))),function(x){as.matrix(param_best[x,]-colMeans(param_best),ncol=1)%*%t(as.matrix(param_best[x,]-colMeans(param_best)))}))/(nrow(param_best)-1)}


# for efficacy prediction -------------------------------------------------

# make a figure of the Ab concentration data, fitted half-life and half-life from the literature:
fun.ab.conc.half.life <- function(ab_name,lin_mod,peak_conc,peak_time,half_life,data_conc){
  # input parameters:
  # ab_name     name of the mAb to plot
  # lin_mod     linear model fit to log-transformed concentration data after peak
  # peak_conc   concentration at peak
  # peak_time   time of peak concentration
  # half_life   half-life of mAb reported in the literature
  # data_conc   concentration data
  lines.plot <- data.frame(group=c("fitted to data","from literature"),
                           x=c(0,0),
                           xend=c(380,380),
                           y=c(exp(lin_mod$coefficients[1]),exp(log(peak_conc)+peak_time*log(2)/half_life)),
                           yend=c(exp(lin_mod$coefficients[1]+380*lin_mod$coefficients[2]),
                                  exp(-380*log(2)/half_life+log(peak_conc)+log(2)/half_life*peak_time)))
  
  my_xlim <- c(1,ifelse(ab_name=="casirivimab + imdevimab",250,
                        ifelse(ab_name=="cilgavimab + tixagevimab",200,
                               ifelse(ab_name=="adintrevimab",370,NA))))
  my_ylim <- c(ifelse(ab_name=="casirivimab + imdevimab",0.2,ifelse(ab_name=="cilgavimab + tixagevimab",5,ifelse(ab_name=="adintrevimab",5,NA))),
               ifelse(ab_name=="casirivimab + imdevimab",110,ifelse(ab_name=="cilgavimab + tixagevimab",26,ifelse(ab_name=="adintrevimab",35,NA))))
  trial <- ifelse(ab_name=="casirivimab + imdevimab","O'Brien & Herman",
                  ifelse(ab_name=="cilgavimab + tixagevimab","Levin",
                         ifelse(ab_name=="adintrevimab","Schmidt",NA)))
  
  ggplot(data=data_conc[data_conc$trial==trial & data_conc$conc>0,],aes(x=day,y=conc)) +
    geom_point() +
    # add model fit and concentration data using the half-lives from the literature:
    geom_segment(data=lines.plot,inherit.aes = FALSE, aes(x=x,xend=xend,y=y,yend=yend,colour=group)) + 
    scale_colour_manual(values = c("red","blue"), name="Fit with half-life") + 
    # axis, theme, etc.:
    scale_y_log10(name = "Concentration [mg/L]", expand = c(2e-2,2e-2)) +
    scale_x_continuous(expand = c(2e-2,2e-3)) +
    coord_cartesian(xlim = my_xlim,ylim=my_ylim) +
    labs(x="Days after treatment", title = glue::glue("Concentration data fit for ",ab_name)) +
    theme_bw() +
    theme(legend.key.size = unit(0.4, 'cm'),panel.grid = element_line(colour="gray95",size = 0.2),legend.position = 'bottom')
}

# time (since treatment) above a efficacy threshold:
time.above.limit <- function(ab_name,variant,eff_thresh,peak_conc,lin_mod,par_eff,Means_table,ic50.fold.change=1,CI=FALSE,eff_IC50=exp(par_eff[3])){
  # input parameters:
  # ab_name           mAb to consider
  # variant           name of the variant category
  # peak_conc         peak Ab concentration
  # lin_mod           linear model fit to log-transformed concentration data after peak
  # eff_thresh        efficacy threshold
  # par_eff           parameters for efficacy function
  # Means_table       table with IC50 of different mAbs against different variants
  # ic50.fold.change  IC50-fold change against the variant specified by "variant"
  # CI                compute a CI bound (TRUE) or the central estimate (FALSE)
  # eff_IC50          IC50 of the efficacy function, if the CI is to be computed, replace with the CI of the efficacy IC50
  if(any(Means_table$Antibodies==ab_name & Means_table$VariantCategory==variant)){
    meta.ic50 <- ic50.fold.change*Means_table$IC50[Means_table$Antibodies==ab_name & Means_table$VariantCategory==variant]
    eff_par <- par_eff
    if(CI){ # if the CI is to be computed, use eff_IC50 as the IC50 for the efficacy function
      if(use_efficacy_model=="logistic"){
        eff_par[3] <- log(eff_IC50)
      }else if(use_efficacy_model=="logistic with max 1"){
        eff_par[2] <- log(eff_IC50)
      }
    }
    peak.efficacy <- efficacy(peak_conc/(meta.ic50*1e-3),eff_par,model = use_efficacy_model)
    # concentration at which the eff_thresh is reached in mg/L:
    thresh.conc <- ifelse(use_efficacy_model=="logistic",inv_eff_logist(eff_thresh,eff_par)*meta.ic50*1e-3,
                          ifelse(use_efficacy_model=="logistic with max 1",inv_eff_logist_max_1(eff_thresh,eff_par)*meta.ic50*1e-3,NA))
    as.numeric(ifelse(peak.efficacy<eff_thresh,0,(log(thresh.conc)-lin_mod$coefficients[1])/lin_mod$coefficients[2]))
  }else{
    NA
  }
}
time.above.limit <- Vectorize(time.above.limit,vectorize.args = "ic50.fold.change")
  

### Chisquared goodness of fit test statistic
chisquared_gof_statistic <- function(par,data_fit,model=use_efficacy_model){sum((((data_fit$n.control*(par[c(1:nrow(data_fit))])-data_fit$e.control)^2)/(data_fit$n.control*(par[c(1:nrow(data_fit))]))+
                                                    ((data_fit$n.treatment*(par[c(1:nrow(data_fit))]*(1-efficacy(data_fit$geomean.conc.ic50,par[(nrow(data_fit)+1):length(par)],model)))-data_fit$e.treatment)^2)/(data_fit$n.treatment*(par[c(1:nrow(data_fit))]*(1-efficacy(data_fit$geomean.conc.ic50,par[(nrow(data_fit)+1):length(par)],model))))))}


# efficacy by time using the fit to the concentration data:
efficacy.by.time <- function(time,ab_name,variant,lin_mod,ic50.fold.change=1,meta_analysis=Means_table,par_eff=param_best){
  # input parameters:
  # time              time for which to compute the efficacy
  # ab_name           name of the mAb
  # variant           SARS-CoV2 variant category
  # lin_mod           linear model fit to the log-transformed concentration data
  # ic50.fold.change  IC50-fold change against the variant specified by "variant"
  # meta_analysis     result of the IC50 meta-analysis
  # par_eff           parameters or the efficacy function
  if(any(meta_analysis$Antibodies==ab_name & meta_analysis$VariantCategory==variant)){
    efficacy(exp(time*lin_mod$coefficients[2]+lin_mod$coefficients[1])/(ic50.fold.change*meta_analysis$IC50[meta_analysis$Antibodies==ab_name & meta_analysis$VariantCategory==variant]*1e-3),par_eff)
  }else{
    NA
  }
}

# visualization of the efficacy over time:
vis.eff.by.time <- function(day_max,ab_name,lin_mod,variants=c("Wild Type","Omicron/BA.1","Omicron/BA.2"),fold.changes=c(5,10,50)){
  # input parameters:
  # day_max       maximum time duration for visualization
  # ab_name           name of the mAb
  # lin_mod           linear model fit to the log-transformed concentration data
  # variants      variants to visualize
  # fold.changes  fold changes of IC50 to visualize (compared to wild type)
  days <- seq(0,day_max,length.out=1e3)
  eff.variants <- sapply(c(1:length(variants)),function(x){efficacy.by.time(time=days,ab_name,variant=variants[x],lin_mod)})
  eff.variants <- as.vector(eff.variants)
  eff.fold <- sapply(c(1:length(fold.changes)),function(x){efficacy.by.time(time=days,ab_name,variant="Wild Type",ic50.fold.change=fold.changes[x],lin_mod)})
  eff.fold <- as.vector(eff.fold)
  data.vis <- data.frame(day=days,
                         efficacy=c(eff.variants,eff.fold),
                         escape=factor(c(rep(variants,rep(length(days),length(variants))),rep(paste("fold.",fold.changes,sep=""),rep(length(days),length(fold.changes))))),levels=c(variants,paste("fold.",fold.changes,sep="")))
  
  # colors:
  colours_red <- colorRampPalette(c("red4", "red"))
  colours_green <- colorRampPalette(c("olivedrab4", "olivedrab1"))
  my.colours.fold.escape <- c("black",colours_red(length(variants)-1),colours_green(length(fold.changes)))[order(unique(data.vis$escape))]
  
  ggplot(data.vis,aes(x=day,y=100*efficacy,color=escape))+
    geom_line() +
    # axis, theme, etc.:
    scale_color_manual(values=my.colours.fold.escape, labels=sort(unique(data.vis$escape)),
                       name="Level of escape [fold IC50 change]:") +
    scale_y_continuous(name = "Efficacy [%]", expand = c(2e-2,2e-2)) +
    scale_x_continuous(expand = c(2e-2,2e-3)) +
    coord_cartesian(ylim = c(1,100)) +
    labs(x="Days after treatment", title = glue::glue(ab_name)) +
    theme_bw() +
    theme(legend.key.size = unit(0.4, 'cm'),panel.grid = element_line(colour="gray95",size = 0.2),legend.position = 'bottom')
}

# fold-in vitro IC50 change to give 30 days of at least 50% protection (compared to Wild Type):
fun.fold.change <- function(ab_name,lin_mod,meta_analysis=Means_table,table_eff_par=table_par_best,CI=FALSE,eff_IC50=table_eff_par$estimate[table_eff_par$parameter=="IC50"]){
  # input parameters:
  # ab_name           name of the mAb
  # lin_mod           linear model fit to the log-transformed concentration data
  # meta_analysis     result of the IC50 meta-analysis
  # table_eff_par     table with efficacy function parameters
  # CI                compute a CI bound (TRUE) or the central estimate (FALSE)
  # eff_IC50          IC50 of the efficacy function, if the CI is to be computed, replace with the CI of the efficacy IC50
  if(CI){
    as.numeric(exp(30*lin_mod$coefficients[2]+lin_mod$coefficients[1])/
                 (eff_IC50*1e-3*meta_analysis$IC50[meta_analysis$Antibodies==ab_name & meta_analysis$VariantCategory=="Wild Type"]))
  }else{
    as.numeric(exp(30*lin_mod$coefficients[2]+lin_mod$coefficients[1])/
                 (table_eff_par$estimate[table_eff_par$parameter=="IC50"]*1e-3*
                    meta_analysis$IC50[meta_analysis$Antibodies==ab_name & meta_analysis$VariantCategory=="Wild Type"]))
  }
}


# for the comparison of vaccine and mAb data ------------------------------

# get the index of parameters that are the same in the vaccine and mAb data:
get.same.par.index <- function(same.par){
  if(use_efficacy_model=="logistic"){
    # parameters of the logistic model: maximum (m), slope (k), IC50 (ic50)
    tmp <- unlist(ifelse(grepl("m",same.par),ifelse(grepl("k",same.par),ifelse(grepl("ic50",same.par),list(c(1,2,3)),list(c(1,2))),
                                             ifelse(grepl("ic50",same.par),list(c(1,3)),list(c(1)))),
                  ifelse(grepl("k",same.par),ifelse(grepl("ic50",same.par),list(c(2,3)),list(c(2))),
                         ifelse(grepl("ic50",same.par),list(c(3)),NA))))
  }else if(use_efficacy_model=="logistic with max 1"){
    # parameters of the logistic model with max 1: slope (k), IC50 (ic50)
    tmp <- unlist(ifelse(grepl("k",same.par),ifelse(grepl("ic50",same.par),list(c(1,2)),1),ifelse(grepl("ic50",same.par),2,NA)))
  }else{
    tmp <- NA
  }
  if(any(is.na(tmp))){
    integer()
  }else{
    tmp
  }
}

# get the index of parameters that are different in the vaccine and the mAb data:
get.different.par.index <- function(same.par){
  if(length(get.same.par.index(same.par))==0){
    c(1:length(LowerB_all_models[[which(names(LowerB_all_models)==use_efficacy_model)]]))
  }else{
    which(!c(1:length(LowerB_all_models[[which(names(LowerB_all_models)==use_efficacy_model)]]))%in%get.same.par.index(same.par))
  }
}

# negative log-likelihood function for fitting the vaccine and mAb data simultaneously:
nllh.vac.mAb <- function(par,same.par,data_fit=combined_data_vac_mAbs_fit,model=use_efficacy_model){
  # parameters for the efficacy of vaccines and mAbs:
  par.vac <- par[nrow(data_fit)+c(1:length(LowerB_all_models[[which(names(LowerB_all_models)==model)]]))]
  par.diff.mAb <- rep(0,length(LowerB_all_models[[which(names(LowerB_all_models)==model)]]))
  if(length(get.different.par.index(same.par))>0){
    par.diff.mAb[get.different.par.index(same.par)] <- tail(par,length(get.different.par.index(same.par[i])))
  }
  par.mAb <- par.vac+par.diff.mAb
  
  # negative log-likelihood:
  -sum(log(dbinom(data_fit$e.control,data_fit$n.control,par[c(1:nrow(data_fit))])))-
    sum(log(dbinom(data_fit$e.treatment[data_fit$type=="vaccine"],data_fit$n.treatment[data_fit$type=="vaccine"],par[which(data_fit$type=="vaccine")]*(1-efficacy(data_fit$fold.conv[data_fit$type=="vaccine"],par.vac)))))-
    sum(log(dbinom(data_fit$e.treatment[data_fit$type=="mAb"],data_fit$n.treatment[data_fit$type=="mAb"],par[which(data_fit$type=="mAb")]*(1-efficacy(data_fit$fold.conv[data_fit$type=="mAb"],par.mAb)))))
}

# model comparison with the likelihood ratio test:
lr.pval <- function(same.par1,same.par2,data_fit=combined_data_vac_mAbs_fit){ # models to be compared need to be nested!
  ind1 <- which(gsub("same.","",names(RandomInitialEstimate_mAb))==same.par1)
  ind2 <- which(gsub("same.","",names(RandomInitialEstimate_mAb))==same.par2)
  
  nllh1 <- min(RandomInitialEstimate_mAb[[ind1]]$nllh)
  nllh2 <- min(RandomInitialEstimate_mAb[[ind2]]$nllh)
  
  # number of parameters:
  # baseline risks plus number of fitted parameters (parameters for vaccines & mAbs minus parameters that are the same)
  n.par1 <- nrow(data_fit)+(2*(ncol(RandomInitialEstimate_mAb[[ind1]])-3)-length(get.same.par.index(same.par1))) 
  n.par2 <- nrow(data_fit)+(2*(ncol(RandomInitialEstimate_mAb[[ind2]])-3)-length(get.same.par.index(same.par2))) 
  
  if(length(grep(paste(get.same.par.index(same.par1),collapse = ""),paste(get.same.par.index(same.par2),collapse = "")))>0){
    lr <- 2*(nllh2-nllh1)
    p.val <- pchisq(lr, df = n.par1-n.par2, lower.tail = FALSE)
  }else if(length(grep(paste(get.same.par.index(same.par2),collapse = ""),paste(get.same.par.index(same.par1),collapse = "")))>0){
    lr <- 2*(nllh1-nllh2)
    p.val <- pchisq(lr, df = n.par2-n.par1, lower.tail = FALSE)
  }else{
    p.val <- NA
  }
  p.val
}

# visualize fit to the vaccine and mAb data:
vis.vac.mAb.fit <- function(same.par,data_vis=combined_data_vac_mAbs){
  title.tmp <- ifelse(same.par=="k.ic50","Same k and c50",ifelse(same.par=="k","Same k, different c50",ifelse(same.par=="ic50","Same c50, different k",ifelse(same.par=="","Different k and c50","Unknown"))))
  
  # visualization:
  my.colours <- c("royalblue3","orangered") #c("goldenrod1","darkorange","red","red3","lightpink","hotpink","mediumpurple","dodgerblue","turquoise1","palegreen2","springgreen4","darkviolet","purple4","black")
  my.shapes <- c(16,17)
  dose.range <- c(min(data_vis$fold.conv.lower),max(data_vis$fold.conv.upper))
  
  ggplot(data_vis,aes(x=fold.conv,y=efficacy,colour=type,shape=type,alpha=as.factor(earliest.time))) +
    # show fit and 95% CI for vaccine and mAb data:
    geom_ribbon(data=models_vac_mAb[[which(gsub("same.","",names(models_vac_mAb))==same.par)]],inherit.aes=FALSE,
                aes(x=dose,ymin=100*fit_vac_lower,ymax=100*fit_vac_upper),fill="orangered", alpha = 0.1) +
    geom_ribbon(data=models_vac_mAb[[which(gsub("same.","",names(models_vac_mAb))==same.par)]],inherit.aes=FALSE,
                aes(x=dose,ymin=100*fit_mAb_lower,ymax=100*fit_mAb_upper),fill="royalblue3", alpha = 0.1) +
    geom_line(data=models_vac_mAb[[which(gsub("same.","",names(models_vac_mAb))==same.par)]],inherit.aes=FALSE,aes(x=dose,y=100*fit_vac),color="orangered") +
    geom_line(data=models_vac_mAb[[which(gsub("same.","",names(models_vac_mAb))==same.par)]],inherit.aes=FALSE,aes(x=dose,y=100*fit_mAb),color="royalblue3") +
    # efficacy estimates and 95% CIs:
    geom_errorbar(aes(ymin=efficacy.lower, ymax=efficacy.upper),width=0.05,size=0.25) +
    geom_errorbarh(aes(xmin=fold.conv.lower, xmax=fold.conv.upper),height=2,size=0.25) +
    geom_point(size=2) +
    scale_shape_manual(values=my.shapes, name="Treatment:") +
    scale_alpha_manual(values=c(1,0.35), labels=c("no","yes"), name="Earliest time interval:") +
    scale_color_manual(values=my.colours, name="Treatment:") +
    # axis, theme, etc.:
    scale_x_log10() +
    scale_y_continuous(minor_breaks = seq(0, 100, 20),breaks=seq(0, 100, 20),labels=seq(0, 100, 20),expand = c(2e-2,1e-1)) +
    coord_cartesian(ylim = c(0,100), xlim = dose.range) + # c(1e-2,100)
    labs(x="Neutralization [fold convalescent]", y="Efficacy [%]", title = title.tmp) +
    theme_bw() +
    theme(legend.key.size = unit(0.4, 'cm'),panel.grid = element_line(colour="gray95",size = 0.2))
}


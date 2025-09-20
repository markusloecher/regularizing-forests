library(ggplot2)

#I give up on pickle, what a massive waste of time, what an inferior package!!
#Instead I will read from directories created by my own custom python `pdDict_tocsv()` function:
  
 
#check to make sure there is no index column in the saved csv files!
list_from_pdDict_tocsv = function(path="sim200_dof_rho035", 
                                  row.names = NULL, 
                                  header=FALSE,
                                  verbose=0){
  cwd= getwd()
  if (verbose >1){
    print(list.files("./"))
    print(cwd)
  }
  #stopifnot(path %in% list.files("./"))
  setwd(path)
  data = list()
  for (f in list.files("./", pattern="*.csv")) {
    ff =gsub(".csv", "", f, fixed=TRUE)
    data[[ff]] = read.csv(f,row.names = row.names, header = header)
    colnames(data[[ff]]) = gsub("X", "", colnames(data[[ff]]))
    data[[ff]]$method = ff
  }
  if (verbose >2) browser()
  setwd(cwd)
  return(data)
}

#copied from https://github.com/syzhou5/randomness-as-regularization/blob/master/DoF/DoF_example.R
plt_figure3_orig = function(){
  #lgd.names <- c("Random Forest: mtry=1","Random Forest: mtry=0.67","Random Forest: mtry=0.33","Random Forest: mtry=0.1")
  set.seed(1)
  n       <- 1000   # training size
  p       <- 10     # feature dimension
  nval    <- 1000   # 
  
  snr     <- 3.52
  maxnd   <- ceiling(seq(2,n/5,length.out = 9))
  # Draw the plot of DoF v.s. maxnodes
  dat <- data.frame(x=rep(maxnd,N),y=matrix(df.all,ncol=1),method=factor(rep(lgd.names,each=lmax)))
  gp<- ggplot(dat,aes(x=x,y=y,color=method))+
    ylab("Degree of Freedom")+
    xlab("maxnodes")+
    ggtitle("MARSadd Model")+
    geom_line(lwd=1) + geom_point(aes(shape = method), size =3 ) +
    geom_line(aes(y=x),linetype="dotted")+
    theme_bw() +
    theme(legend.justification = c(1,0),legend.position = c(0.99,0.01))+
    theme(plot.title = element_text(hjust = 0.5,size=rel(2.2)), 
          legend.key.size = unit(0.06, "npc"), 
          legend.spacing.y = unit(.02, "npc"),
          legend.text = element_text(size=rel(1.5)),
          legend.title = element_blank(),
          axis.title = element_text(hjust=0.5,size=rel(2)),
          axis.text = element_text(size=rel(1.7)))
}

#modified version to fit our data structure
plt_figure3_from_python = function(
    data_list, ##<< list from a call to list_from_pdDict_tocsv()
    scales = "free", #determines the y scales
    verbose = 1
){
  
  data_list = do.call("rbind.data.frame", data_list)
  dat = pivot_longer(data_list, 1:5, names_to = "lambda", values_to = "dof")
  dat$method = gsub("lin_func", "Low", dat$method)
  
  if (verbose>1) return(dat) #browser()
  gp = ggplot(dat,aes(x=maxnd,y=dof,color=lambda))+geom_line(lwd=1)+ 
    geom_point(aes(shape = lambda), size =1 )+ 
    facet_wrap(~ method, scales = scales) + 
    ylab("Degree of Freedom")+ xlab("maxnodes") + 
    geom_abline(intercept = 0, slope = 1, linewidth = 0.5,linetype="dotted")+
    theme_bw()
    # theme(legend.justification = c(1,0),legend.position = c(0.99,0.01))+
    # theme(plot.title = element_text(hjust = 0.5,size=rel(2.2)), 
    #       legend.key.size = unit(0.06, "npc"), 
    #       legend.spacing.y = unit(.02, "npc"),
    #       legend.text = element_text(size=rel(1.5)),
    #       legend.title = element_blank(),
    #       axis.title = element_text(hjust=0.5,size=rel(2)),
    #       axis.text = element_text(size=rel(1.7)))
  return(gp)
}
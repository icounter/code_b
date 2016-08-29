##############################################################compare mean-variance wrt utility
####################first scenario no transcation cost no bayesian net######################################################
############################################################################################################################
mean_variance<-function(w,mu=mu1,A=A1,cov=asset_cov1){
  rett<-w%*%mu-1/2*A*w%*%cov%*%w
  return(-rett)
}
mean_variance_grad<-function(w,mu=mu1,A=A1,cov=asset_cov1) nl.grad(w,mean_variance,heps=1e-10)
mean_variance2<-function(w,mu,beta1,cov,w_now,trans_cost,finan_cost,haircut,real_finance_weight,principal1){
  print("mean_variance initialize")
  tcost<-cost(w_now,w,trans_cost,finan_cost,haircut,real_finance_weight,principal1)
  rett<-w%*%mu-1/2*beta1*w%*%cov%*%w+tcost
  return(-rett)
}
####using this to generate N random variables
generate_N_rand<-function(N_p,corr,margins_r,paramMargins_r,N){
  #set.seed(100)
  ####only generate half of the total required paths 
  N<-as.integer(N/2)
  ret<-c()
  sd1<-c()
  for(i in 1:length(paramMargins_r)){
    ret<-c(ret,paramMargins_r[i][[1]]$mean)
    sd1<-c(sd1,paramMargins_r[i][[1]]$sd)
  }
  ####generate normal copula
  myCop <- normalCopula(param=corr, dim = N_p, dispstr = "un")
  ####function to generate random paths
  set.seed(100)
  gmvd2<-mvdc(myCop,margins_r,paramMargins_r)
  Z2 <- rMvdc(N,gmvd2)
  ######################################################variance reduction technique make mean and variance consistent 
  Z3<-sweep(-sweep(Z2,2,ret),2,-ret)
  Z2<-rbind(Z2,Z3)
  a<-sweep(Z2,2,ret)
  b<-apply(a,2,function(x) sd(x))
  for(i in 1:nrow(a)){
    a[i,]=a[i,]/b*sd1
  }
  a[,ncol(a)]<-rep(0,nrow(a))
  Z2<-sweep(a,2,-ret)
  #colnames(Z2) <-parameter_name
  return(Z2)
}
expo_u<-function(x,lambda){
  u<- (1-exp(-lambda*x))/lambda
  return(u)
}
get_expo_u<-function(ret,lambda){
  return(expo_u(ret,lambda))
}
expo_find<-function(w1,randu=Z1,lambda=A1,sample_number=N_sample){
  # print(w1)
  exp_expo_ut<-mean(get_expo_u(randu %*% w1,lambda))
  return(-(exp_expo_ut))
}
get_rounds<-function(strings){
  return(as.double(gsub("iteration: ","",toString(strings[length(strings)]))))
}
#############################################################################################################
####this func is used to find suboptimal solutions when global search fails 
###########################################################################################################
find_optimal<-function(buffer,signal){
  solutionrounds<-get_rounds(buffer[grepl("iteration:", buffer)])
  solutionPathParamRaw = buffer[grepl("^\tx", buffer)]
  solutionPathParam<-c()
  for(i in 1:solutionrounds){
    temp<-eval(parse(text=paste0("c",gsub("\tx \\= ","",toString(solutionPathParamRaw[i])))))
    if(i==1)  solutionPathParam<-temp
    else solutionPathParam<-rbind(solutionPathParam,temp)
  }
  solutionPathParam<-as.matrix(solutionPathParam)
  solutionPathParamRaw_f_v = buffer[grepl("^\tf", buffer)]
  solutionPathParam_f_v<-c()
  for(i in 1:solutionrounds){
    temp<-eval(parse(text=gsub("\tf\\(x\\) \\= ","",toString(solutionPathParamRaw_f_v[i]))))
    if(i==1)  solutionPathParam_f_v<-temp
    else solutionPathParam_f_v<-c(solutionPathParam_f_v,temp)
  }
  if(signal[1]==1){
    solutionPathParamRaw_g_v = buffer[grepl("^\tg", buffer)]
    solutionPathParam_g_v<-c()
    for(i in 1:solutionrounds){
      if(substr(toString(solutionPathParamRaw_g_v[1]),nchar(toString(solutionPathParamRaw_g_v[1])),nchar(toString(solutionPathParamRaw_g_v[1])))==")"){
        temp<-eval(parse(text=paste0("c",gsub("\tg\\(x\\) \\= ","",toString(solutionPathParamRaw_g_v[i])))))
      }else{
        temp<-eval(parse(text=gsub("\tg\\(x\\) \\= ","",toString(solutionPathParamRaw_g_v[i]))))
      }
      if(i==1)  solutionPathParam_g_v<-temp
      else solutionPathParam_g_v<-rbind(solutionPathParam_g_v,temp)
    }
    solutionPathParam_g_v<-as.matrix(solutionPathParam_g_v)
  }else{
    solutionPathParam_g_v<-as.matrix((rep(0,solutionrounds))) 
  }
  if(signal[2]==1){
    solutionPathParamRaw_h_v = buffer[grepl("^\th", buffer)]
    solutionPathParam_h_v<-c()
    for(i in 1:solutionrounds){
      if(substr(toString(solutionPathParamRaw_h_v[1]),nchar(toString(solutionPathParamRaw_h_v[1])),nchar(toString(solutionPathParamRaw_h_v[1])))!=")"){
        temp<-eval(parse(text=gsub("\th\\(x\\) \\= ","",toString(solutionPathParamRaw_h_v[i]))))
      }else{
        temp<-eval(parse(text=paste0("c",gsub("\th\\(x\\) \\= ","",toString(solutionPathParamRaw_h_v[i])))))
      }
      if(i==1)  {
        solutionPathParam_h_v<-temp
      }else{
        solutionPathParam_h_v<-rbind(solutionPathParam_h_v,temp)
      }
    }
    solutionPathParam_h_v<-as.matrix(solutionPathParam_h_v)
  }else{
    solutionPathParam_h_v<-as.matrix(rep(0,solutionrounds))
  }
  list_h<-as.vector(which(apply(solutionPathParam_h_v,1,function(x) sum(ifelse(x<=tol_constraints_eq&&x>=-tol_constraints_eq,0,1)))==0))
  list_g<-as.vector(which(apply(solutionPathParam_g_v,1,function(x) sum(ifelse(x<=tol_constraints_ineq,0,1)))==0))
  if(length(intersect(list_h,list_g))==0){
    M<-which(solutionPathParam_f_v!=0)
    solutionPathParam_f_v<-solutionPathParam_f_v[M]
    solutionPathParam<-solutionPathParam[M,]
    w_temp<-solutionPathParam[which.min(solutionPathParam_f_v),]
    return(w_temp)
  }else{
    M<-intersect(list_h,list_g)
    M<-intersect(which(solutionPathParam_f_v!=0),M)
    if(length(M)==0){
      M<-which(solutionPathParam_f_v!=0)
      solutionPathParam_f_v<-solutionPathParam_f_v[M]
      solutionPathParam<-solutionPathParam[M,]
      w_temp<-solutionPathParam[which.min(solutionPathParam_f_v),]
    }else{
    solutionPathParam_f_v<-solutionPathParam_f_v[M]
    solutionPathParam<-solutionPathParam[M,]
    w_temp<-solutionPathParam[which.min(solutionPathParam_f_v),]
    }

    return(w_temp)
  }
}

expo_find_grad<-function(w1,randu=Z1,lambda=A1,sample_number=N_sample) nl.grad(w1,expo_find,heps=1e-10)
scenario_no_cost_no_BN<-function(assets_num,lambda,asset_ret,asset_var,asset_corr,sample_number,lower_bound,upper_bound){
  ####input variable:
  ###
  ###
  ##
  #########
  N<-assets_num
  A<-lambda
  mu<<-asset_ret
  asset_cor<-matrix(1,N,N)
  kk=1
  for(i in 1:(N-1)){
    for(j in (i+1):N){
       asset_cor[i,j]=asset_corr[kk]
       kk=kk+1
       asset_cor[j,i]=asset_cor[i,j]
    }
  }
  ###################################these global variables below are not redunant becuase we need them to calcualte the gradient and jacobian matrix(if needed) for functions
  asset_corr1<<-asset_corr
  asset_cor<<-asset_cor
  asset_cov<-r2cov(sd =sqrt(asset_var),R = asset_cor)
  asset_cov<<-asset_cov
  init_val<-rep(0,N)
  lower1<<-lower_bound
  upper1<<-upper_bound
  mu1<<-mu
  A1<<-A
  asset_cov1<<-asset_cov
  ###################################set optimize conditions use AUGLAG+SLSQP to do mean variance search(convex optimal local opimize is enough)
  local_opts <<- list( "algorithm" = "NLOPT_LD_SLSQP")
  opts <<- list( "algorithm" = "NLOPT_LD_AUGLAG",
                 "maxeval" = maxeval,
                 "local_opts" = local_opts,
                 "print_level" = print_level)
  res1<<-nloptr(x0=init_val,eval_f=mean_variance,eval_grad_f = mean_variance_grad,lb=lower1,ub=upper1,opts=opts,mu=mu,A=A,cov=asset_cov)
  #####define what's the copula and marginal distr
  margins_r<-rep("norm",N)
  paramMargins_r <- list()
  # Now the new experiments
  for(i in 1:N){
    paramMargins_r[[length(paramMargins_r)+1]] <- list(mean =asset_ret[i], sd =sqrt(asset_var[i]))
  }
  paramMargins_r<<-paramMargins_r
  margins_r<<-margins_r
  N_sample<<-sample_number
  Z<<-generate_N_rand(N,asset_corr,margins_r,paramMargins_r,N_sample)
  Z1<<-Z
  ###start expo optimize
  ###################################set optimize conditions first global search AUGLAG+MLSL 
  ################################and all the other conditions
  local_opts <<- list( "algorithm" = "NLOPT_GD_MLSL","xtol_abs"=rep(1e-20,NN))
  opts <<- list( "algorithm" = "NLOPT_LD_AUGLAG",
                 "xtol_abs" = xtol_abs,
                 "ftol_abs"=ftol_abs,
                 "maxeval" = maxeval,
                 "local_opts" = local_opts,
                 "print_level" = print_level)
  #######################buffer is a variable which stores the optimization information###################
  ############we need to use buffer to find the optimial global search's result point#########################
  # buffer<<-capture.output(res3<<-nloptr(x0=rep(0,N),eval_f=expo_find,eval_grad_f = expo_find_grad,lb=lower1,ub=upper1,opts=opts,randu=Z,lambda=A,sample_number=N_sample))
  opts$local_opts$algorithm="NLOPT_LD_SLSQP"
  res2<<-nloptr(x0=init_val,eval_f=expo_find,eval_grad_f = expo_find_grad,lb=lower1,ub=upper1,opts=opts,randu=Z,lambda=A,sample_number=N_sample)
  weights_mean=res1$solution
  weights_exp=res2$solution
  asset_name1=rep("asset_",N)
  for(i in 1:N){
    asset_name1[i]=paste0(asset_name1[i],i)
  }
  ##################calculate result of 2 different solutions(just to verify their utility values is very close)
  util_1<-expo_find(res1$solution,Z,A,sample_number)
  util_2<-expo_find(res2$solution,Z,A,sample_number)
  methods1=c(rep(paste0("mean_variance,utility=",sprintf("%.6f",-util_1)),N),rep(paste0("exponential utility,utility=",sprintf('%.6f',-util_2)),N))
  dimen<-2*length(weights_mean)
  weights1<<-data.frame(
    methods=methods1,
    asset_name=rep(asset_name1,2),
    weight=c(weights_mean,weights_exp))
  return(weights1)
}
#############################this is a function to preprocess some of the variables need to input into scenario_no_cost_no_BN function
#####################same idea for the other functions
call_scenario_no_cost_no_BN<-function(assets_num,lambda,asset_ret,asset_vol,asset_corr,sample_number,lower_bound,upper_bound){
  asset_ret=as.double(unlist(strsplit(asset_ret,',')))
  asset_vol=as.double(unlist(strsplit(asset_vol,',')))
  asset_var=asset_vol*asset_vol
  asset_corr=as.double(unlist(strsplit(asset_corr,',')))
  lower_bound=as.double(unlist(strsplit(lower_bound,',')))
  upper_bound=as.double(unlist(strsplit(upper_bound,',')))
  weights=scenario_no_cost_no_BN(assets_num,lambda,asset_ret,asset_var,asset_corr,sample_number,lower_bound,upper_bound)
  return(weights)
}
#################################build conditional probability matrix
build_cond<-function(affect_relation,prob_table){
  affect<-strsplit(affect_relation,";")[[1]]
  table<-strsplit(prob_table,";")[[1]]
  yn<-c("1","0")
  N=length(affect)
  cond_table=list()
  for(i in 1:N){
    bb<-paste0("~",affect[i])
    aa<-paste0("c(",table[i],")")
    cond_table[[i]]<-cptable(eval(parse(text=bb)),values=eval(parse(text=aa)),levels=yn)
  }
  return(cond_table)
}
####################basedon condition table build bayesian matrix
bayesian_matrix<-function(cond_table=NULL,asset_name){
  ###incase we don't have any cond_table initialize
  if(length(cond_table)==0){
    yn<-c("1","0")
    c<-cptable(~CLO,values=c(4,96),levels=yn)
    cm<-cptable(~CMBS|CLO,values = c(4,6,7,93),levels=yn)
    rm.cm<-cptable(~RMBS|CMBS,values=c(3,97,4,96),levels=yn)
    n=3
    cond_table=list(c,cm,rm.cm)
  }
  ####compile CPT table
  plist <- compileCPT(cond_table)
  ###########build net work
  net <- grain(plist)
  ###########query the net to get joint probability table
  bbb<-querygrain(net,nodes=asset_name, type="joint")
  ddd1<-data.frame(bbb)
  n=length(asset_name)
  pro_dict<-matrix(0,nrow=2^n,ncol=n+1)
  k<-1
  #convert joint probability table to formalized matrix 
  for(i in 1:ncol(ddd1)){
    for(j in 1:2){
      string<-paste0(row.names(ddd1)[j],".",substr(colnames(ddd1)[i],2,nchar(colnames(ddd1)[i])))
      value<-ddd1[j,i]
      pro_dict[k,1]=value
      pro_dict[k,2:(n+1)]=stringsplit(string)
      k=k+1 
    }
  }
  pro_dict<-data.frame(pro_dict)
  joint_table<-aggregate((pro_dict[,1]), as.list(pro_dict[,2:(n+1)]), FUN = sum)
  joint_table=t(joint_table[,c(n+1,seq(1,n))])
  colname=c()
  for(i in seq(0,2^n-1)){
    colname[i+1]<-paste0("p",i)
  }
  colnames(joint_table)<-colname
  rownames(joint_table)<-c("probability",asset_name)
  pro_dict<-data.frame(joint_table)
  pro_dict_all<<-pro_dict
  pro_dict<-pro_dict[,which(pro_dict[1,]!=0)]
  pro_dict<<-pro_dict
  return(pro_dict)
}
#####a help function for the transaction cost and finance cost to get the range of piecewise linear
get_number<-function(string){
  string=str_replace_all(string,">"," ")
  string=str_replace_all(string,"<"," ")
  string=str_replace_all(string,"="," ")
  string=str_replace_all(string,"&"," ")
  string=str_replace_all(string,"!"," ")
  string=strsplit(string," ")[[1]]
  string<-as.numeric(string[which(string!='')])
  string=na.omit(string)
  return(string)
}
###############################this is a help function to help get the accumulate cost for the trans and finance
get_accum_cost<-function(cost3,w_index,principal1,i,j,cond1){
  if(j==1){
    return(cost3[i,j+1]*w_index[i])
  }else{
    ret<-0
    for(kk in 1:(j-1)){
      range<-get_number(cond1[kk])
      if(length(range)==1) len<-(range[1]-0)/principal1
      else{
        len<-(range[2]-range[1])/principal1
      }
      ret<-ret+cost3[i,kk+1]*len
    }
    range<-get_number(cond1[j])
    ret<-ret+cost3[i,j+1]*(w_index[i]-range[1]/principal1)
    return(ret)
  }
}
####################calcualte finance cost and trans cost given a vector of weights
object_fun<-function(w_index,principal1,cond1,finan_cost3,leverage1){
  w<-w_index*principal1
  ww<-matrix(0,nrow=length(w),ncol=length(cond1))
  ret<-0
  for(i in 1:length(w)){
    for(j in 1:length(cond1)){
      if(w_index[i]==0) {
        next
      }else if(j==1||j==length(cond1)){
        ww[i,j]=paste0(w[i],cond1[j])
      }else{
        ww[i,j]=paste0(w[i],cond1[j],w[i])
      }
      if(eval(parse(text=ww[i,j]))){
        ret<-ret+get_accum_cost(finan_cost3,w_index,principal1,i,j,cond1)
      }
    }
  }
  return(-ret)
}
################################################display the cost in a foramt of vector
object_fun_vector<-function(w_index,principal1,cond1,finan_cost3,leverage1){
  w<-w_index*principal1
  ww<-matrix(0,nrow=length(w),ncol=length(cond1))
  ret<-rep(0,length(w_index))
  for(i in 1:length(w)){
    for(j in 1:length(cond1)){
      if(w_index[i]==0) {
        next
      }else if(j==1||j==length(cond1)){
        ww[i,j]=paste0(w[i],cond1[j])
      }else{
        ww[i,j]=paste0(w[i],cond1[j],w[i])
      }
      if(eval(parse(text=ww[i,j]))){
        ret[i]<-get_accum_cost(finan_cost3,w_index,principal1,i,j,cond1)
      }
    }
  }
  return(ret)
}
# equal_h0<<-function(w_index,leverage3=leverage1,principal2=principal1,con3=cond1,finan_cost3=finan_cost3){
#   return(sum(w_index)-leverage3)
# }
# equal_jac<<-function(x) nl.jacobian(x, equal_h0, heps = 1e-2)
#####################################################################################################################
######################################
######################################This function is used to calucalte the inequality constraints
cost_vector<-function(w_1,finan_cost,haircut,real_finance_weight,principal1=1720){
  ##cash position at MM
  #print(w_1)
  N<-length(w_1)
  N_fina_con<-ncol(finan_cost)
  if(length(N_fina_con)==0) {
    return(rep(0,N))
  }else if(is.na(w_1)[1]){
    return(rep(-1,N))
  }else{
    flag=0
    w=0
    list1=which(real_finance_weight!=1)
    ############first examine if have assets in short positions
    if(length(list1)!=0){
      w=sum(abs(((1-haircut)*real_finance_weight*w_1)[list1]))+sum(abs((haircut*w_1)[list1]))
      if(w_1[length(w_1)]<w){
        flag=1
        w_1[length(w_1)]=w_1[length(w_1)]-w
      }else{
        w_1[length(w_1)]=w_1[length(w_1)]-w
      }
    }
    i=0
    M=MAXIMUM_LOSS
    ####can't consider cash fee
    con2<-colnames(finan_cost)[2:N_fina_con]
    list=which(real_finance_weight!=1)
    # list=unique(c(list,which(real_finance_weight==0)))
    if(length(list)!=0){
      leverage1<<-sum(haircut[-list]*w_1[-list]+(1-haircut[-list])*real_finance_weight[-list]*w_1[-list])-(1-w)
    }else{
      leverage1<<-sum(haircut*w_1+(1-haircut)*real_finance_weight*w_1)-(1-w)
    }
    finan_cost2=finan_cost[,2:N_fina_con] 
    # finan_cost2[which(grepl("cash",finan_cost[,1])==TRUE),2:N_fina_con]=rep(MAXIMUM_LOSS,N_fina_con-1)
    haircut1=(1-haircut)*real_finance_weight*w_1
    w_index=seq(1,N)
    if(leverage1<=0){  #################don't need to leverage return
      return(rep(0,N))
    }else if(leverage1>round(sum(haircut1[intersect(which(w_1>0),which(real_finance_weight!=0))]),6)){   ########doesn't hold haircut constraints return 
      NN<-which(w_1!=0)[1]
      return(c(rep(0,NN-1),rep((sum(haircut1[intersect(which(w_1>0),which(real_finance_weight!=0))])-leverage1),N-NN+1))) 
    }else{
      ##########################convert piecewise linear to linear programming 
      mm<-length(w_index)
      if(length(con2)==1){ 
        finan_cost3<-finan_cost2
        lower_b<-rep(0,mm)
        upper_b<-ifelse(haircut1>0,haircut1,0)
        f.obj<-finan_cost3
        f.con<-rbind(diag(mm),diag(mm),rep(1,mm))
        f.dir<-c(rep("<=",mm),rep(">=",mm),"==")
        f.rhs<-c(upper_b,lower_b,leverage1)
        return(upper_b-lp("max", f.obj, f.con, f.dir, f.rhs)$solution)  ###convex piece wise linear
      }else if(length(con2)!=1){
        finan_cost3<--finan_cost2
        lower_b<-rep(0,mm)
        upper_b<-ifelse(haircut1>0,haircut1,0)
        f.obj<-c(rep(0,mm),rep(1,mm))
        f.con<-rbind(cbind(diag(mm),matrix(0,mm,mm)),cbind(diag(mm),matrix(0,mm,mm)),c(rep(1,mm),rep(0,mm)))
        for(kk in 1:mm){
          temp<-matrix(0,nrow=length(con2),ncol=2*mm)
          for(jj in 1:length(con2)){
            temp[jj,kk]<-finan_cost3[kk,jj]
            temp[jj,mm+kk]<--1
          }
          f.con<-rbind(f.con,temp)
        }
        f.dir<-c(rep("<=",mm),rep(">=",mm),"==",rep("<=",length(con2)*mm))
        f.rhs<-c(upper_b,lower_b,leverage1)
        for(kk in 1:mm){
          temp<-c(0)
          for(jj in 2:length(con2)){
            range<-get_number(con2[jj])
            if(length(range)==2){
              range1<-(range[2]-range[1])/principal1
            }else{
              range1<-range[1]/principal1
            } 
            temp<-c(temp,temp[length(temp)]+range1*(finan_cost3[kk,jj]-finan_cost3[kk,jj-1]))
          }
          f.rhs<-c(f.rhs,temp)
        }
        return(upper_b-lp("min", f.obj, f.con, f.dir, f.rhs)$solution[1:mm])   ###convex piece wise linear
      }
    }
  }
}
#####################################################
cost_vector2<-function(w_1,finan_cost,haircut,real_finance_weight,principal1=1720){
  ##cash position at MM
  #print(w_1)
  N<-length(w_1)
  N_fina_con<-ncol(finan_cost)
  if(length(N_fina_con)==0) {
    return(rep(0,N))
  }else if(is.na(w_1)[1]){
    return(rep(-1,N))
  }else{
    flag=0
    w=0
    list1=which(real_finance_weight!=1)
    ############first examine if have assets in short positions
    if(length(list1)!=0){
      w=sum(abs(((1-haircut)*real_finance_weight*w_1)[list1]))+sum(abs((haircut*w_1)[list1]))
      if(w_1[length(w_1)]<w){
        flag=1
        w_1[length(w_1)]=w_1[length(w_1)]-w
      }else{
        w_1[length(w_1)]=w_1[length(w_1)]-w
      }
    }
    i=0
    M=MAXIMUM_LOSS
    ####can't consider cash fee
    con2<-colnames(finan_cost)[2:N_fina_con]
    list=which(real_finance_weight!=1)
    # list=unique(c(list,which(real_finance_weight==0)))
    if(length(list)!=0){
      leverage1<<-sum(haircut[-list]*w_1[-list]+(1-haircut[-list])*real_finance_weight[-list]*w_1[-list])-(1-w)
    }else{
      leverage1<<-sum(haircut*w_1+(1-haircut)*real_finance_weight*w_1)-(1-w)
    }
    finan_cost2=finan_cost[,2:N_fina_con] 
    # finan_cost2[which(grepl("cash",finan_cost[,1])==TRUE),2:N_fina_con]=rep(MAXIMUM_LOSS,N_fina_con-1)
    haircut1=(1-haircut)*real_finance_weight*w_1
    w_index=seq(1,N)
    if(leverage1<=0){  #################don't need to leverage return
      return(rep(0,N))
    }else if(leverage1>round(sum(haircut1[intersect(which(w_1>0),which(real_finance_weight!=0))]),6)){   ########doesn't hold haircut constraints return 
      return(rep((sum(haircut1[intersect(which(w_1>0),which(real_finance_weight!=0))])-leverage1),N)) 
    }else{
      ##########################convert piecewise linear to linear programming 
      mm<-length(w_index)
      if(length(con2)==1){ 
        finan_cost3<-finan_cost2
        lower_b<-rep(0,mm)
        upper_b<-ifelse(haircut1>0,haircut1,0)
        f.obj<-finan_cost3
        f.con<-rbind(diag(mm),diag(mm),rep(1,mm))
        f.dir<-c(rep("<=",mm),rep(">=",mm),"==")
        f.rhs<-c(upper_b,lower_b,leverage1)
        return(lp("max", f.obj, f.con, f.dir, f.rhs)$solution)  ###convex piece wise linear
      }else if(length(con2)!=1){
        finan_cost3<--finan_cost2
        lower_b<-rep(0,mm)
        upper_b<-ifelse(haircut1>0,haircut1,0)
        f.obj<-c(rep(0,mm),rep(1,mm))
        f.con<-rbind(cbind(diag(mm),matrix(0,mm,mm)),cbind(diag(mm),matrix(0,mm,mm)),c(rep(1,mm),rep(0,mm)))
        for(kk in 1:mm){
          temp<-matrix(0,nrow=length(con2),ncol=2*mm)
          for(jj in 1:length(con2)){
            temp[jj,kk]<-finan_cost3[kk,jj]
            temp[jj,mm+kk]<--1
          }
          f.con<-rbind(f.con,temp)
        }
        f.dir<-c(rep("<=",mm),rep(">=",mm),"==",rep("<=",length(con2)*mm))
        f.rhs<-c(upper_b,lower_b,leverage1)
        for(kk in 1:mm){
          temp<-c(0)
          for(jj in 2:length(con2)){
            range<-get_number(con2[jj])
            if(length(range)==2){
              range1<-(range[2]-range[1])/principal1
            }else{
              range1<-range[1]/principal1
            } 
            temp<-c(temp,temp[length(temp)]+range1*(finan_cost3[kk,jj]-finan_cost3[kk,jj-1]))
          }
          f.rhs<-c(f.rhs,temp)
        }
        return(lp("min", f.obj, f.con, f.dir, f.rhs)$solution[1:mm])   ###convex piece wise linear
      }
    }
  }
}
####################################################this function returns the cost of our financing strategy and trans cost
cost<-function(w_now,w_1,trans_cost,finan_cost,haircut,real_finance_weight,principal1=1720){
  ##cash position at MM
  N<-length(w_now)
  N_tran_con<-ncol(trans_cost)
  N_fina_con<-ncol(finan_cost)
  if(length(N_tran_con)==0&&length(N_fina_con)==0) {
    return(0)
  }else if(is.na(w_1)[1]){
    return(MAXIMUM_LOSS)
  }else{                             #################################the first part calcualte transcation cost matrix
    con<-colnames(trans_cost)[2:N_tran_con]
    delta_w<-(w_1-w_now)*principal1
    ret_sum<-0
    for(j in 1:N){
      for(i in 1:(N_tran_con-1)){
        if(i==1||i==(N_tran_con-1)) {
          aa<-paste0(abs(delta_w[j]),con[i])
          if(eval(parse(text=aa))){
            ret_sum=ret_sum+get_accum_cost(trans_cost,abs(w_now-w_1),principal1,j,i,colnames(trans_cost)[-1]) ########add transaction cost to our total return
            
          }}
        else{
          aa<-paste0(abs(delta_w[j]),con[i],abs(delta_w[j]))
          if(eval(parse(text=aa))){
            ret_sum=ret_sum+get_accum_cost(trans_cost,abs(w_now-w_1),principal1,j,i,colnames(trans_cost)[-1]) ########add transaction cost to our total return
          }}
      }
    }
    ####################################################################calcualte finance strategy
    flag=0   ##########its a flag to see wheher we have enough cash to maintain margin
    w=0      ############################amount need to be subtract from the total amount
    list1=which(real_finance_weight!=1)      #################################find synthetic assets
   if(length(list1)!=0){        ##########if there exits assets synthetic
       w=sum(abs(((1-haircut)*real_finance_weight*w_1)[list1]))+sum(abs((haircut*w_1)[list1]))    #####################calcualte how much cash we should hold
         if(w_1[length(w_1)]<w){  ##############if we don;t have enough cash flag =1
           flag=1
           w_1[length(w_1)]=w_1[length(w_1)]-w     ##############cash's position be updated
         }else{                                ######else just update the cash's position
        w_1[length(w_1)]=w_1[length(w_1)]-w
   }
    }
      i=0
      ####initialize M
      M=MAXIMUM_LOSS
      ####can't consider cash fee
      con2<-colnames(finan_cost)[2:N_fina_con]
      list=which(real_finance_weight!=1)
      # list=unique(c(list,which(real_finance_weight==0)))
      if(length(list)!=0){
        leverage1<<-sum(haircut[-list]*w_1[-list]+(1-haircut[-list])*real_finance_weight[-list]*w_1[-list])-(1-w)
      }else{
        leverage1<<-sum(haircut*w_1+(1-haircut)*real_finance_weight*w_1)-(1-w)
      }
      finan_cost2=finan_cost[,2:N_fina_con] 
      haircut1=(1-haircut)*real_finance_weight*w_1
      w_index=seq(1,N)
      if(leverage1<=0){ #########dont need do finance just return
        return(ret_sum)
      }else if(leverage1>sum(haircut1[intersect(which(w_1>0),which(real_finance_weight!=0))])){ ######need return but can't hold haircut conditions return the largest cost(point is just make it continuous)
        return(ret_sum+min(finan_cost[,N_fina_con])*leverage1) 
      }else{
        w_index2<-w_index[which(haircut1>0)]
        mm<-length(w_index2)
        haircut1<-haircut1[w_index2]
        if(length(con2)==1){     ##################################only one condtition linear programming
        finan_cost3<-finan_cost2[w_index2]
        lower_b<-rep(0,mm)
        upper_b<-haircut1
        f.obj<-finan_cost3
        f.con<-rbind(diag(mm),diag(mm),rep(1,mm))
        f.dir<-c(rep("<=",mm),rep(">=",mm),"==")
        f.rhs<-c(upper_b,lower_b,leverage1)
        ret_sum<-ret_sum+lp("max", f.obj, f.con, f.dir, f.rhs)$objval   ###convex piece wise linear
        }else if(length(con2)!=1){     ############################# cobnvert the piecewise linear to linear programming
          finan_cost3<--finan_cost2[w_index2,]
          lower_b<-rep(0,mm)
          upper_b<-haircut1
          f.obj<-c(rep(0,mm),rep(1,mm))
          f.con<-rbind(cbind(diag(mm),matrix(0,mm,mm)),cbind(diag(mm),matrix(0,mm,mm)),c(rep(1,mm),rep(0,mm)))
         for(kk in 1:mm){
           temp<-matrix(0,nrow=length(con2),ncol=2*mm)
           for(jj in 1:length(con2)){
             temp[jj,kk]<-finan_cost3[kk,jj]
             temp[jj,mm+kk]<--1
           }
           f.con<-rbind(f.con,temp)
         }
          f.dir<-c(rep("<=",mm),rep(">=",mm),"==",rep("<=",length(con2)*mm))
          f.rhs<-c(upper_b,lower_b,leverage1)
          for(kk in 1:mm){
            temp<-c(0)
            for(jj in 2:length(con2)){
              range<-get_number(con2[jj])
              if(length(range)==2){
                range1<-(range[2]-range[1])/principal1
              }else{
                range1<-range[1]/principal1
              } 
              temp<-c(temp,temp[length(temp)]+range1*(finan_cost3[kk,jj]-finan_cost3[kk,jj-1]))
            }
            f.rhs<-c(f.rhs,temp)
          }
          ret_sum<-ret_sum-lp("min", f.obj, f.con, f.dir, f.rhs)$objval   ###convex piece wise linear
        }
        if(flag==1&&ret_sum!=0){
          return((ret_sum+min(finan_cost[,N_fina_con])*w))  
        }else{
          return(ret_sum)
        }
      }
  }
}
trans_cost_vector<-function(w_now,w_1,trans_cost,principal1){
  N<-length(w_now)
  cost<-rep(0,N)
  N_tran_con<-ncol(trans_cost)
  con<-colnames(trans_cost)[2:N_tran_con]
  delta_w<-(w_1-w_now)*principal1
  for(j in 1:N){
    for(i in 1:(N_tran_con-1)){
      if(i==1||i==(N_tran_con-1)) {
        aa<-paste0(abs(delta_w[j]),con[i])
        if(eval(parse(text=aa))){
          cost[j]=cost[j]+get_accum_cost(trans_cost,abs(w_now-w_1),principal1,j,i,colnames(trans_cost)[-1])
          
        }}
      else{
        aa<-paste0(abs(delta_w[j]),con[i],abs(delta_w[j]))
        if(eval(parse(text=aa))){
          cost[j]=cost[j]+get_accum_cost(trans_cost,abs(w_now-w_1),principal1,j,i,colnames(trans_cost)[-1])
        }}
    }
  }
  return(cost)
}
finance_cost_vector<-function(w_1,finan_cost,principal1){
  N<-length(w_1)
  cost<-rep(0,N)
  N_finan_con<-ncol(finan_cost)
  con<-colnames(finan_cost)[2:N_finan_con]
  delta_w<-w_1*principal1
  for(j in 1:N){
    for(i in 1:(N_finan_con-1)){
      if(i==1||i==(N_finan_con-1)) {
        aa<-paste0(abs(delta_w[j]),con[i])
        if(eval(parse(text=aa))){
          cost[j]=cost[j]+get_accum_cost(finan_cost,w_1,principal1,j,i,colnames(finan_cost)[-1])
          
        }}
      else{
        aa<-paste0(abs(delta_w[j]),con[i],abs(delta_w[j]))
        if(eval(parse(text=aa))){
          cost[j]=cost[j]+get_accum_cost(finan_cost,w_1,principal1,j,i,colnames(finan_cost)[-1])
        }}
    }
  }
  return(cost)
}
####################################################this function returns the vector of our financing strategy
cost2<-function(w_now,w_1,trans_cost,finan_cost,haircut,real_finance_weight,principal1=1720){
  ##cash position at MM
  N<-length(w_now)
  mat<-matrix(0,nrow=N,ncol=7)
  mat[,4]=haircut
  mat[,5]=w_1
  mat[,6]=w_now
  mat[,7]=w_1-w_now
  N_tran_con<-ncol(trans_cost)
  N_fina_con<-ncol(finan_cost)
  if(length(N_tran_con)==0&&length(N_fina_con)==0) {
    return(mat)
  }else{
    con<-colnames(trans_cost)[2:N_tran_con]
    delta_w<-(w_1-w_now)*principal1
    ret_sum<-0
    for(j in 1:N){
      for(i in 1:(N_tran_con-1)){
        if(i==1||i==(N_tran_con-1)) {
          aa<-paste0(abs(delta_w[j]),con[i])
          if(eval(parse(text=aa))){
            mat[j,1]=mat[j,1]+get_accum_cost(trans_cost,abs(w_now-w_1),principal1,j,i,colnames(trans_cost)[-1])
            
          }}
        else{
          aa<-paste0(abs(delta_w[j]),con[i],abs(delta_w[j]))
          if(eval(parse(text=aa))){
            mat[j,1]=mat[j,1]+get_accum_cost(trans_cost,abs(w_now-w_1),principal1,j,i,colnames(trans_cost)[-1])
          }}
      }
    }
    list1=which(real_finance_weight!=1)
    if(length(list1)!=0){
      w=sum(abs(((1-haircut)*real_finance_weight*w_1)[list1]))+sum(abs((haircut*w_1)[list1]))
    }else{
      w=0
    }
    w_1[length(w_1)]=w_1[length(w_1)]-w
    i=0
    M=MAXIMUM_LOSS
    ####can't consider cash fee
    haircut1=(1-haircut)*real_finance_weight*w_1
    con2<-colnames(finan_cost)[2:N_fina_con]
    list=which(real_finance_weight!=1)
    # list=unique(c(list,which(real_finance_weight==0)))
    if(length(list)!=0){
      leverage1<<-sum(haircut[-list]*w_1[-list]+(1-haircut[-list])*real_finance_weight[-list]*w_1[-list])-(1-w)
    }else{
      leverage1<<-sum(haircut*w_1+(1-haircut)*real_finance_weight*w_1)-(1-w)
    }
    finan_cost2=finan_cost
    # finan_cost2[which(grepl("cash",finan_cost[,1])==TRUE),2:N_fina_con]=rep(MAXIMUM_LOSS,N_fina_con-1)
    w_index=seq(1,N)
    haircut1<-haircut1
    finan_cost3<-finan_cost2[,2:N_fina_con]
    w_init<-rep(0,length(haircut1))
    lower_b<-rep(0,length(haircut1))
    upper_b<-ifelse(haircut1>0,haircut1,0)
    if(leverage1<=0){
      return(mat)
    }else{
      mm<-length(haircut1)
      haircut1<-haircut1
      if(length(con2)==1){ 
        lower_b<-rep(0,mm)
        upper_b<-ifelse(haircut1>0,haircut1,0)
        f.obj<-finan_cost3
        f.con<-rbind(diag(mm),diag(mm),rep(1,mm))
        f.dir<-c(rep("<=",mm),rep(">=",mm),"==")
        f.rhs<-c(upper_b,lower_b,leverage1)
        ress<-lp("max", f.obj, f.con, f.dir, f.rhs)$solution   ###convex piece wise linear
        mat[,3]<-ress
        mat[,2]=object_fun_vector(ress,principal1,con2,finan_cost2,leverage1)
        return(mat)
      }else if(length(con2)!=1){
        finan_cost3<--finan_cost3
        lower_b<-rep(0,mm)
        upper_b<-ifelse(haircut1>0,haircut1,0)
        f.obj<-c(rep(0,mm),rep(1,mm))
        f.con<-rbind(cbind(diag(mm),matrix(0,mm,mm)),cbind(diag(mm),matrix(0,mm,mm)),c(rep(1,mm),rep(0,mm)))
        for(kk in 1:mm){
          temp<-matrix(0,nrow=length(con2),ncol=2*mm)
          for(jj in 1:length(con2)){
            temp[jj,kk]<-finan_cost3[kk,jj]
            temp[jj,mm+kk]<--1
          }
          f.con<-rbind(f.con,temp)
        }
        f.dir<-c(rep("<=",mm),rep(">=",mm),"==",rep("<=",length(con2)*mm))
        f.rhs<-c(upper_b,lower_b,leverage1)
        for(kk in 1:mm){
          temp<-c(0)
          for(jj in 2:length(con2)){
            range<-get_number(con2[jj])
            if(length(range)==2){
              range1<-(range[2]-range[1])/principal1
            }else{
              range1<-range[1]/principal1
            } 
            temp<-c(temp,temp[length(temp)]+range1*(finan_cost3[kk,jj]-finan_cost3[kk,jj-1]))
          }
          f.rhs<-c(f.rhs,temp)
        }
        ress<-lp("min", f.obj, f.con, f.dir, f.rhs)  ###convex piece wise linear
      mat[,3]<-ress$solution[1:mm]
      mat[,2]=-ress$solution[(mm+1):(2*mm)]
      return(mat)
    }
    }
  }
}
#################################################utility functions
power_utility<-function(x,beta){
u<-ifelse(x<=MINIMUM_RET,1/(1-beta)*(MINIMUM_RET^(1-beta)-1)+MINIMUM_RET^(-beta)*(x-MINIMUM_RET),1/(1-beta)*(x^(1-beta)-1))
  return(u)
}
log_utility<-function(x,beta){
  u<-ifelse(x<=MINIMUM_RET,(log(MINIMUM_RET)+1/MINIMUM_RET*(x-MINIMUM_RET)),log(x))
  return(u)
}
expo_utility<-function(x,lambda){
  u<- (1-exp(-lambda*x))/lambda
  return(u)
}
##############################################call utility functions
get_expo_ut<-function(ret,beta,w_now,w1,tcost){
  return(expo_utility(ret+tcost,beta))
}
get_power_ut<-function(ret,beta,w_now,w1,tcost){
  return(power_utility(1+ret+tcost,beta))
}
get_log_ut<-function(ret,beta,w_now,w1,tcost){
  return(log_utility(1+beta*(ret+tcost),beta))
}
####################################extreme part utility value just matrix multiply
get_extreme_expo_uti<<-function(pro_dict,loss,beta,w_now,w,tcost,mu,rand_sub){
  rett<-0
  u<-rep(1,nrow(pro_dict)-1)
  pro_dict<-as.matrix(pro_dict)
  rett<-sum(1/nrow(rand_sub)*pro_dict[1,2:ncol(pro_dict)]*expo_utility((rep(t(pro_dict)[2:ncol(pro_dict),2:nrow(pro_dict)]%*%(loss*w),nrow(rand_sub))+as.vector((u-t(pro_dict)[2:ncol(pro_dict),2:nrow(pro_dict)])%*%(apply(rand_sub,1,function(x) x*w))))+tcost,beta))
  return(rett)
}
get_extreme_power_uti<<-function(pro_dict,loss,beta,w_now,w,tcost,mu,rand_sub){
  rett<-0
  u<-rep(1,nrow(pro_dict)-1)
  pro_dict<-as.matrix(pro_dict)
  rett<-sum(1/nrow(rand_sub)*pro_dict[1,2:ncol(pro_dict)]*power_utility((1+rep(t(pro_dict)[2:ncol(pro_dict),2:nrow(pro_dict)]%*%(loss*w), nrow(rand_sub))+as.vector((u-t(pro_dict)[2:ncol(pro_dict),2:nrow(pro_dict)])%*%(apply(rand_sub,1,function(x) x*w))))+tcost,beta))
  return(rett)
}
#####
get_extreme_log_uti<<-function(pro_dict,loss,beta,w_now,w,tcost,mu,rand_sub){
  rett<-0
  u<-rep(1,nrow(pro_dict)-1)
  pro_dict<-as.matrix(pro_dict)
  rett<-sum(1/nrow(rand_sub)*pro_dict[1,2:ncol(pro_dict)]*log_utility(1+beta*(rep(t(pro_dict)[2:ncol(pro_dict),2:nrow(pro_dict)]%*%(loss*w),nrow(rand_sub))+as.vector((u-t(pro_dict)[2:ncol(pro_dict),2:nrow(pro_dict)])%*%(apply(rand_sub,1,function(x) x*w)))+tcost),beta))
  return(rett)
}
#########################################################total functions to combine normal part and extreme part
power_find_w<-function(w1,w_now=w_now1,beta1=beta11,trans_cost=trans_cost1,finan_cost=finan_cost1,haircut=haircut2,real_finance_weight=real_finance_weight1,principal1=principal11,rand2=rand21,rand_sub=rand_sub1,rand2_moment=rand21_moment,loss1=loss11,pro_dict=pro_dictt,k=k1,mu=mu1){
  tcost<-cost(w_now,w1,trans_cost,finan_cost,haircut,real_finance_weight,principal1)
  exp_pow_ut<-(1-k)*mean(na.omit(get_power_ut(rand2 %*% w1,beta1,w_now,w1,tcost)))+(k/(1-pro_dict$p0[1]))*get_extreme_power_uti(pro_dict,loss1,beta1,w_now,w1,tcost,mu,rand_sub)
  return(-exp_pow_ut)
}
####gradient function
power_find_w_grad<-function(w1,w_now=w_now1,beta1=beta11,trans_cost=trans_cost1,finan_cost=finan_cost1,haircut=haircut2,real_finance_weight=real_finance_weight1,principal1=principal11,rand2=rand21,rand_sub=rand_sub1,rand2_moment=rand21_moment,loss1=loss11,pro_dict=pro_dictt,k=k1,mu=mu1)  nl.grad(w1,power_find_w,heps=1e-10)

log_find_w<-function(w1,w_now=w_now1,beta1=beta11,trans_cost=trans_cost1,finan_cost=finan_cost1,haircut=haircut2,real_finance_weight=real_finance_weight1,principal1=principal11,rand2=rand21,rand_sub=rand_sub1,rand2_moment=rand21_moment,loss1=loss11,pro_dict=pro_dictt,k=k1,mu=mu1){
  tcost<-cost(w_now,w1,trans_cost,finan_cost,haircut,real_finance_weight,principal1)
  exp_log_ut<-(1-k)*mean(na.omit(get_log_ut(rand2 %*% w1,beta1,w_now,w1,tcost)))+k/(1-pro_dict$p0[1])*get_extreme_log_uti(pro_dict,loss1,beta1,w_now,w1,tcost,mu,rand_sub)
  return(-exp_log_ut)  
}
log_find_w_grad<-function(w1,w_now=w_now1,beta1=beta11,trans_cost=trans_cost1,finan_cost=finan_cost1,haircut=haircut2,real_finance_weight=real_finance_weight1,principal1=principal11,rand2=rand21,rand_sub=rand_sub1,rand2_moment=rand21_moment,loss1=loss11,pro_dict=pro_dictt,k=k1,mu=mu1)  nl.grad(w1,log_find_w,heps=1e-10)

expo_find_w<-function(w1,w_now=w_now1,beta1=beta11,trans_cost=trans_cost1,finan_cost=finan_cost1,haircut=haircut2,real_finance_weight=real_finance_weight1,principal1=principal11,rand2=rand21,rand_sub=rand_sub1,rand2_moment=rand21_moment,loss1=loss11,pro_dict=pro_dictt,k=k1,mu=mu1){
  tcost<-cost(w_now,w1,trans_cost,finan_cost,haircut,real_finance_weight,principal1)
  exp_expo_ut<-(1-k)*mean(na.omit(get_expo_ut(rand2 %*% w1,beta1,w_now,w1,tcost)))+k/(1-pro_dict$p0[1])*get_extreme_expo_uti(pro_dict,loss1,beta1,w_now,w1,tcost,mu,rand_sub)
  return(-exp_expo_ut)
}
expo_find_w_grad<-function(w1,w_now=w_now1,beta1=beta11,trans_cost=trans_cost1,finan_cost=finan_cost1,haircut=haircut2,real_finance_weight=real_finance_weight1,principal1=principal11,rand2=rand21,rand_sub=rand_sub1,rand2_moment=rand21_moment,loss1=loss11,pro_dict=pro_dictt,k=k1,mu=mu1)  nl.grad(w1,expo_find_w,heps=1e-10)
function_expo<-function(x){
  return(exp(-x*mid1)-prob1*exp(-x*min1)-(1-prob1)*exp(-x*max1))
}
function_power<-function(x){
  return((1+mid1)^(1-x)-prob1*(1+min1)^(1-x)-(1-prob1)*(1+max1)^(1-x))
}
function_log<-function(x){
  return(log(1+x*mid1)-prob1*log(1+x*min1)-(1-prob1)*log(1+x*max1))
}
calculate_l<-function(method,mid,min,max,prob){
  mid1<<-mid
  min1<<-min
  max1<<-max
  prob1<<-prob
  
  if(method=='log'){
    a<-uniroot.all(function_log,c(0,30))
    return(a[which(a>0)])
  }else if(method=='expo'){
   a<-uniroot.all(function_expo,c(0,30))
   return(a[which(a>0)])
  }else{
    a<-uniroot.all(function_power,c(0,30))
    return(a[which(a>1.00001)])
  }
}
stringsplit<-function(string){
  b<-seq(1,2)
  a<-strsplit(string,NULL)[[1]]
  ll<-(length(a)+1)/2
  for(i in 1:ll){
    b[i]<-as.numeric(a[2*i-1])
  }
  return(b)
}
generate_moment_matching<-function(pro_dict1,mu1,k,asset_cov,N_sample){
  p0<-pro_dict1$p0[1]
  pro_dict_matrix<-as.matrix(pro_dict1[2:nrow(pro_dict1),])
  pro_dict1<-as.matrix(pro_dict1)
  mu_adjust<-mu1
  P<-mu1
  P1<-mu1
  P2<-mu1
  var_T<-mu1
  var_s<-mu1
  for(idx in 1:length(mu1)){
    mu_adjust[idx]<-(1-k+k*(1-sum(pro_dict1[1,which(pro_dict1[idx+1,]==1)])-p0)/(1-p0))*mu1[idx]+k*sum(pro_dict1[1,which(pro_dict1[idx+1,]==1)])/(1-p0)*loss1[idx]
    P[idx]<-sum(pro_dict1[1,which(pro_dict1[idx+1,]==1)])
    P1[idx]<-k-k/(1-p0)*(1-P[idx]-p0)
    var_s[idx]<-P[idx]*(1-P[idx])
    P2[idx]<-k/(1-p0)*P[idx]
    var_1<-P1[idx]*(1-P1[idx])
    var_2<-P2[idx]*(1-P2[idx])  ##(k-k/(1-p0)*(P[idx]-p0))*(1-(k-k/(1-p0)*(P[idx]-p0)))
    var_T[idx]=(1-P1[idx])^2*asset_cov[idx,idx]+P1[idx]*(1-P1[idx])*(asset_cov[idx,idx]+(mu1[idx])^2)+var_2*(loss1[idx])^2+P2[idx]^2*0-2*(1-P1[idx])*P2[idx]*mu1[idx]*(loss1[idx])
     var_T[idx]=ifelse(var_T[idx]<0,0,var_T[idx])
    }
  #covariance matrix
  cov_adjust<-matrix(0,nrow=length(mu1),ncol=length(mu1))
  E_matrix<-matrix(0,nrow=length(mu1),ncol=length(mu1))
  for(i in 1:length(mu1)){
    cov_adjust[i,i]=var_T[i]
  }
  P_00<-matrix(0,nrow=length(mu1),ncol=length(mu1))
  P_10<-matrix(0,nrow=length(mu1),ncol=length(mu1))
  P_01<-matrix(0,nrow=length(mu1),ncol=length(mu1))
  P_11<-matrix(0,nrow=length(mu1),ncol=length(mu1))
  P1_00<-matrix(0,nrow=length(mu1),ncol=length(mu1))
  P1_10<-matrix(0,nrow=length(mu1),ncol=length(mu1))
  P1_01<-matrix(0,nrow=length(mu1),ncol=length(mu1))
  P1_11<-matrix(0,nrow=length(mu1),ncol=length(mu1))
  P2_00<-matrix(0,nrow=length(mu1),ncol=length(mu1))
  P2_10<-matrix(0,nrow=length(mu1),ncol=length(mu1))
  P2_01<-matrix(0,nrow=length(mu1),ncol=length(mu1))
  P2_11<-matrix(0,nrow=length(mu1),ncol=length(mu1))
  ######initialize matrix
  for(idx1 in 1:length(mu1)){
    for(idx2 in 1:length(mu1)){
      if(idx1!=idx2){
        M_00<-intersect(which(pro_dict1[(idx1+1),]==0),which(pro_dict1[(idx2+1),]==0))
        p_ij_00=ifelse(length(M_00)==0,0,sum(pro_dict1[1,M_00]))
        M_01<-intersect(which(pro_dict1[(idx1+1),]==0),which(pro_dict1[(idx2+1),]==1))
        p_ij_01=ifelse(length(M_01)==0,0,sum(pro_dict1[1,M_01]))
        M_10<-intersect(which(pro_dict1[(idx1+1),]==1),which(pro_dict1[(idx2+1),]==0))
        p_ij_10=ifelse(length(M_10)==0,0,sum(pro_dict1[1,M_10]))
        M_11<-intersect(which(pro_dict1[(idx1+1),]==1),which(pro_dict1[(idx2+1),]==1))
        p_ij_11=ifelse(length(M_11)==0,0,sum(pro_dict1[1,M_11]))
        ##calculate e_ij
        e_ij=asset_cov[idx1,idx2]+mu1[idx1]*mu1[idx2]
        E_matrix[idx1,idx2]=(1-k+k/(1-p0)*(p_ij_00-p0))*e_ij+k/(1-p0)*p_ij_01*mu1[idx1]*loss1[idx2]+k/(1-p0)*p_ij_10*mu1[idx2]*loss1[idx1]+k/(1-p0)*p_ij_11*loss1[idx1]*loss1[idx2] 
        P_00[idx1,idx2]=p_ij_00
        P_10[idx1,idx2]=p_ij_10
        P_01[idx1,idx2]=p_ij_01
        P_11[idx1,idx2]=p_ij_11
        P1_00[idx1,idx2]=k-k/(1-p0)*(1-p_ij_00-p0)
        P1_10[idx1,idx2]=k-k/(1-p0)*(1-p_ij_10-p0)
        P1_01[idx1,idx2]=k-k/(1-p0)*(1-p_ij_01-p0)
        P1_11[idx1,idx2]=k-k/(1-p0)*(1-p_ij_11-p0)
        P2_00[idx1,idx2]=k/(1-p0)*p_ij_00
        P2_10[idx1,idx2]=k/(1-p0)*p_ij_10
        P2_01[idx1,idx2]=k/(1-p0)*p_ij_01
        P2_11[idx1,idx2]=k/(1-p0)*p_ij_11
      }
    }
  }
  for(idx1 in 1:length(mu1)){
    for(idx2 in 1:length(mu1)){
      if(idx1!=idx2){
        e_ij=asset_cov[idx1,idx2]+mu1[idx1]*mu1[idx2]   #
        cov_adjust[idx1,idx2]=(1-k+k/(1-p0)*(P_00[idx1,idx2]-p0))*e_ij+P2_01[idx1,idx2]*mu1[idx1]*loss1[idx2]+P2_10[idx1,idx2]*loss1[idx1]*mu1[idx2]+P2_11[idx1,idx2]*loss1[idx1]*loss1[idx2]-(1-P1[idx1])*(1-P1[idx2])*mu1[idx1]*mu1[idx2]-P2[idx1]*P2[idx2]*loss1[idx1]*loss1[idx2]-(1-P1[idx1])*P2[idx2]*mu1[idx1]*loss1[idx2]-(1-P1[idx2])*P2[idx1]*mu1[idx2]*loss1[idx1] 
      }
    }
  }
  ###########################generate_copula
  N<-length(mu1)
  asset_corr_adjust<-rep(0,N*(N-1)/2)
  idx=1
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      if(cov_adjust[i,i]!=0&cov_adjust[j,j]!=0){
        asset_corr_adjust[idx]=cov_adjust[i,j]/sqrt(cov_adjust[i,i]*cov_adjust[j,j])
        # asset_corr_adjust[idx]=ifelse((asset_corr_adjust[idx]<1)&(asset_corr_adjust[idx]>-1),asset_corr_adjust[idx],ifelse(asset_corr_adjust[idx]>1,1,-1))
      }else{
        asset_corr_adjust[idx]=0
      }
      idx=idx+1
    }
  }
  paramMargins_r_adjust <- list()
  # Now the new experiments
  for(i in 1:N){
    paramMargins_r_adjust[[length(paramMargins_r_adjust)+1]] <- list(mean =mu_adjust[i], sd =sqrt(cov_adjust[i,i]))
  }
  if(length(which(var_T!=0))==0){
    rand2_moment<-t(matrix(rep(mu_adjust,N_sample),nrow=length(mu_adjust)))
  }else{
    margins_r<-rep("norm",N)
    rand2_moment<-generate_N_rand(N,asset_corr_adjust,margins_r,paramMargins_r_adjust,N_sample)
  }
  return(rand2_moment)
}
expo_find_reduce<-function(w1,w_now=w_now1,beta1=beta11,trans_cost=trans_cost1,finan_cost=finan_cost1,haircut=haircut2,real_finance_weight=real_finance_weight1,principal1=principal11,rand2=rand21,rand_sub=rand_sub1,rand2_moment=rand21_moment,loss1=loss11,pro_dict=pro_dictt,k=k1,mu=mu1){
  tcost<-cost(w_now,w1,trans_cost,finan_cost,haircut,real_finance_weight,principal1)
  exp_expo_ut<-mean(get_expo_ut(rand2_moment %*% w1,beta1,w_now,w1,tcost))
  return(-(exp_expo_ut))
}
log_find_reduce<-function(w1,w_now=w_now1,beta1=beta11,trans_cost=trans_cost1,finan_cost=finan_cost1,haircut=haircut2,real_finance_weight=real_finance_weight1,principal1=principal11,rand2=rand21,rand_sub=rand_sub1,rand2_moment=rand21_moment,loss1=loss11,pro_dict=pro_dictt,k=k1,mu=mu1){
  tcost<-cost(w_now,w1,trans_cost,finan_cost,haircut,real_finance_weight,principal1)
  exp_expo_ut<-mean(get_log_ut(rand2_moment %*% w1,beta1,w_now,w1,tcost))
  return(-(exp_expo_ut))
}
power_find_reduce<-function(w1,w_now=w_now1,beta1=beta11,trans_cost=trans_cost1,finan_cost=finan_cost1,haircut=haircut2,real_finance_weight=real_finance_weight1,principal1=principal11,rand2=rand21,rand_sub=rand_sub1,rand2_moment=rand21_moment,loss1=loss11,pro_dict=pro_dictt,k=k1,mu=mu1){
  tcost<-cost(w_now,w1,trans_cost,finan_cost,haircut,real_finance_weight,principal1)
  exp_expo_ut<-mean(get_power_ut(rand2_moment %*% w1,beta1,w_now,w1,tcost))
  return(-(exp_expo_ut))
}

combo_find_reduce<-function(w1,w_now=w_now1,x_downturning=x_downturning1,AA=AA1,ll=ll1,k_2=k_21,k_1=k_11,x_1=x_11,trans_cost=trans_cost1,finan_cost=finan_cost1,haircut=haircut2,real_finance_weight=real_finance_weight1,principal1=principal11,rand2=rand21,rand_sub=rand_sub1,rand2_moment=rand21_moment,loss1=loss11,pro_dict=pro_dictt,k=k1,mu=mu1){
  tcost<-cost(w_now,w1,trans_cost,finan_cost,haircut,real_finance_weight,principal1)
  exp_pow_ut<-mean((get_combo_ut(rand2_moment %*% w1,x_downturning,AA,ll,k_2,k_1,x_1,w_now,w1,tcost)))
  return(-exp_pow_ut)
}
expo_find_reduce_grad<-function(w1,w_now=w_now1,beta1=beta11,trans_cost=trans_cost1,finan_cost=finan_cost1,haircut=haircut2,real_finance_weight=real_finance_weight1,principal1=principal11,rand2=rand21,rand_sub=rand_sub1,rand2_moment=rand21_moment,loss1=loss11,pro_dict=pro_dictt,k=k1,mu=mu1)  nl.grad(w1,expo_find_reduce,heps=1e-10)
log_find_reduce_grad<-function(w1,w_now=w_now1,beta1=beta11,trans_cost=trans_cost1,finan_cost=finan_cost1,haircut=haircut2,real_finance_weight=real_finance_weight1,principal1=principal11,rand2=rand21,rand_sub=rand_sub1,rand2_moment=rand21_moment,loss1=loss11,pro_dict=pro_dictt,k=k1,mu=mu1)  nl.grad(w1,log_find_reduce,heps=1e-10)
power_find_reduce_grad<-function(w1,w_now=w_now1,beta1=beta11,trans_cost=trans_cost1,finan_cost=finan_cost1,haircut=haircut2,real_finance_weight=real_finance_weight1,principal1=principal11,rand2=rand21,rand_sub=rand_sub1,rand2_moment=rand21_moment,loss1=loss11,pro_dict=pro_dictt,k=k1,mu=mu1)  nl.grad(w1,power_find_reduce,heps=1e-10)
combo_find_reduce_grad<-function(w1,w_now=w_now1,x_downturning=x_downturning1,AA=AA1,ll=ll1,k_2=k_21,k_1=k_11,x_1=x_11,trans_cost=trans_cost1,finan_cost=finan_cost1,haircut=haircut2,real_finance_weight=real_finance_weight1,principal1=principal11,rand2=rand21,rand_sub=rand_sub1,rand2_moment=rand21_moment,loss1=loss11,pro_dict=pro_dictt,k=k1,mu=mu1) nl.grad(w1,combo_find_reduce,heps=1e-10)

initial_point<-function(method,moment,start,lower_bound,upper_bound,w_now,lambda,trans_cost,finan_cost,haircut,real_finance_weight,principal1,rand2,rand_sub,rand2_moment,loss1,pro_dict1,k,mu,signal){
  if(moment==TRUE){
    globalsearch<<-(-1)
   if(method=='power'){
     res3<<-nloptr(x0=w_now,eval_f = power_find_reduce,eval_grad_f = power_find_reduce_grad,eval_g_ineq = eval_g0,eval_jac_g_ineq = eval_g0_jac,
                eval_g_eq = eval_h0,eval_jac_g_eq = eval_h0_jac,lb=lower_bound,ub=upper_bound,opts=opts_moment,w_now=w_now,beta1=lambda,
                trans_cost=trans_cost,finan_cost=finan_cost,haircut=haircut,real_finance_weight=real_finance_weight,principal1=principal1,rand2=rand2,rand_sub=rand_sub,rand2_moment=rand2_moment,loss1=loss1,pro_dict=pro_dict1,k=k,mu=mu)
   }else if(method=='log'){
     res3<<-nloptr(x0=w_now,eval_f = log_find_reduce,eval_grad_f = log_find_reduce_grad,eval_g_ineq = eval_g0,eval_jac_g_ineq = eval_g0_jac,
                eval_g_eq = eval_h0,eval_jac_g_eq = eval_h0_jac,lb=lower_bound,ub=upper_bound,opts=opts_moment,w_now=w_now,beta1=lambda,
                trans_cost=trans_cost,finan_cost=finan_cost,haircut=haircut,real_finance_weight=real_finance_weight,principal1=principal1,rand2=rand2,rand_sub=rand_sub,rand2_moment=rand2_moment,loss1=loss1,pro_dict=pro_dict1,k=k,mu=mu)
     
   }else if(method=='expo'){
     res3<<-nloptr(x0=w_now,eval_f = expo_find_reduce,eval_grad_f = expo_find_reduce_grad,eval_g_ineq = eval_g0,eval_jac_g_ineq = eval_g0_jac,
                eval_g_eq = eval_h0,eval_jac_g_eq = eval_h0_jac,lb=lower_bound,ub=upper_bound,opts=opts_moment,w_now=w_now,beta1=lambda,
                trans_cost=trans_cost,finan_cost=finan_cost,haircut=haircut,real_finance_weight=real_finance_weight,principal1=principal1,rand2=rand2,rand_sub=rand_sub,rand2_moment=rand2_moment,loss1=loss1,pro_dict=pro_dict1,k=k,mu=mu)
   }
    initial_point<-res3$solution
  }else{
    if(method=='power'){
      ##########################################first global search based on start point
      buffer<<-capture.output(res3<<-nloptr(x0=start,eval_f = power_find_w,eval_grad_f = power_find_w_grad,eval_g_ineq = eval_g0,eval_jac_g_ineq = eval_g0_jac,
                                            eval_g_eq = eval_h0,eval_jac_g_eq = eval_h0_jac,lb=lower_bound,ub=upper_bound,opts=opts_global,w_now=w_now,beta1=lambda,
                                            trans_cost=trans_cost,finan_cost=finan_cost,haircut=haircut,real_finance_weight=real_finance_weight,principal1=principal1,rand2=rand2,rand_sub=rand_sub,rand2_moment=rand2_moment,loss1=loss1,pro_dict=pro_dict1,k=k,mu=mu))
      ################if res3 return is -1 which means fail or solution is equal to start point which means optimize fails use the best point we can find
      ##############same idea for log and exponential utility 
      if(sum(ifelse(res3$solution==start,0,1))==0||res3$status==-1){
        globalsearch<<-1 
        initial_point<-find_optimal(buffer,signal) 
      }else{
        globalsearch<<-0
        initial_point<-res3$solution
      }}else if(method=='log'){
        buffer<<-capture.output(res3<<-nloptr(x0=start,eval_f = log_find_w,eval_grad_f = log_find_w_grad,eval_g_ineq = eval_g0,eval_jac_g_ineq = eval_g0_jac,
                                               eval_g_eq = eval_h0,eval_jac_g_eq = eval_h0_jac,lb=lower_bound,ub=upper_bound,opts=opts_global,w_now=w_now,beta1=lambda,
                                               trans_cost=trans_cost,finan_cost=finan_cost,haircut=haircut,real_finance_weight=real_finance_weight,principal1=principal1,rand2=rand2,rand_sub=rand_sub,rand2_moment=rand2_moment,loss1=loss1,pro_dict=pro_dict1,k=k,mu=mu))
        if(sum(ifelse(res3$solution==start,0,1))==0||res3$status==-1){
          globalsearch<<-1
          initial_point<-find_optimal(buffer,signal) 
        }else{
          globalsearch<<-0
          initial_point<-res3$solution
        }}else if(method=='expo'){
          buffer<<-capture.output(res3<<-nloptr(x0=start,eval_f = expo_find_w,eval_grad_f = expo_find_w_grad,eval_g_ineq = eval_g0,eval_jac_g_ineq = eval_g0_jac,
                                                eval_g_eq = eval_h0,eval_jac_g_eq = eval_h0_jac,lb=lower_bound,ub=upper_bound,opts=opts_global,w_now=w_now,beta1=lambda,
                                                trans_cost=trans_cost,finan_cost=finan_cost,haircut=haircut,real_finance_weight=real_finance_weight,principal1=principal1,rand2=rand2,rand_sub=rand_sub,rand2_moment=rand2_moment,loss1=loss1,pro_dict=pro_dict1,k=k,mu=mu))
          if(sum(ifelse(res3$solution==start,0,1))==0||res3$status==-1){
          globalsearch<<-1
          initial_point<-find_optimal(buffer,signal) 
          }else{
            globalsearch<<-0 
            initial_point<-res3$solution
          }}}
  return(initial_point)
}
scenario_cost_BN<-function(method,assets_num,lambda,asset_ret,asset_var,asset_corr,sample_number,extreme_stress_loss,pro_dict,principal1,
                           trans_cost,finan_cost,haircut,real_finance_weight,w_now,lower_bound,upper_bound,subjective_k,
                           inequality_constraints,equality_constraints,maxeval_global,maxeval_local,asset_name1=NULL,moment,smallprob){
  ##trans_cost matrix like
  N=assets_num
  A=lambda
  mu<<-asset_ret
  asset_cor<-matrix(1,N,N)
  kk=1
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      asset_cor[i,j]=asset_corr[kk]
      kk=kk+1
      asset_cor[j,i]=asset_cor[i,j]
    }
  }
  asset_cov<<-r2cov(sd =sqrt(asset_var),R = asset_cor)
  if(length(asset_name1)==0){
    asset_name1=rep("asset_",N)
    for(i in 1:N){
      asset_name1[i]=paste0(asset_name1[i],i)
    }
  }
  margins_r<<-rep("norm",N)
  paramMargins_r <- list()
  # Now the new experiments
  for(i in 1:N){
    paramMargins_r[[length(paramMargins_r)+1]] <- list(mean =asset_ret[i], sd =sqrt(asset_var[i]))
  }
  paramMargins_r<<-paramMargins_r
  asset_corr<<-asset_corr
  N_sample<-sample_number
  k<<-subjective_k
  rand2<<-generate_N_rand(N,asset_corr,margins_r,paramMargins_r,N_sample) 
  ####below is a way to reduce computation time,when the total scenario number is above 2000,we only pick out the scenarios happen with more than 0.01% probability
# if(ncol(pro_dict)>2000){
  pro_dict1<<-pro_dict[,which(pro_dict[1,]>smallprob)]
  rand_sub<<-generate_N_rand(N,asset_corr,margins_r,paramMargins_r,sub_sample) 
#   }else{
#     pro_dict1<<-pro_dict    
#   }
  # pro_dict1<<-pro_dict 
  #####################################global variable initialization
  loss1<<-extreme_stress_loss
  trans_cost<<-trans_cost
  finan_cost<<-finan_cost
  haircut<<-haircut
  real_finance_weight<<-real_finance_weight
  w_now<<-w_now
  w1<<-w_now
  eval_g0 <<- inequality_constraints
  eval_h0 <<- equality_constraints
  lower_bound<<-lower_bound
  upper_bound<<-upper_bound
  w_now1<<-w_now
  beta11<<-lambda
  trans_cost1<<-trans_cost
  finan_cost1<<-finan_cost
  haircut2<<-haircut
  real_finance_weight1<<-real_finance_weight
  principal11<<-principal1
  rand21<<-rand2
  rand_sub1<<-rand_sub
  loss11<<-loss1
  pro_dictt<<-pro_dict1
  k1<<-k
  mu1<<-mu
  # maxeval_global<<-maxeval_global
  signal<-rep(0,2)
  # rand_sub1<<-rand_sub
  rand2_moment<<-generate_moment_matching(pro_dict1,mu1,k,asset_cov,N_sample)
  rand21_moment<<-rand2_moment
  ###############################calcualte equality and inequality jacobian matrix
  if(is.null(eval_g0)){
    eval_g0_jac<<-NULL
    signal[1]<-0
  }else{
    eval_g0_jac<<-function(w1,w_now,beta1,rand2,rand_sub,rand2_moment,loss1,pro_dict,k,trans_cost,finan_cost,haircut,real_finance_weight,principal1,mu) nl.jacobian(w1, eval_g0, heps = 1e-10)
    signal[1]<-1
    }
  if(is.null(eval_h0)){
    eval_h0_jac<<-NULL
    signal[2]<-0
  }else{
    eval_h0_jac<<-function(w1,w_now,beta1,rand2,rand_sub,rand2_moment,loss1,pro_dict,k,trans_cost,finan_cost,haircut,real_finance_weight,principal1,mu) nl.jacobian(w1, eval_h0, heps = 1e-10)
    signal[1]<-1
  }
  ##########################################initial point is our present position
  start<-w_now
  ####initialize default parameters
  ###############################################
  #######global optimizer AUGLAG+MLSL
  local_opts_global<<-list("algorithm" = "NLOPT_GD_MLSL",
                           "xtol_abs" = xtol_abs_global,
                           "ftol_abs"=ftol_abs_global)
    opts_global <<- list( "algorithm" = "NLOPT_LD_AUGLAG",
                         "xtol_abs" = xtol_abs_global,
                        "ftol_abs"=ftol_abs_global,
                        "maxeval" = maxeval_global,
                        "local_opts" = local_opts_global,
                        "print_level" = print_level)
  ###############################################
  #######local optimizer AUGLAG+SLSQP
   local_opts <<- list( "algorithm" = "NLOPT_LD_SLSQP",
                        "xtol_abs" =xtol_abs,
                        "ftol_abs"=ftol_abs)
    opts <<- list( "algorithm" = "NLOPT_LD_AUGLAG",
                   "xtol_abs" = xtol_abs,
                   "ftol_abs"=ftol_abs,
                   "maxeval" = maxeval_local,
                   "local_opts" = local_opts,
                   "print_level" = print_level)
    opts_moment <<- list( "algorithm" = "NLOPT_LD_SLSQP",
                         "xtol_abs" = 1e-10,
                         "ftol_abs"=1e-10,
                         "maxeval" = maxeval_global,
                         "print_level" = print_level)
  w_start<-initial_point(method,moment,start,lower_bound,upper_bound,w_now,lambda,trans_cost,finan_cost,haircut,real_finance_weight,principal1,rand2,rand_sub,rand2_moment,loss1,pro_dict1,k,mu,signal)
  if(method=='power'){
  w1<-nloptr(x0=w_start,eval_f = power_find_w,eval_grad_f = power_find_w_grad,eval_g_ineq = eval_g0,eval_jac_g_ineq = eval_g0_jac,
                 eval_g_eq = eval_h0,eval_jac_g_eq = eval_h0_jac,lb=lower_bound,ub=upper_bound,opts=opts,w_now=w_now,beta1=lambda,
                 trans_cost=trans_cost,finan_cost=finan_cost,haircut=haircut,real_finance_weight=real_finance_weight,principal1=principal1,rand2=rand2,rand_sub=rand_sub,rand2_moment=rand2_moment,loss1=loss1,pro_dict=pro_dict1,k=k,mu=mu)
    }else if(method=='log'){
        w1<-nloptr(x0=w_start,eval_f = log_find_w,eval_grad_f = log_find_w_grad,eval_g_ineq = eval_g0,eval_jac_g_ineq = eval_g0_jac,
                   eval_g_eq = eval_h0,eval_jac_g_eq = eval_h0_jac,lb=lower_bound,ub=upper_bound,opts=opts,w_now=w_now,beta1=lambda,
                   trans_cost=trans_cost,finan_cost=finan_cost,haircut=haircut,real_finance_weight=real_finance_weight,principal1=principal1,rand2=rand2,rand_sub=rand_sub,rand2_moment=rand2_moment,loss1=loss1,pro_dict=pro_dict1,k=k,mu=mu)
      }else if(method=='expo'){
          w1<-nloptr(x0=w_start,eval_f = expo_find_w,eval_grad_f = expo_find_w_grad,eval_g_ineq = eval_g0,eval_jac_g_ineq = eval_g0_jac,
                     eval_g_eq = eval_h0,eval_jac_g_eq = eval_h0_jac,lb=lower_bound,ub=upper_bound,opts=opts,w_now=w_now,beta1=lambda,
                     trans_cost=trans_cost,finan_cost=finan_cost,haircut=haircut,real_finance_weight=real_finance_weight,principal1=principal1,rand2=rand2,rand_sub=rand_sub,rand2_moment=rand2_moment,loss1=loss1,pro_dict=pro_dict1,k=k,mu=mu)
      }
  w111<<-w1
  w1<-w1$solution
  #####################calcualte total amount of borrow weights
  w1[N+1]=sum(cost2(w_now,w1,trans_cost,finan_cost,haircut,real_finance_weight,principal1)[,3])
#   if(sum(real_finance_weight[which(w1>0)]*(1-haircut[which(w1>0)])*w1[which(w1>0)])+sum(haircut[which(w1>0)]*w1[which(w1>0)])>1){   ###borrow weights larger than 1
#     w1[N+1]=sum(real_finance_weight[which(w1>0)]*(1-haircut[which(w1>0)])*w1[which(w1>0)])+sum(haircut[which(w1>0)]*w1[which(w1>0)])-1
#   }else{                 ###borrow weights smaller than 1
#     w1[N+1]=0
#   }
  method<-rep(method,N+1)  
  weights2<<-data.frame(
    method=method,
    asset_name=c(asset_name1,"borrow"),
    weights=w1
  )
  return(weights2)
}
##########################a function to pre process data for scenario_cost_BN
call_scenario_cost_BN<-function(uti,assets_num1,lambda1,asset_ret1,asset_vol1,
                                asset_corr1,sample_number1,extreme_stress_loss,pro_dict,
                                principal1,trans_cost,finan_cost,haircut,real_finance_weight,w_now,lower_bound,upper_bound,subjective_k,
                                inequality_constraints,equality_constraints,maxeval_global,maxeval_local,asset_name1,moment,smallprob){
  uti<-uti
  assets_num1<-as.double(assets_num1)
  lambda1<-lambda1
  asset_ret1<-as.double(strsplit(asset_ret1,",")[[1]])
  asset_vol1<-as.double(strsplit(asset_vol1,",")[[1]])
  asset_var1<-(asset_vol1)^2
  asset_corr1<-as.double(strsplit(asset_corr1,",")[[1]])
  sample_number1<-as.double(sample_number1)
  extreme_stress_loss<-as.double(strsplit(extreme_stress_loss,",")[[1]])
  principal1<<-as.double(principal1)
  w_now<-as.double(strsplit(w_now,",")[[1]])
  lower_bound<-as.double(strsplit(lower_bound,",")[[1]])
  upper_bound<-as.double(strsplit(upper_bound,",")[[1]])
  k<-as.double(subjective_k)
  if(length(inequality_constraints)==0){
    inequality_constraints=NULL
  }else{
    inequality_constraints<-eval(parse(text=inequality_constraints))}
  if(length(equality_constraints)==0){
    equality_constraints=NULL
  }else{
    equality_constraints<-eval(parse(text=equality_constraints))}
  asset_name1<-strsplit(asset_name1,",")[[1]]
  weights=scenario_cost_BN(uti,assets_num1,lambda1,asset_ret1,asset_var1,
                           asset_corr1,sample_number1,extreme_stress_loss,pro_dict,
                           principal1,trans_cost,finan_cost,haircut,real_finance_weight,w_now,lower_bound,upper_bound,k,
                           inequality_constraints,equality_constraints,maxeval_global,maxeval_local,asset_name1,moment,smallprob)
  return(weights)
}
#####vectorize function for plot
utility<-function(x,x_extreme1=x_extreme,x_downturning1=x_downturning,k_21=k_2,k_11=k_1,x_11=x_1){
  return((ifelse((x<x_downturning1),1,0)*(l+AA*log_utility(1+k_21*x))+ifelse(x>=x_downturning1,1,0)*(1/(1+exp(-k_11*(x-x_11))))))
}
vutility<-Vectorize(utility,vectorize.args="x")
#####################below are functions for combination utility function
function_logistic<-function(x){
  return(1/(1+exp(-x*(x_mid-x_1)))-prob2/(1+exp(-x*(x_downturning-x_1)))-(1-prob2)/(1+exp(-x*(x_upturning-x_1))))
}

utility_solve<-function(x_extreme,x_downturning,x_upturning,x_mid,prob2){
  x_extreme<<-x_extreme
  x_downturning<<-x_downturning
  x_upturning<<-x_upturning
  x_mid<<-x_mid
  prob2<<-prob2
  x_1<<-x_upturning
  k_1<<-uniroot.all(function_logistic,c(0,100))[-1]
  k_2<<--1/x_extreme
  AA<<-k_1*exp(-k_1*(x_downturning-x_1))/((1+exp(-k_1*(x_downturning-x_1)))^2)*(1+k_2*x_downturning)/k_2
  l<<-1/(1+exp(-k_1*(x_downturning-x_1)))-AA*log(1+k_2*x_downturning)
  result<-c(x_downturning,x_1,k_1,k_2,l,AA)
  return(result)
}
combo_utility<-function(x,x_downturning=-0.05,AA= 0.161206,l= 0.0452684,k_2=3.33,k_1=41.28,x_1=0.05){
  return(ifelse(x<x_downturning,1,0)*(l+AA*log_utility(1+k_2*x))+ifelse(x>=x_downturning,1,0)*(1/(1+exp(-k_1*(x-x_1)))))
}
get_combo_ut<-function(ret,x_downturning,AA,l,k_2,k_1,x_1,w_now,w1,tcost){
  return(combo_utility(ret+tcost,x_downturning,AA,l,k_2,k_1,x_1))
}
get_extreme_combo_uti<<-function(pro_dict,loss,x_downturning,AA,l,k_2,k_1,x_1,w_now,w,tcost,mu,rand_sub){
  rett<-0
  u<-rep(1,nrow(pro_dict)-1)
  pro_dict<-as.matrix(pro_dict)
  rett<-sum(1/nrow(rand_sub)*pro_dict[1,2:ncol(pro_dict)]*combo_utility(rep(t(pro_dict)[2:ncol(pro_dict),2:nrow(pro_dict)]%*%(loss*w),nrow(rand_sub))+as.vector((u-t(pro_dict)[2:ncol(pro_dict),2:nrow(pro_dict)])%*%(apply(rand_sub,1,function(x) x*w)))+tcost,x_downturning,AA,l,k_2,k_1,x_1))
  return(rett)
}
combo_find_w<-function(w1,w_now=w_now1,x_downturning=x_downturning1,AA=AA1,ll=ll1,k_2=k_21,k_1=k_11,x_1=x_11,trans_cost=trans_cost1,finan_cost=finan_cost1,haircut=haircut2,real_finance_weight=real_finance_weight1,principal1=principal11,rand2=rand21,rand_sub=rand_sub1,rand2_moment=rand21_moment,loss1=loss11,pro_dict=pro_dictt,k=k1,mu=mu1){
  tcost<-cost(w_now,w1,trans_cost,finan_cost,haircut,real_finance_weight,principal1)
  exp_pow_ut<-(1-k)*mean(na.omit(get_combo_ut(rand2 %*% w1,x_downturning,AA,ll,k_2,k_1,x_1,w_now,w1,tcost)))+(k/(1-pro_dict$p0[1]))*get_extreme_combo_uti(pro_dict,loss1,x_downturning,AA,ll,k_2,k_1,x_1,w_now,w1,tcost,mu,rand_sub)
  return(-exp_pow_ut)
}
combo_find_w_grad<-function(w1,w_now=w_now1,x_downturning=x_downturning1,AA=AA1,ll=ll1,k_2=k_21,k_1=k_11,x_1=x_11,trans_cost=trans_cost1,finan_cost=finan_cost1,haircut=haircut2,real_finance_weight=real_finance_weight1,principal1=principal11,rand2=rand21,rand_sub=rand_sub1,rand2_moment=rand21_moment,loss1=loss11,pro_dict=pro_dictt,k=k1,mu=mu1) nl.grad(w1,combo_find_w,heps=1e-10)
######################################################
###########a function to do optimization based on risk seeking combo utility function
scenario_cost_BN2<-function(assets_num,x_1,k_1,k_2,AA,l,x_downturning,asset_ret,asset_var,asset_corr,sample_number,extreme_stress_loss,pro_dict,principal1,
                           trans_cost,finan_cost,haircut,real_finance_weight,w_now,lower_bound,upper_bound,subjective_k,
                           inequality_constraints,equality_constraints,maxeval_global,maxeval_local,asset_name1=NULL,moment,smallprob){
  method="combo"
  N=assets_num
  mu<<-asset_ret
  asset_cor<-matrix(1,N,N)
  kk=1
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      asset_cor[i,j]=asset_corr[kk]
      kk=kk+1
      asset_cor[j,i]=asset_cor[i,j]
    }
  }
  asset_cov<-r2cov(sd =sqrt(asset_var),R = asset_cor)
  if(length(asset_name1)==0){
    asset_name1=rep("asset_",N)
    for(i in 1:N){
      asset_name1[i]=paste0(asset_name1[i],i)
    }
  }
  margins_r<-rep("norm",N)
  paramMargins_r <- list()
  # Now the new experiments
  for(i in 1:N){
    paramMargins_r[[length(paramMargins_r)+1]] <- list(mean =asset_ret[i], sd =sqrt(asset_var[i]))
  }
  N_sample<-sample_number
  k<<-subjective_k
  rand2<<-generate_N_rand(N,asset_corr,margins_r,paramMargins_r,N_sample) 
  rand_sub<<-generate_N_rand(N,asset_corr,margins_r,paramMargins_r,sub_sample) 
  # if(ncol(pro_dict)>2000){
    pro_dict1<<-pro_dict[,which(pro_dict[1,]>smallprob)]
#   }else{
#     pro_dict1<<-pro_dict    
#   }
  # pro_dict1<<-pro_dict    
  loss1<<-extreme_stress_loss
  trans_cost<<-trans_cost
  finan_cost<<-finan_cost
  haircut<<-haircut
  real_finance_weight<<-real_finance_weight
  w_now<<-w_now
  w1<<-w_now
  eval_g0 <<- inequality_constraints
  eval_h0 <<- equality_constraints
  lower_bound<<-lower_bound
  upper_bound<<-upper_bound
  signal<-rep(0,2)
  if(is.null(eval_g0)){
    eval_g0_combo<<-NULL
    eval_g0_combo_jac<<-NULL
  }else{
    eval_g0_combo<<-function(w1,w_now=w_now1,x_downturning=x_downturning1,AA=AA1,ll=ll1,k_2=k_21,k_1=k_11,x_1=x_11,trans_cost=trans_cost1,finan_cost=finan_cost1,haircut=haircut2,real_finance_weight=real_finance_weight1,principal1=principal11,rand2=rand21,rand_sub=rand_sub1,rand2_moment=rand21_moment,loss1=loss11,pro_dict=pro_dictt,k=k1,mu=mu1)   eval_g0(w1,w_now,beta1,rand2,rand_sub,loss1,pro_dict,k,trans_cost,finan_cost,haircut,real_finance_weight,principal1,mu)
    eval_g0_combo_jac<<-function(w1,w_now,x_downturning,AA,ll,k_2,k_1,x_1,trans_cost,finan_cost,haircut,real_finance_weight,principal1,rand2,rand_sub,rand2_moment,loss1,pro_dict,k,mu) nl.jacobian(w1, eval_g0_combo, heps = 1e-10)
    signal[1]<-1
    }
  if(is.null(eval_h0)){
    eval_h0_combo<<-NULL
    eval_h0_combo_jac<<-NULL
  }else{
    eval_h0_combo<<-function(w1,w_now=w_now1,x_downturning=x_downturning1,AA=AA1,ll=ll1,k_2=k_21,k_1=k_11,x_1=x_11,trans_cost=trans_cost1,finan_cost=finan_cost1,haircut=haircut2,real_finance_weight=real_finance_weight1,principal1=principal11,rand2=rand21,rand_sub=rand_sub1,rand2_moment=rand21_moment,loss1=loss11,pro_dict=pro_dictt,k=k1,mu=mu1)   eval_h0(w1,w_now,beta1,rand2,rand_sub,loss1,pro_dict,k,trans_cost,finan_cost,haircut,real_finance_weight,principal1,mu)
    eval_h0_combo_jac<<-function(w1,w_now,x_downturning,AA,ll,k_2,k_1,x_1,trans_cost,finan_cost,haircut,real_finance_weight,principal1,rand2,rand_sub,rand2_moment,loss1,pro_dict,k,mu) nl.jacobian(w1,eval_h0_combo, heps = 1e-10)
    signal[2]<-1
    }
  w_now1<<-w_now
  x_downturning1<<-x_downturning
  AA1<<-AA
  ll1<<-l
  k_21<<-k_2
  k_11<<-k_1
  x_11<<-x_1
  trans_cost1<<-trans_cost
  finan_cost1<<-finan_cost
  haircut2<<-haircut
  real_finance_weight1<<-real_finance_weight
  principal11<<-principal1
  rand21<<-rand2
  rand_sub1<<-rand_sub
  loss11<<-loss1
  pro_dictt<<-pro_dict
  k1<<-k
  mu1<<-mu
  rand2_moment<<-generate_moment_matching(pro_dict1,mu1,k,asset_cov,N_sample)
  rand21_moment<<-rand2_moment
  start<-w_now
  ###########################opts initialize select parameters
  if(moment==TRUE){
    globalsearch<<--1
    opts_moment <<- list( "algorithm" = "NLOPT_LD_SLSQP",
                   "xtol_abs" = 1e-10,
                   "ftol_abs"=1e-10,
                   "maxeval" = maxeval_global,
                   "print_level" = print_level)
    res3<<-nloptr(x0=start,eval_f = combo_find_reduce,eval_grad_f = combo_find_reduce_grad,eval_g_ineq = eval_g0_combo,eval_jac_g_ineq = eval_g0_combo_jac,
               eval_g_eq = eval_h0_combo,eval_jac_g_eq=eval_h0_combo_jac,lb=lower_bound,ub=upper_bound,opts=opts_moment,w_now=w_now,x_downturning=x_downturning,AA=AA,ll=l,k_2=k_2,k_1=k_1,x_1=x_1,
               trans_cost=trans_cost,finan_cost=finan_cost,haircut=haircut,real_finance_weight=real_finance_weight,principal1=principal1,rand2=rand2,rand_sub=rand_sub,rand2_moment=rand2_moment,loss1=loss1,pro_dict=pro_dict1,k=k,mu=mu)
    w_start<-res3$solution
  }else{
    local_opts <<- list( "algorithm" = "NLOPT_LD_SLSQP",
                         "xtol_abs" =xtol_abs,
                         "ftol_abs"=ftol_abs)
    local_opts_global<<-list("algorithm" = "NLOPT_GD_MLSL",
                             "xtol_abs" = xtol_abs_global,
                             "ftol_abs"=ftol_abs_global)
    opts_global <<- list( "algorithm" = "NLOPT_LD_AUGLAG",
                          "xtol_abs" = xtol_abs_global,
                          "ftol_abs"=ftol_abs_global,
                          "maxeval" = maxeval_global,
                          "local_opts" = local_opts_global,   
                          "print_level" = print_level)
    opts <<- list( "algorithm" = "NLOPT_LD_AUGLAG",
                   "xtol_abs" = xtol_abs,
                   "ftol_abs"=ftol_abs,
                   "maxeval" = maxeval_local,
                   "local_opts" = local_opts,
                   "print_level" = print_level)
    buffer<<-capture.output(res3<<-nloptr(x0=start,eval_f = combo_find_w,eval_grad_f = combo_find_w_grad,eval_g_ineq = eval_g0_combo,eval_jac_g_ineq = eval_g0_combo_jac,
                                          eval_g_eq = eval_h0_combo,eval_jac_g_eq = eval_h0_combo_jac,lb=lower_bound,ub=upper_bound,opts=opts_global,w_now=w_now,x_downturning=x_downturning,AA=AA,ll=l,k_2=k_2,k_1=k_1,x_1=x_1,
                                          trans_cost=trans_cost,finan_cost=finan_cost,haircut=haircut,real_finance_weight=real_finance_weight,principal1=principal1,rand2=rand2,rand_sub=rand_sub,rand2_moment=rand2_moment,loss1=loss1,pro_dict=pro_dict1,k=k,mu=mu))
    if(sum(ifelse(res3$solution==start,0,1))==0||res3$status==-1){
      globalsearch<<-1
      w_start<-find_optimal(buffer,signal)
    }else{
      globalsearch<<-0
      w_start<-res3$solution
    }
  }   
  w1<-nloptr(x0=w_start,eval_f = combo_find_w,eval_grad_f = combo_find_w_grad,eval_g_ineq = eval_g0_combo,eval_jac_g_ineq = eval_g0_combo_jac,
               eval_g_eq = eval_h0_combo,eval_jac_g_eq=eval_h0_combo_jac,lb=lower_bound,ub=upper_bound,opts=opts,w_now=w_now,x_downturning=x_downturning,AA=AA,ll=l,k_2=k_2,k_1=k_1,x_1=x_1,
               trans_cost=trans_cost,finan_cost=finan_cost,haircut=haircut,real_finance_weight=real_finance_weight,principal1=principal1,rand2=rand2,rand_sub=rand_sub,rand2_moment=rand2_moment,loss1=loss1,pro_dict=pro_dict1,k=k,mu=mu)

  # stopCluster(cl)
  w1<-w1$solution
  w1[N+1]=sum(cost2(w_now,w1,trans_cost,finan_cost,haircut,real_finance_weight,principal1)[,3])
# #   w=ifelse(length(which(w1<0)),abs(sum(((1-haircut)*real_finance_weight*w1)[which(w1<0)]))+abs(sum((haircut*w1)[which(w1<0)])),0)
# #   w1[length(w1)]=ifelse(w1[length(w1)]>w,w1[length(w1)],w)
#   if(sum(real_finance_weight[which(w1>0)]*(1-haircut[which(w1>0)])*w1[which(w1>0)])+sum(haircut[which(w1>0)]*w1[which(w1>0)])>1){
#     w1[N+1]=sum(real_finance_weight[which(w1>0)]*(1-haircut[which(w1>0)])*w1[which(w1>0)])+sum(haircut[which(w1>0)]*w1[which(w1>0)])-1
#   }else{
#     w1[N+1]=0
#     # w1[N]=w1[N]+1-sum(real_finance_weight[which(w1>0)]*w1[which(w1>0)])
#   }
  method<-rep(method,N+1)  
  weights2<<-data.frame(
    method=method,
    asset_name=c(asset_name1,"borrow"),
    weights=w1
  )
  return(weights2)
}
call_scenario_cost_BN2<-function(assets_num1,x_extreme,x_downturning,x_mid2,x_upturning,prob2,asset_ret1,asset_vol1,
                                 asset_corr1,sample_number1,extreme_stress_loss,pro_dict,
                                 principal1,trans_cost,finan_cost,haircut,real_finance_weight,w_now,lower_bound,upper_bound,subjective_k,
                                 inequality_constraints,equality_constraints,maxeval_global,maxeval_local,asset_name1,moment,smallprob){
  assets_num1<-as.double(assets_num1)
  x_extreme<<-x_extreme
  x_downturning<<-x_downturning
  x_upturning<<-x_upturning
  x_mid<<-x_mid2
  prob2<<-prob2
  x_1<<-x_upturning
  k_1<<-uniroot.all(function_logistic,c(0,100))[-1]
  k_2<<--1/x_extreme
  AA<<-k_1*exp(-k_1*(x_downturning-x_1))/((1+exp(-k_1*(x_downturning-x_1)))^2)*(1+k_2*x_downturning)/k_2
  l<<-1/(1+exp(-k_1*(x_downturning-x_1)))-AA*log(1+k_2*x_downturning)
  asset_ret1<-as.double(strsplit(asset_ret1,",")[[1]])
  asset_vol1<-as.double(strsplit(asset_vol1,",")[[1]])
  asset_var1<-(asset_vol1)^2
  asset_corr1<-as.double(strsplit(asset_corr1,",")[[1]])
  sample_number1<-as.double(sample_number1)
  extreme_stress_loss<-as.double(strsplit(extreme_stress_loss,",")[[1]])
  # cond_table<-build_cond(affect_relation,prob_table)
  principal1<<-as.double(principal1)
  w_now<-as.double(strsplit(w_now,",")[[1]])
  lower_bound<-as.double(strsplit(lower_bound,",")[[1]])
  upper_bound<-as.double(strsplit(upper_bound,",")[[1]])
  k<-as.double(subjective_k)
  if(length(inequality_constraints)==0){
    inequality_constraints=NULL
  }else{
    inequality_constraints<-eval(parse(text=inequality_constraints))}
  if(length(equality_constraints)==0){
    equality_constraints=NULL
  }else{
    equality_constraints<-eval(parse(text=equality_constraints))}
  asset_name1<-strsplit(asset_name1,",")[[1]]
  weights=scenario_cost_BN2(assets_num1,x_1,k_1,k_2,AA,l,x_downturning,asset_ret1,asset_var1,
                           asset_corr1,sample_number1,extreme_stress_loss,pro_dict,
                           principal1,trans_cost,finan_cost,haircut,real_finance_weight,w_now,lower_bound,upper_bound,k,
                           inequality_constraints,equality_constraints,maxeval_global,maxeval_local,asset_name1,moment,smallprob)
  return(weights)
}
scenario_cost_BN_matrix<-function(method,assets_num,lambdas,asset_ret,asset_var,asset_corr,sample_number,extreme_stress_loss,pro_dict,principal1,
                           trans_cost,finan_cost,haircut,real_finance_weight,w_now,lower_bound,upper_bound,subjective_ks,
                           inequality_constraints,equality_constraints,maxeval_global,maxeval_local,asset_name1=NULL,moment,smallprob){
  ##
  #
  ##trans_cost matrix like
  N=assets_num
  # A=lambda
  mu<<-asset_ret
  asset_cor<-matrix(1,N,N)
  kk=1
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      asset_cor[i,j]=asset_corr[kk]
      kk=kk+1
      asset_cor[j,i]=asset_cor[i,j]
    }
  }
#   asset_cor[upper.tri(asset_cor)]=asset_corr
#   asset_cor[lower.tri(asset_cor)]=t(asset_cor)[lower.tri(asset_cor)]
  asset_cov<-r2cov(sd =sqrt(asset_var),R = asset_cor)
  if(length(asset_name1)==0){
    asset_name1=rep("asset_",N)
    for(i in 1:N){
      asset_name1[i]=paste0(asset_name1[i],i)
    }
  }
  margins_r<-rep("norm",N)
  paramMargins_r <- list()
  # Now the new experiments
  for(i in 1:N){
    paramMargins_r[[length(paramMargins_r)+1]] <- list(mean =asset_ret[i], sd =sqrt(asset_var[i]))
  }
  N_sample<-sample_number
  k<<-subjective_ks
  rand2<<-generate_N_rand(N,asset_corr,margins_r,paramMargins_r,N_sample) 
  rand_sub<<-generate_N_rand(N,asset_corr,margins_r,paramMargins_r,sub_sample) 
  ###########################reduce total scenario to scenarios happen with probability larger than 0.01%
 # if(ncol(pro_dict)>2000){
    pro_dict1<<-pro_dict[,which(pro_dict[1,]>smallprob)]
#   }else{
#     pro_dict1<<-pro_dict    
#   }
  # pro_dict1<<-pro_dict 
  loss1<<-extreme_stress_loss
  trans_cost<<-trans_cost
  finan_cost<<-finan_cost
  haircut<<-haircut
  real_finance_weight<<-real_finance_weight
  w_now<<-w_now
  w1<<-w_now
  eval_g0 <<- inequality_constraints
  eval_h0 <<- equality_constraints
  lower_bound<<-lower_bound
  upper_bound<<-upper_bound
  signal<-rep(0,2)
  if(is.null(eval_g0)){
    eval_g0_jac<<-NULL
  }else{
    eval_g0_jac<<-function(w1,w_now,beta1,rand2,rand_sub,rand2_moment,loss1,pro_dict,k,trans_cost,finan_cost,haircut,real_finance_weight,principal1,mu) nl.jacobian(w1, eval_g0, heps = 1e-10)
   signal[1]<-1
    }
  if(is.null(eval_h0)){
    eval_h0_jac<<-NULL
  }else{
    eval_h0_jac<<-function(w1,w_now,beta1,rand2,rand_sub,rand2_moment,loss1,pro_dict,k,trans_cost,finan_cost,haircut,real_finance_weight,principal1,mu) nl.jacobian(w1, eval_h0, heps = 1e-10)
   signal[2]<-1
    }
  w_now1<<-w_now
  trans_cost1<<-trans_cost
  finan_cost1<<-finan_cost
  haircut2<<-haircut
  real_finance_weight1<<-real_finance_weight
  principal11<<-principal1
  rand21<<-rand2
  rand_sub1<<-rand_sub
  loss11<<-loss1
  pro_dictt<<-pro_dict1
  mu1<<-mu
  local_opts <<- list( "algorithm" = "NLOPT_LD_SLSQP",
                       "xtol_abs" =xtol_abs,
                       "ftol_abs"=ftol_abs)
  local_opts_global<<-list("algorithm" = "NLOPT_GD_MLSL",
                           "xtol_abs" = xtol_abs_global,
                           "ftol_abs"=ftol_abs_global)
  opts_global <<- list( "algorithm" = "NLOPT_LD_AUGLAG",
                        "xtol_abs" = xtol_abs_global,
                        "ftol_abs"=ftol_abs_global,
                        "maxeval" = maxeval_global,
                        "local_opts" = local_opts_global,    
                        "print_level" = print_level)
  opts_moment <<- list( "algorithm" = "NLOPT_LD_SLSQP",
                 "xtol_abs" = 1e-10,
                 "ftol_abs"=1e-10,
                 "maxeval" = maxeval_global,
                 "print_level" = print_level)
  opts <<- list( "algorithm" = "NLOPT_LD_AUGLAG",
                 "xtol_abs" = xtol_abs,
                 "ftol_abs"=ftol_abs,
                 "maxeval" = maxeval_local,
                 "local_opts" = local_opts,
                 "print_level" = print_level)
  start<-w_now
    for(kk in 1:length(method)){
    if(method[kk]=='power'){
      for(i in 1:length(lambdas)){
        for(j in 1:length(subjective_ks)){
          rand2_moment<<-generate_moment_matching(pro_dict1,mu1,subjective_ks[j],asset_cov,N_sample)
          rand21_moment<<-rand2_moment
          beta11<<-lambdas[i]
          k1<<-subjective_ks[j]
          if(moment==TRUE){
            globalsearch<<-(-1)
           res3<<-nloptr(x0=start,eval_f = power_find_reduce,eval_grad_f = power_find_reduce_grad,eval_g_ineq = eval_g0,eval_jac_g_ineq = eval_g0_jac,
                         eval_g_eq = eval_h0,eval_jac_g_eq = eval_h0_jac,lb=lower_bound,ub=upper_bound,opts=opts_moment,w_now=w_now,beta1=lambdas[i],
                         trans_cost=trans_cost,finan_cost=finan_cost,haircut=haircut,real_finance_weight=real_finance_weight,principal1=principal1,rand2=rand2,rand_sub=rand_sub,rand2_moment=rand2_moment,loss1=loss1,pro_dict=pro_dict1,k=subjective_ks[j],mu=mu)
           w_temp<-nloptr(x0=res3$solution,eval_f = power_find_w,eval_grad_f = power_find_w_grad,eval_g_ineq = eval_g0,eval_jac_g_ineq = eval_g0_jac,
                          eval_g_eq = eval_h0,eval_jac_g_eq = eval_h0_jac,lb=lower_bound,ub=upper_bound,opts=opts,w_now=w_now,beta1=lambdas[i],
                          trans_cost=trans_cost,finan_cost=finan_cost,haircut=haircut,real_finance_weight=real_finance_weight,principal1=principal1,rand2=rand2,rand_sub=rand_sub,rand2_moment=rand2_moment,loss1=loss1,pro_dict=pro_dict1,k=subjective_ks[j],mu=mu)
           }else{
             buffer<<-capture.output(res3<<-nloptr(x0=start,eval_f = power_find_w,eval_grad_f = power_find_w_grad,eval_g_ineq = eval_g0,eval_jac_g_ineq = eval_g0_jac,
                                                   eval_g_eq = eval_h0,eval_jac_g_eq = eval_h0_jac,lb=lower_bound,ub=upper_bound,opts=opts_global,w_now=w_now,beta1=lambdas[i],
                                                   trans_cost=trans_cost,finan_cost=finan_cost,haircut=haircut,real_finance_weight=real_finance_weight,principal1=principal1,rand2=rand2,rand_sub=rand_sub,rand2_moment=rand2_moment,loss1=loss1,pro_dict=pro_dict1,k=subjective_ks[j],mu=mu))
             if(sum(ifelse(res3$solution==start,0,1))==0||res3$status==-1){
               globalsearch<<-1
               w_temp<-nloptr(x0=find_optimal(buffer,signal),eval_f = power_find_w,eval_grad_f = power_find_w_grad,eval_g_ineq = eval_g0,eval_jac_g_ineq = eval_g0_jac,
                              eval_g_eq = eval_h0,eval_jac_g_eq = eval_h0_jac,lb=lower_bound,ub=upper_bound,opts=opts,w_now=w_now,beta1=lambdas[i],
                              trans_cost=trans_cost,finan_cost=finan_cost,haircut=haircut,real_finance_weight=real_finance_weight,principal1=principal1,rand2=rand2,rand_sub=rand_sub,rand2_moment=rand2_moment,loss1=loss1,pro_dict=pro_dict1,k=subjective_ks[j],mu=mu)
             }
             w_temp<-nloptr(x0=res3$solution,eval_f = power_find_w,eval_grad_f = power_find_w_grad,eval_g_ineq = eval_g0,eval_jac_g_ineq = eval_g0_jac,
                            eval_g_eq = eval_h0,eval_jac_g_eq = eval_h0_jac,lb=lower_bound,ub=upper_bound,opts=opts,w_now=w_now,beta1=lambdas[i],
                            trans_cost=trans_cost,finan_cost=finan_cost,haircut=haircut,real_finance_weight=real_finance_weight,principal1=principal1,rand2=rand2,rand_sub=rand_sub,rand2_moment=rand2_moment,loss1=loss1,pro_dict=pro_dict1,k=subjective_ks[j],mu=mu)
             
            }
           w_temp<-w_temp$solution
           w_temp[N+1]=sum(cost2(w_now,w_temp,trans_cost,finan_cost,haircut,real_finance_weight,principal1)[,3])
          if(kk==1&&i==1&&j==1){
            w_temp<-c(lambdas[i],subjective_ks[j],w_temp)
            weights_matrix<-data.frame(t(w_temp))
            colnames(weights_matrix)<-c("lambda","k",asset_name1,"borrow")
          }else{
            w_temp<-c(lambdas[i],subjective_ks[j],w_temp)
            weights_matrix<-rbind(weights_matrix,w_temp)
          }
        }
      }
    }else if(method[kk]=='log'){
      for(i in 1:length(lambdas)){
        for(j in 1:length(subjective_ks)){
          beta11<<-lambdas[i]
          k1<<-subjective_ks[j]
          rand2_moment<<-generate_moment_matching(pro_dict1,mu1,subjective_ks[j],asset_cov,N_sample)
          rand21_moment<<-rand2_moment
          start<-w_now
          if(moment==TRUE){
            res3<<-nloptr(x0=start,eval_f = log_find_reduce,eval_grad_f = log_find_reduce_grad,eval_g_ineq = eval_g0,eval_jac_g_ineq = eval_g0_jac,
                           eval_g_eq = eval_h0,eval_jac_g_eq = eval_h0_jac,lb=lower_bound,ub=upper_bound,opts=opts_moment,w_now=w_now,beta1=lambdas[i],
                           trans_cost=trans_cost,finan_cost=finan_cost,haircut=haircut,real_finance_weight=real_finance_weight,principal1=principal1,rand2=rand2,rand_sub=rand_sub,rand2_moment=rand2_moment,loss1=loss1,pro_dict=pro_dict1,k=subjective_ks[j],mu=mu)
            w_temp<-nloptr(x0=res3$solution,eval_f = log_find_w,eval_grad_f = log_find_w_grad,eval_g_ineq = eval_g0,eval_jac_g_ineq = eval_g0_jac,
                           eval_g_eq = eval_h0,eval_jac_g_eq = eval_h0_jac,lb=lower_bound,ub=upper_bound,opts=opts,w_now=w_now,beta1=lambdas[i],
                           trans_cost=trans_cost,finan_cost=finan_cost,haircut=haircut,real_finance_weight=real_finance_weight,principal1=principal1,rand2=rand2,rand_sub=rand_sub,rand2_moment=rand2_moment,loss1=loss1,pro_dict=pro_dict1,k=subjective_ks[j],mu=mu)
          }else{
           buffer<<-capture.output(res3<<-nloptr(x0=start,eval_f = log_find_w,eval_grad_f = log_find_w_grad,eval_g_ineq = eval_g0,eval_jac_g_ineq = eval_g0_jac,
                        eval_g_eq = eval_h0,eval_jac_g_eq = eval_h0_jac,lb=lower_bound,ub=upper_bound,opts=opts_global,w_now=w_now,beta1=lambdas[i],
                        trans_cost=trans_cost,finan_cost=finan_cost,haircut=haircut,real_finance_weight=real_finance_weight,principal1=principal1,rand2=rand2,rand_sub=rand_sub,rand2_moment=rand2_moment,loss1=loss1,pro_dict=pro_dict1,k=subjective_ks[j],mu=mu))
          if(sum(ifelse(res3$solution==start,0,1))==0||res3$status==-1){
            globalsearch<<-1
            w_temp<-nloptr(x0=find_optimal(buffer,signal),eval_f = log_find_w,eval_grad_f = log_find_w_grad,eval_g_ineq = eval_g0,eval_jac_g_ineq = eval_g0_jac,
                           eval_g_eq = eval_h0,eval_jac_g_eq = eval_h0_jac,lb=lower_bound,ub=upper_bound,opts=opts,w_now=w_now,beta1=lambdas[i],
                           trans_cost=trans_cost,finan_cost=finan_cost,haircut=haircut,real_finance_weight=real_finance_weight,principal1=principal1,rand2=rand2,rand_sub=rand_sub,rand2_moment=rand2_moment,loss1=loss1,pro_dict=pro_dict1,k=subjective_ks[j],mu=mu)
          }
          w_temp<-nloptr(x0=res3$solution,eval_f = log_find_w,eval_grad_f = log_find_w_grad,eval_g_ineq = eval_g0,eval_jac_g_ineq = eval_g0_jac,
                         eval_g_eq = eval_h0,eval_jac_g_eq = eval_h0_jac,lb=lower_bound,ub=upper_bound,opts=opts,w_now=w_now,beta1=lambdas[i],
                         trans_cost=trans_cost,finan_cost=finan_cost,haircut=haircut,real_finance_weight=real_finance_weight,principal1=principal1,rand2=rand2,rand_sub=rand_sub,rand2_moment=rand2_moment,loss1=loss1,pro_dict=pro_dict1,k=subjective_ks[j],mu=mu)
          
          }
           w_temp<-w_temp$solution
           w_temp[N+1]=sum(cost2(w_now,w_temp,trans_cost,finan_cost,haircut,real_finance_weight,principal1)[,3])
#           if( sum(real_finance_weight[which(w_temp>0)]*(1-haircut[which(w_temp>0)])*w_temp[which(w_temp>0)])+sum(haircut[which(w_temp>0)]*w_temp[which(w_temp>0)])>1){
#             w_temp[N+1]= sum(real_finance_weight[which(w_temp>0)]*(1-haircut[which(w_temp>0)])*w_temp[which(w_temp>0)])+sum(haircut[which(w_temp>0)]*w_temp[which(w_temp>0)])-1
#           }else{
#             w_temp[N+1]=0
#             # w_temp[N]=w_temp[N]+1-sum(real_finance_weight[which(w_temp>0)]*w_temp[which(w_temp>0)])
#           }
          if(kk==1&&i==1&&j==1){
            w_temp<-c(lambdas[i],subjective_ks[j],w_temp)
            weights_matrix<-data.frame(t(w_temp))
            colnames(weights_matrix)<-c("lambda","k",asset_name1,"borrow")
          }else{
            w_temp<-c(lambdas[i],subjective_ks[j],w_temp)
            weights_matrix<-rbind(weights_matrix,w_temp)
          }
        }
      }
      
    }else if(method[kk]=='expo'){
      for(i in 1:length(lambdas)){
        for(j in 1:length(subjective_ks)){
          beta11<<-lambdas[i]
          k1<<-subjective_ks[j]
          rand2_moment<<-generate_moment_matching(pro_dict1,mu1,subjective_ks[j],asset_cov,N_sample)
          rand21_moment<<-rand2_moment
          start<-w_now
          if(moment==TRUE){
            res3<<-nloptr(x0=start,eval_f = expo_find_reduce,eval_grad_f = expo_find_reduce_grad,eval_g_ineq = eval_g0,eval_jac_g_ineq = eval_g0_jac,
                          eval_g_eq = eval_h0,eval_jac_g_eq = eval_h0_jac,lb=lower_bound,ub=upper_bound,opts=opts_moment,w_now=w_now,beta1=lambdas[i],
                          trans_cost=trans_cost,finan_cost=finan_cost,haircut=haircut,real_finance_weight=real_finance_weight,principal1=principal1,rand2=rand2,rand_sub=rand_sub,rand2_moment=rand2_moment,loss1=loss1,pro_dict=pro_dict1,k=subjective_ks[j],mu=mu)
            w_temp<-nloptr(x0=res3$solution,eval_f = expo_find_w,eval_grad_f = expo_find_w_grad,eval_g_ineq = eval_g0,eval_jac_g_ineq = eval_g0_jac,
                           eval_g_eq = eval_h0,eval_jac_g_eq = eval_h0_jac,lb=lower_bound,ub=upper_bound,opts=opts,w_now=w_now,beta1=lambdas[i],
                           trans_cost=trans_cost,finan_cost=finan_cost,haircut=haircut,real_finance_weight=real_finance_weight,principal1=principal1,rand2=rand2,rand_sub=rand_sub,rand2_moment=rand2_moment,loss1=loss1,pro_dict=pro_dict1,k=subjective_ks[j],mu=mu)
          }else{
           buffer<<-capture.output(res3<<-nloptr(x0=start,eval_f = expo_find_w,eval_grad_f = expo_find_w_grad,eval_g_ineq = eval_g0,eval_jac_g_ineq = eval_g0_jac,
                        eval_g_eq = eval_h0,eval_jac_g_eq = eval_h0_jac,lb=lower_bound,ub=upper_bound,opts=opts_global,w_now=w_now,beta1=lambdas[i],
                        trans_cost=trans_cost,finan_cost=finan_cost,haircut=haircut,real_finance_weight=real_finance_weight,principal1=principal1,rand2=rand2,rand_sub=rand_sub,rand2_moment=rand2_moment,loss1=loss1,pro_dict=pro_dict1,k=subjective_ks[j],mu=mu))
          if(sum(ifelse(res3$solution==start,0,1))==0||res3$status==-1){
            globalsearch<<-1
            w_temp<-nloptr(x0=find_optimal(buffer,signal),eval_f = expo_find_w,eval_grad_f = expo_find_w_grad,eval_g_ineq = eval_g0,eval_jac_g_ineq = eval_g0_jac,
                           eval_g_eq = eval_h0,eval_jac_g_eq = eval_h0_jac,lb=lower_bound,ub=upper_bound,opts=opts,w_now=w_now,beta1=lambdas[i],
                           trans_cost=trans_cost,finan_cost=finan_cost,haircut=haircut,real_finance_weight=real_finance_weight,principal1=principal1,rand2=rand2,rand_sub=rand_sub,rand2_moment=rand2_moment,loss1=loss1,pro_dict=pro_dict1,k=subjective_ks[j],mu=mu)
            
          }
          w_temp<-nloptr(x0=res3$solution,eval_f = expo_find_w,eval_grad_f = expo_find_w_grad,eval_g_ineq = eval_g0,eval_jac_g_ineq = eval_g0_jac,
                         eval_g_eq = eval_h0,eval_jac_g_eq = eval_h0_jac,lb=lower_bound,ub=upper_bound,opts=opts,w_now=w_now,beta1=lambdas[i],
                         trans_cost=trans_cost,finan_cost=finan_cost,haircut=haircut,real_finance_weight=real_finance_weight,principal1=principal1,rand2=rand2,rand_sub=rand_sub,rand2_moment=rand2_moment,loss1=loss1,pro_dict=pro_dict1,k=subjective_ks[j],mu=mu)
 
          }
           w_temp<-w_temp$solution
           w_temp[N+1]=sum(cost2(w_now,w_temp,trans_cost,finan_cost,haircut,real_finance_weight,principal1)[,3])
#            if(sum(real_finance_weight[which(w_temp>0)]*(1-haircut[which(w_temp>0)])*w_temp[which(w_temp>0)])+sum(haircut[which(w_temp>0)]*w_temp[which(w_temp>0)])>1){
#              w_temp[N+1]= sum(real_finance_weight[which(w_temp>0)]*(1-haircut[which(w_temp>0)])*w_temp[which(w_temp>0)])+sum(haircut[which(w_temp>0)]*w_temp[which(w_temp>0)])-1
#            }else{
#              w_temp[N+1]=0
#              # w_temp[N]=w_temp[N]+1-sum(real_finance_weight[which(w_temp>0)]*w_temp[which(w_temp>0)])
#            }
          if(kk==1&&i==1&&j==1){
            w_temp<-c(lambdas[i],subjective_ks[j],w_temp)
            weights_matrix<-data.frame(t(w_temp))
            colnames(weights_matrix)<-c("lambda","k",asset_name1,"borrow")
          }else{
            w_temp<-c(lambdas[i],subjective_ks[j],w_temp)
            weights_matrix<-rbind(weights_matrix,w_temp)
          }
        }
      }  }
    # stopCluster(cl)
  }
  methods<-c()
  for(kk in 1:length(method)){
    methods<-c(methods,rep(method[kk],length(lambdas)*length(subjective_ks)))
  }
  weights_matrix<-cbind(methods,weights_matrix)
  return(weights_matrix)
}
call_scenario_cost_BN_matrix<-function(uti2,assets_num1,lambdas,asset_ret1,asset_vol1,
                             asset_corr1,sample_number2,extreme_stress_loss,pro_dict,
                             principal1,trans_cost,finan_cost,haircut,real_finance_weight,w_now,lower_bound,upper_bound,ks,
                             inequality_constraints,equality_constraints,maxeval_global,maxeval_local,asset_name1,moment,smallprob){
  uti<-uti2
  assets_num1<-as.double(assets_num1)
  asset_ret1<-as.double(strsplit(asset_ret1,",")[[1]])
  asset_vol1<-as.double(strsplit(asset_vol1,",")[[1]])
  asset_var1<-(asset_vol1)^2
  asset_corr1<-as.double(strsplit(asset_corr1,",")[[1]])
  sample_number2<-as.double(sample_number2)
  extreme_stress_loss<-as.double(strsplit(extreme_stress_loss,",")[[1]])
  # cond_table<-build_cond(affect_relation,prob_table)
  principal1<<-as.double(principal1)
  w_now<-as.double(strsplit(w_now,",")[[1]])
  lower_bound<-as.double(strsplit(lower_bound,",")[[1]])
  upper_bound<-as.double(strsplit(upper_bound,",")[[1]])
  lambdas=as.double(strsplit(lambdas,",")[[1]])
  ks<-as.double(strsplit(ks,",")[[1]])
  if(length(inequality_constraints)==0){
    inequality_constraints=NULL
  }else{
    inequality_constraints<-eval(parse(text=inequality_constraints))}
  if(length(equality_constraints)==0){
    equality_constraints=NULL
  }else{
    equality_constraints<-eval(parse(text=equality_constraints))}
  asset_name1<-strsplit(asset_name1,",")[[1]]
  weights=scenario_cost_BN_matrix(uti,assets_num1,lambdas,asset_ret1,asset_var1,
                           asset_corr1,sample_number2,extreme_stress_loss,pro_dict,
                           principal1,trans_cost,finan_cost,haircut,real_finance_weight,w_now,lower_bound,upper_bound,ks,
                           inequality_constraints,equality_constraints,maxeval_global,maxeval_local,asset_name1,moment,smallprob)
  return(weights)
}
call_scenario_cost_BN_matrix2<-function(assets_num1,x_extreme,x_downturning,x_mid2,x_upturning,prob2,asset_ret1,asset_vol1,
                                        asset_corr1,sample_number1,extreme_stress_loss,pro_dict,
                                        principal1,trans_cost,finan_cost,haircut,real_finance_weight,w_now,lower_bound,upper_bound,ks,
                                        inequality_constraints,equality_constraints,maxeval_global,maxeval_local,asset_name1,moment,smallprob){
  assets_num1<-as.double(assets_num1)
  x_extreme<<-x_extreme
  x_downturning<<-x_downturning
  x_upturning<<-x_upturning
  x_mid<<-x_mid2
  prob2<<-prob2
  x_1<<-x_upturning
  k_1<<-uniroot.all(function_logistic,c(0,100))[-1]
  k_2<<--1/x_extreme
  AA<<-k_1*exp(-k_1*(x_downturning-x_1))/((1+exp(-k_1*(x_downturning-x_1)))^2)*(1+k_2*x_downturning)/k_2
  l<<-1/(1+exp(-k_1*(x_downturning-x_1)))-AA*log(1+k_2*x_downturning)
  asset_ret1<-as.double(strsplit(asset_ret1,",")[[1]])
  asset_vol1<-as.double(strsplit(asset_vol1,",")[[1]])
  asset_var1<-(asset_vol1)^2
  asset_corr1<-as.double(strsplit(asset_corr1,",")[[1]])
  sample_number1<-as.double(sample_number1)
  extreme_stress_loss<-as.double(strsplit(extreme_stress_loss,",")[[1]])
  # cond_table<-build_cond(affect_relation,prob_table)
  principal1<<-as.double(principal1)
  w_now<-as.double(strsplit(w_now,",")[[1]])
  lower_bound<-as.double(strsplit(lower_bound,",")[[1]])
  upper_bound<-as.double(strsplit(upper_bound,",")[[1]])
  ks<-as.double(strsplit(ks,",")[[1]])
  if(length(inequality_constraints)==0){
    inequality_constraints=NULL
  }else{
    inequality_constraints<-eval(parse(text=inequality_constraints))}
  if(length(equality_constraints)==0){
    equality_constraints=NULL
  }else{
    equality_constraints<-eval(parse(text=equality_constraints))}
  asset_name1<-strsplit(asset_name1,",")[[1]]
  weights<-matrix(0,nrow=length(ks),ncol=assets_num1+7)
  for(jj in 1:length(ks)){
    k=ks[jj]
    weights3=scenario_cost_BN2(assets_num1,x_1,k_1,k_2,AA,l,x_downturning,asset_ret1,asset_var1,
                            asset_corr1,sample_number1,extreme_stress_loss,pro_dict,
                            principal1,trans_cost,finan_cost,haircut,real_finance_weight,w_now,lower_bound,upper_bound,k,
                            inequality_constraints,equality_constraints,maxeval_global,maxeval_local,asset_name1,moment,smallprob)
    weights[jj,]=c(ks[jj],x_extreme,x_downturning,x_upturning,x_mid,prob2,weights3$weights)
  }
  colnames(weights)<-c("k","extreme","downside","upside","mid","prob",as.character(weights2$asset_name))
  return(weights)
}
getwarningmessage<-function(w1,w_now,w_finan,haircut,real_finance_weight,finan_cost,lower_bound,upper_bound){
  if(length(which(ifelse(w1>=lower_bound,0,1)==1))!=0){
    M<-which(ifelse(w1>=lower_bound,0,1)==1)
    return(-1-0.0001*M[1])
  }
  if(length(which(ifelse(w1<=upper_bound,0,1)==1))!=0){
    M<-which(ifelse(w1<=upper_bound,0,1)==1)
    return(-2-0.0001*M[1])
  }
  if(max(eval_g0(w1))>=tolerance){
    M<-which.max(eval_g0(w1))
    return(-3-0.0001*M[1])
  }  
  if(length(eval_h0)!=0){
    if(max(eval_h0(w1))>tolerance||min(eval_h0(w1))<(-tolerance)) return(-4)
  } 
  list1=which(real_finance_weight!=1)
  haircut1=(1-haircut)*real_finance_weight*w1
  h=ifelse(length(list1)!=0,w1[length(w1)]-sum(abs(((1-haircut2)*real_finance_weight*w1)[list1]))-sum(abs((haircut*w1)[list1]))+tolerance,1)
  # if(w_finan[1]<0) return(-2)
  if(h<0) {
    return(1)
  }
  if(length(list1)!=0){
     M<-which(ifelse(w_finan[-list1]<=(1-haircut[-list1])*w1[-list1]*real_finance_weight[-list1]+tolerance,0,1)==1)+1
  }else{
    M<-which(ifelse(w_finan<=(1-haircut)*w1*real_finance_weight+tolerance,0,1)==1)+1 
  }
  if(length(M)==0) {
    N_fina_con<-ncol(finan_cost)
    w=0
    list1=which(real_finance_weight!=1)
    if(length(list1)!=0){
      w=sum(abs(((1-haircut)*real_finance_weight*w1)[list1]))+sum(abs((haircut*w1)[list1]))
      if(w1[length(w1)]<w){
        flag=1
        w1[length(w1)]=w1[length(w1)]-w
      }else{
        w1[length(w1)]=w1[length(w1)]-w
      }
    }
    con2<-colnames(finan_cost)[2:N_fina_con]
    list=which(real_finance_weight!=1)
    # list=unique(c(list,which(real_finance_weight==0)))
    if(length(list)!=0){
      leverage1<-sum(haircut[-list]*w1[-list]+(1-haircut[-list])*real_finance_weight[-list]*w1[-list])-(1-w)
    }else{
      leverage1<-sum(haircut*w1+(1-haircut)*real_finance_weight*w1)-(1-w)
    }
    if(leverage1<0){
      return(0)
    }else if(sum(w_finan)<=leverage1+tolerance&&sum(w_finan)>=leverage1-tolerance){
      return(0)
      }else{
        return(-1)
      }}else{
      return(M[1])
    }
}
####this function is to calculate the mean,variance and sharpe for user input weights
get_mean_vol_sharpe<-function(pro_dict1=pro_dict_all,w_1,tcost,k,mu=mu1,loss=loss1,asset_cov=asset_cov){ 
  p0<-pro_dict1$p0[1]
  # pro_dict_matrix<-as.matrix(pro_dict1[2:nrow(pro_dict1),])
  pro_dict1<-as.matrix(pro_dict1)
  mu_adjust<-mu1
  P<-mu1
  P1<-mu1
  P2<-mu1
  var_T<-mu1
  var_s<-mu1
  for(idx in 1:length(mu1)){
    mu_adjust[idx]<-(1-k+k*(1-sum(pro_dict1[1,which(pro_dict1[idx+1,]==1)])-p0)/(1-p0))*mu1[idx]+k*sum(pro_dict1[1,which(pro_dict1[idx+1,]==1)])/(1-p0)*loss1[idx]
    P[idx]<-sum(pro_dict1[1,which(pro_dict1[idx+1,]==1)])
    P1[idx]<-k-k/(1-p0)*(1-P[idx]-p0)
    var_s[idx]<-P[idx]*(1-P[idx])
    P2[idx]<-k/(1-p0)*P[idx]
    var_1<-P1[idx]*(1-P1[idx])
    var_2<-P2[idx]*(1-P2[idx])  ##(k-k/(1-p0)*(P[idx]-p0))*(1-(k-k/(1-p0)*(P[idx]-p0)))
    var_T[idx]=(1-P1[idx])^2*asset_cov[idx,idx]+P1[idx]*(1-P1[idx])*(asset_cov[idx,idx]+(mu1[idx])^2)+var_2*(loss1[idx])^2+P2[idx]^2*0-2*(1-P1[idx])*P2[idx]*mu1[idx]*(loss1[idx])
    var_T[idx]=ifelse(var_T[idx]<0,0,var_T[idx])
  }
  #covariance matrix
  cov_adjust<-matrix(0,nrow=length(mu1),ncol=length(mu1))
  E_matrix<-matrix(0,nrow=length(mu1),ncol=length(mu1))
  for(i in 1:length(mu1)){
    cov_adjust[i,i]=var_T[i]
  }
  P_00<-matrix(0,nrow=length(mu1),ncol=length(mu1))
  P_10<-matrix(0,nrow=length(mu1),ncol=length(mu1))
  P_01<-matrix(0,nrow=length(mu1),ncol=length(mu1))
  P_11<-matrix(0,nrow=length(mu1),ncol=length(mu1))
  P1_00<-matrix(0,nrow=length(mu1),ncol=length(mu1))
  P1_10<-matrix(0,nrow=length(mu1),ncol=length(mu1))
  P1_01<-matrix(0,nrow=length(mu1),ncol=length(mu1))
  P1_11<-matrix(0,nrow=length(mu1),ncol=length(mu1))
  P2_00<-matrix(0,nrow=length(mu1),ncol=length(mu1))
  P2_10<-matrix(0,nrow=length(mu1),ncol=length(mu1))
  P2_01<-matrix(0,nrow=length(mu1),ncol=length(mu1))
  P2_11<-matrix(0,nrow=length(mu1),ncol=length(mu1))
  ######initialize matrix
  for(idx1 in 1:length(mu1)){
    for(idx2 in 1:length(mu1)){
      if(idx1!=idx2){
        M_00<-intersect(which(pro_dict1[(idx1+1),]==0),which(pro_dict1[(idx2+1),]==0))
        p_ij_00=ifelse(length(M_00)==0,0,sum(pro_dict1[1,M_00]))
        M_01<-intersect(which(pro_dict1[(idx1+1),]==0),which(pro_dict1[(idx2+1),]==1))
        p_ij_01=ifelse(length(M_01)==0,0,sum(pro_dict1[1,M_01]))
        M_10<-intersect(which(pro_dict1[(idx1+1),]==1),which(pro_dict1[(idx2+1),]==0))
        p_ij_10=ifelse(length(M_10)==0,0,sum(pro_dict1[1,M_10]))
        M_11<-intersect(which(pro_dict1[(idx1+1),]==1),which(pro_dict1[(idx2+1),]==1))
        p_ij_11=ifelse(length(M_11)==0,0,sum(pro_dict1[1,M_11]))
        ##calculate e_ij
        e_ij=asset_cov[idx1,idx2]+mu1[idx1]*mu1[idx2]
        E_matrix[idx1,idx2]=(1-k+k/(1-p0)*(p_ij_00-p0))*e_ij+k/(1-p0)*p_ij_01*mu1[idx1]*loss1[idx2]+k/(1-p0)*p_ij_10*mu1[idx2]*loss1[idx1]+k/(1-p0)*p_ij_11*loss1[idx1]*loss1[idx2] 
        P_00[idx1,idx2]=p_ij_00
        P_10[idx1,idx2]=p_ij_10
        P_01[idx1,idx2]=p_ij_01
        P_11[idx1,idx2]=p_ij_11
        P1_00[idx1,idx2]=k-k/(1-p0)*(1-p_ij_00-p0)
        P1_10[idx1,idx2]=k-k/(1-p0)*(1-p_ij_10-p0)
        P1_01[idx1,idx2]=k-k/(1-p0)*(1-p_ij_01-p0)
        P1_11[idx1,idx2]=k-k/(1-p0)*(1-p_ij_11-p0)
        P2_00[idx1,idx2]=k/(1-p0)*p_ij_00
        P2_10[idx1,idx2]=k/(1-p0)*p_ij_10
        P2_01[idx1,idx2]=k/(1-p0)*p_ij_01
        P2_11[idx1,idx2]=k/(1-p0)*p_ij_11
      }
    }
  }
  for(idx1 in 1:length(mu1)){
    for(idx2 in 1:length(mu1)){
      if(idx1!=idx2){
        e_ij=asset_cov[idx1,idx2]+mu1[idx1]*mu1[idx2]   #
        cov_adjust[idx1,idx2]=(1-k+k/(1-p0)*(P_00[idx1,idx2]-p0))*e_ij+P2_01[idx1,idx2]*mu1[idx1]*loss1[idx2]+P2_10[idx1,idx2]*loss1[idx1]*mu1[idx2]+P2_11[idx1,idx2]*loss1[idx1]*loss1[idx2]-(1-P1[idx1])*(1-P1[idx2])*mu1[idx1]*mu1[idx2]-P2[idx1]*P2[idx2]*loss1[idx1]*loss1[idx2]-(1-P1[idx1])*P2[idx2]*mu1[idx1]*loss1[idx2]-(1-P1[idx2])*P2[idx1]*mu1[idx2]*loss1[idx1] 
      }
    }
  }
  mean_1=mu_adjust%*%w_1+tcost
  sd_1=sqrt(w_1%*%cov_adjust%*%w_1)
  shapre_1=(mean_1-RISK_FREE_RATE)/sd_1
  rett<-c(mean_1,sd_1,shapre_1)
  return(rett)

}
#####this function is a help function to automatically compute finance strategy for a given assets weight
assign_finance_weights<-function(w_1,finan_cost,haircut,real_finance_weight,principal1=1720){
  N<-length(w_1)
  N_fina_con<-ncol(finan_cost)
  flag=0
  w=0
  list1=which(real_finance_weight!=1)
  if(length(list1)!=0){
    #####################################select out assets in short positions,cash nee to meet the margin 
    w=sum(abs(((1-haircut)*real_finance_weight*w_1)[list1]))+sum(abs((haircut*w_1)[list1]))
    if(w_1[length(w_1)]<w){
      flag=1
      w_1[length(w_1)]=w_1[length(w_1)]-w
    }else{
      w_1[length(w_1)]=w_1[length(w_1)]-w
    }
  }
  i=0
  M=MAXIMUM_LOSS
  ####can't consider cash fee
  con2<-colnames(finan_cost)[2:N_fina_con]
  list=which(real_finance_weight!=1)
  # list=unique(c(list,which(real_finance_weight==0)))
  if(length(list)!=0){
    leverage1<-sum(haircut[-list]*w_1[-list]+(1-haircut[-list])*real_finance_weight[-list]*w_1[-list])-(1-w)
  }else{
    leverage1<-sum(haircut*w_1+(1-haircut)*real_finance_weight*w_1)-(1-w)
  }
  finan_cost2=finan_cost[,2:N_fina_con]
  haircut1=(1-haircut)*real_finance_weight*w_1
  w_index=seq(1,N)
  if(leverage1<=0){
    return(rep(0,N))
  }else if(leverage1>sum(haircut1[intersect(which(w_1>0),which(real_finance_weight!=0))])){
    return(rep(-1,N)) 
  }else{
    finan_cost2=finan_cost
    w_index=seq(1,N)
    haircut1<-haircut1
    finan_cost3<-finan_cost2[,2:N_fina_con]
    w_init<-rep(0,length(haircut1))
    lower_b<-rep(0,length(haircut1))
    upper_b<-ifelse(haircut1>0,haircut1,0)
    mm<-length(haircut1)
    haircut1<-haircut1
    if(length(con2)==1){ 
      lower_b<-rep(0,mm)
      upper_b<-ifelse(haircut1>0,haircut1,0)
      f.obj<-finan_cost3
      f.con<-rbind(diag(mm),diag(mm),rep(1,mm))
      f.dir<-c(rep("<=",mm),rep(">=",mm),"==")
      f.rhs<-c(upper_b,lower_b,leverage1)
      ress<-lp("max", f.obj, f.con, f.dir, f.rhs)$solution   ###convex piece wise linear
      return(ress)
    }else if(length(con2)!=1){
      finan_cost3<--finan_cost3
      lower_b<-rep(0,mm)
      upper_b<-ifelse(haircut1>0,haircut1,0)
      f.obj<-c(rep(0,mm),rep(1,mm))
      f.con<-rbind(cbind(diag(mm),matrix(0,mm,mm)),cbind(diag(mm),matrix(0,mm,mm)),c(rep(1,mm),rep(0,mm)))
      for(kk in 1:mm){
        temp<-matrix(0,nrow=length(con2),ncol=2*mm)
        for(jj in 1:length(con2)){
          temp[jj,kk]<-finan_cost3[kk,jj]
          temp[jj,mm+kk]<--1
        }
        f.con<-rbind(f.con,temp)
      }
      f.dir<-c(rep("<=",mm),rep(">=",mm),"==",rep("<=",length(con2)*mm))
      f.rhs<-c(upper_b,lower_b,leverage1)
      for(kk in 1:mm){
        temp<-c(0)
        for(jj in 2:length(con2)){
          range<-get_number(con2[jj])
          if(length(range)==2){
            range1<-(range[2]-range[1])/principal1
          }else{
            range1<-range[1]/principal1
          } 
          temp<-c(temp,temp[length(temp)]+range1*(finan_cost3[kk,jj]-finan_cost3[kk,jj-1]))
        }
        f.rhs<-c(f.rhs,temp)
      }
      ress<-lp("min", f.obj, f.con, f.dir, f.rhs)  ###convex piece wise linear
      ress2<-ress$solution[1:mm]
      return(ress2)
    }
  }
}
calculate_total_weights<-function(w_1,finan_cost,haircut,real_finance_weight){
  list1=which(real_finance_weight!=1)
  ifelse(length(list1)!=0,sum(w_1[-list1]),sum(w_1))
}

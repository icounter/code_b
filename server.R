# install.packages(c("ggplot2","gRbase","MASS","gRain","copula","psych","nloptr","rpsychi","plyr","shinysky","DT","rootSolve","stringr","bnlearn","lpSolve","shinydashboard","data.table"))
#source("http://bioconductor.org/biocLite.R"); biocLite(c("graph","RBGL","Rgraphviz"))
# install.packages("devtools")
# library(devtools)
# install_github("shinyAce", "trestletech")
# devtools::install_github("AnalytixWare/ShinySky")
library(ggplot2)
library(MASS)
#library(gRbase)
library(gRain)
#library(LambertW)
library(copula)
library(psych)
library(nloptr)
library(rpsychi)
library(plyr)
library(shinysky)
library(DT)
library(rootSolve)
library(stringr)
library(graph)
library(Rgraphviz)
library(bnlearn)
library(lpSolve)
library(shinydashboard)
library(data.table)
# library(Rsolnp)
# library(parallel)
source("script.R")
shinyServer(
  function(input, output,session){
  output$contents <- renderTable({
  inFile <- input$file1
  if(is.null(inFile)){
    return(NULL)
  }else{
    input_data<<-read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote)
  }})
##########################################################################################################
############################################read in file&update###################################
##########################################################################################################
  observe({
     if(!is.null(input$file1)){
       NN<<-nrow(input_data)
       updateNumericInput(session,"assets_num1",value=nrow(input_data))
       updateNumericInput(session,"assets_num",value=nrow(input_data))
     name<-as.character(input_data[,1])
     for(i in 1:length(name)){
       name[i]<-paste0(strsplit(name[i]," ")[[1]],collapse = "")
     }
     updateTextInput(session,"asset_name1",value=paste0(name,collapse = ","))
     updateTextInput(session,"asset_vol",value=paste0(input_data[,5],collapse = ","))
     updateTextInput(session,"asset_vol1",value=paste0(input_data[,5],collapse = ","))
     updateTextInput(session,"asset_ret1",value=paste0(input_data[,3],collapse = ","))
     updateTextInput(session,"asset_ret",value=paste0(input_data[,3],collapse = ","))
     updateTextInput(session,"extreme_stress_loss",value=paste0(input_data[,4],collapse = ","))
     N<-nrow(input_data)
     a=c()
     for( i in 1:(N-1)){
       a<-c(a,as.numeric(input_data[i,(6+i):(5+N)]))
     }
     updateTextInput(session,"asset_corr1",value=paste0(a,collapse = ","))
     updateTextInput(session,"asset_corr",value=paste0(a,collapse = ","))
     updateTextInput(session,"w_now",value=paste0(input_data[,2],collapse = ","))
     updateTextInput(session,"w1",value=paste0(input_data[,2],collapse = ","))
     updateTextInput(session,"lower_bound",value=paste0(input_data[,which(colnames(input_data)=="lower_bounds")],collapse = ","))
     updateTextInput(session,"lower_bound1",value=paste0(input_data[,which(colnames(input_data)=="lower_bounds")],collapse = ","))
     updateTextInput(session,"upper_bound",value=paste0(input_data[,which(colnames(input_data)=="upper_bounds")],collapse = ","))
     updateTextInput(session,"upper_bound1",value=paste0(input_data[,which(colnames(input_data)=="upper_bounds")],collapse = ","))
     bn<-as.character(input_data[,which(colnames(input_data)=='BN')])
     bn_hidden<-as.character(input_data[,which(colnames(input_data)=='BN_hidden')])
     BN<-paste0(bn_hidden[which(nchar(bn_hidden)!=0)],collapse = ";")
     BN2<-paste0(name[which(nchar(bn)!=0)],"+",bn[which(nchar(bn)!=0)])
     BN<-paste0(BN,";",paste0(BN2,collapse = ";"),";")
     BN<-paste0(BN,name[which(nchar(bn)==0)],collapse = ";")
     updateTextInput(session,"affect_relation",value=BN)
     con=as.character(input_data[,which(colnames(input_data)=='cond')])
     con_hidden=as.character(input_data[,which(colnames(input_data)=='cond_Hidden')])
     condd<-paste0(con_hidden[which(nchar(con_hidden)!=0)],collapse=";")
     condd<-paste0(condd,";",paste0(con[which(nchar(bn)!=0)],collapse = ";"),";")
     condd<-paste0(condd,con[which(nchar(bn)==0)],collapse = ";")
     updateTextInput(session,"prob_table",value=condd)
       trans_generaed<<-1
       N1<-which(colnames(input_data)=="haircut")
       N4<-which(colnames(input_data)=="transaction_cost")
       dddd<-input_data[,(N4+1):(N1-1)]
       N6<-which(colnames(dddd)=="finance_cost")
       dddd<-dddd[,-N6]
       N3<-(ncol(dddd)+1)-N6
       N2<-(ncol(dddd))-N3
       trans_cost1<-data.frame((matrix(0,nrow=N,ncol=N2+1)))
       finan_cost1<-data.frame((matrix(0,nrow=N,ncol=N3+1)))
       if(N2==1){
         string1<-rep(name,1)
         string1=paste0(string1,"_trans_cost")
         trans_cost1[,1]<-string1
         trans_cost1[,2]=dddd[,1]
         colnames(trans_cost1)<-c("Asset_name","!=0")
       }else{
         string1<-rep(name,1)
         string1=paste0(string1,"_trans_cost")
         trans_cost1[,1]<-string1
         trans_cost1[,2:(N2+1)]=dddd[,1:N2]
         eq<-colnames(dddd)[1:N2]
         for(i in 1:N2){
           eq[i]=str_replace_all(eq[i],"and","&")
           eq[i]=str_replace_all(eq[i],"s","<")
           eq[i]=str_replace_all(eq[i],"l",">")
           eq[i]=str_replace_all(eq[i],"t","")
           eq[i]=str_replace_all(eq[i],"e","=")
           eq[i]=str_replace_all(eq[i],"ineq","!=")
         }
         colnames(trans_cost1)<-c("Asset_name",eq)
       }
       if(N3==1){
         string1<-rep(name,1)
         string1=paste0(string1,"_finance_cost")
         finan_cost1[,1]<-string1
         finan_cost1[,2]=dddd[,(N2+1)]
         colnames(finan_cost1)<-c("Asset_name","!=0")
       }else{
         string1<-rep(name,1)
         string1=paste0(string1,"_finance_cost")
         finan_cost1[,1]<-string1
         finan_cost1[,2:(N3+1)]=dddd[,(N2+1):(N3+N2)]
         eq<-colnames(dddd)[(N2+1):(N2+N3)]
         for(i in 1:N3){
           eq[i]=str_replace_all(eq[i],"and","&")
           eq[i]=str_replace_all(eq[i],"s","<")
           eq[i]=str_replace_all(eq[i],"l",">")
           eq[i]=str_replace_all(eq[i],"f","")
           eq[i]=str_replace_all(eq[i],"e","=")
           eq[i]=str_replace_all(eq[i],"ineq","!=")
         }
         colnames(finan_cost1)<-c("Asset_name",eq)
       }
       haircut<<-input_data[,which(colnames(input_data)=="haircut")]
       real_finance_weight<<-input_data[,which(colnames(input_data)=="real_finance_weight")]
       trans_cost1<<-trans_cost1
       finan_cost1<<-finan_cost1}   
   })
########################################################################################################################
##########################################out plot weights## ####################################################
########################################################################################################################
  weights1<-eventReactive(input$go1,{
   call_scenario_no_cost_no_BN(input$assets_num,input$lambda,input$asset_ret,input$asset_vol,input$asset_corr,input$sample_number,input$lower_bound1,input$upper_bound1)
  })#####this function is to get weights
  output$ggplot <- renderPlot({
  weights<-weights1()
  ggplot(data=weights, aes(x=asset_name, y=weight, fill=methods))+geom_bar(stat="identity", position=position_dodge())+geom_text(aes(label=round(weight,digits = 4)), vjust=1.6,position = position_dodge(0.9), size=5)
                             })
#####################################################out put rescaled probability matrix  
con_tabe<-eventReactive(input$go_BN_cond,{
      bayesian_matrix1<-data.frame(t(bayesian_matrix1))
      rescale<-rep(0,nrow(bayesian_matrix1))
      if(bayesian_matrix1[1,1]!=1){
      rescale[1]<-1-input$tentative_k
      rescale[2:length(rescale)]<-bayesian_matrix1[2:nrow(bayesian_matrix1),1]/(1-bayesian_matrix1[1,1])*input$tentative_k
      }
      bayesian_matrix1<-cbind(bayesian_matrix1,rescale)
      colnames(bayesian_matrix1)[ncol(bayesian_matrix1)]<-"Rescaled prob"
      bayesian_matrix1<-bayesian_matrix1[,c(c(1,ncol(bayesian_matrix1)),seq(2,ncol(bayesian_matrix1)-1))]
      bayesian_matrix1[,1]<-sprintf("%.6f",bayesian_matrix1[,1])
      bayesian_matrix1[,2]<-sprintf("%.6f",bayesian_matrix1[,2])
      bayesian_matrix1[,3:ncol(bayesian_matrix1)]<-as.data.frame(lapply(bayesian_matrix1[,3:ncol(bayesian_matrix1)], format_num))
      return(bayesian_matrix1)
    })
#####################################################conditional probability matrixA happens B happens probability 
generate_bayesian_cor<-eventReactive(input$go_BN_cond,{
      pro_dict3<-as.matrix(pro_dict)
      pro_dict3<-pro_dict3[,-1]
      if(is.null(nrow(pro_dict3))){
        pro_dict_prob<-pro_dict3[1]
        pro_dict_matrix<-as.matrix(pro_dict3[2:length(pro_dict3)])
      }else{
        pro_dict_prob<-pro_dict3[1,]
        pro_dict_matrix<-pro_dict3[2:nrow(pro_dict3),]
      }
      name<-rownames(pro_dict_matrix)
      N<-length(name)
      cor<-diag(N)
      if(ncol(pro_dict_matrix)==1){
       cor<-matrix(1,N,N)
       cor[,N]=0
       cor[N,]=0
      }else{
      for(i in 1:N){
        for(j in 1:N){
          if(i==j) {
            cor[i,j]=1
          }else{
           temp_pro<-which(pro_dict_matrix[i,]==1)
           if(length(temp_pro)==0){
             cor[i,j]=0
           }else{
             temp<-pro_dict_matrix[,temp_pro]
             temp_pro_sum<-pro_dict_prob[temp_pro]
             temp_pro2<-which(temp[j,]==1)
             if(length(temp_pro2)==0){
               cor[i,j]=0
             }else{
               cor[i,j]=(sum(temp_pro_sum[temp_pro2]))/sum(temp_pro_sum)
             }
           } 
          }
        }
      }}
      cor<-data.frame(cor)
      colnames(cor)<-name
      cor<-as.data.frame(lapply(cor, format_num9))
      cor<-cbind(name,cor)
      colnames(cor)[1]<-c("conditional probability matrix")
      return(cor) })
#####################################################correlation matrix based on simulation   
generate_extreme_cor<-eventReactive(input$go_BN_cond,{
      pro_dict3<-as.matrix(pro_dict)
      pro_dict3<-pro_dict3[,-1]
      pro_dict_prob<-pro_dict3[1,]
      pro_dict_matrix<-pro_dict3[2:nrow(pro_dict3),]
      M<-50000
      extreme_stress_loss<-as.double(strsplit(input$extreme_stress_loss,",")[[1]])
      name<-rownames(pro_dict_matrix)
      N<-length(name)
      asset_ret<-as.double(strsplit(input$asset_ret1,",")[[1]])
      asset_corr<-as.double(strsplit(input$asset_corr1,",")[[1]])
      asset_vol1<-as.double(strsplit(input$asset_vol1,",")[[1]])
      asset_var<-(asset_vol1)^2
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
      margins_r<-rep("norm",N)
      paramMargins_r <- list()
      for(i in 1:N){
        paramMargins_r[[length(paramMargins_r)+1]] <- list(mean =asset_ret[i], sd =sqrt(asset_var[i]))
      }
      randl<<-generate_N_rand(N,asset_corr,margins_r,paramMargins_r,M)
      rand<-randl
      order<-as.integer(M*pro_dict_prob/sum(pro_dict_prob))
      for(i in 1:length(order)){
        if(i==1&&order[i]!=0){
          rand[1:order[i],]<-t(((1-pro_dict_matrix[,i])*t(rand[1:order[i],])))
          rand[1:order[i],]<-t(apply(rand[1:order[i],],1,function(x) x+as.vector(pro_dict_matrix[,i]*extreme_stress_loss)))
        }else if(order[i]!=0){
          rand[(sum(order[1:(i-1)])+1):(sum(order[1:i])),]<-t(((1-pro_dict_matrix[,i])*t(rand[(sum(order[1:(i-1)])+1):(sum(order[1:i])),])))
          if((sum(order[1:(i-1)])+1)==(sum(order[1:i]))){
            rand[(sum(order[1:(i-1)])+1):(sum(order[1:i])),]<-rand[(sum(order[1:(i-1)])+1):(sum(order[1:i])),]+as.vector(pro_dict_matrix[,i]*extreme_stress_loss)
          }else{
            rand[(sum(order[1:(i-1)])+1):(sum(order[1:i])),]<-t(apply(rand[(sum(order[1:(i-1)])+1):(sum(order[1:i])),],1,function(x) x+as.vector(pro_dict_matrix[,i]*extreme_stress_loss)))
          }
        }
      }
      pro_dict_matrix<-rand
      cor<-diag(N)
      cor<-cor(pro_dict_matrix)
      for(i in 1:ncol(cor)){
        for(j in 1:nrow(cor)){
          if(is.na(cor[i,j]))  cor[i,j]=0
        }
      }
      cor<-data.frame(cor)
      colnames(cor)<-name
      cor<-as.data.frame(lapply(cor, format_num9))
      cor<-cbind(name,cor)
      colnames(cor)[1]<-c("correlation matrix")
      return(cor)
    })
#####################################################rescaled correlation matrix based on simulation       
    generate_bayesian_cor2<-eventReactive(input$go_BN_cond,{
      pro_dict3<-as.matrix(t(bayesian_matrix1))
      rescale<-rep(0,nrow(pro_dict3))
      if(pro_dict3[1,1]!=1){
        rescale[1]<-1-input$tentative_k
        rescale[2:length(rescale)]<-pro_dict3[2:nrow(pro_dict3),1]/(1-pro_dict3[1,1])*input$tentative_k
      }
      pro_dict3[,1]<-rescale
      pro_dict3<-t(pro_dict3)
      pro_dict_prob<-pro_dict3[1,]
      pro_dict_matrix<-pro_dict3[2:nrow(pro_dict3),]
      name<-rownames(pro_dict_matrix)
      N<-length(name)
      extreme_stress_loss<-as.double(strsplit(input$extreme_stress_loss,",")[[1]])
       M<-50000
       name<-rownames(pro_dict_matrix)
       N<-length(name)
       if(exists("randl")&&ncol(randl)==N){
         rand<-randl
       }else{
         name<-rownames(pro_dict_matrix)
         N<-length(name)
         asset_ret<-as.double(strsplit(input$asset_ret1,",")[[1]])
         asset_corr<-as.double(strsplit(input$asset_corr1,",")[[1]])
         asset_vol1<-as.double(strsplit(input$asset_vol1,",")[[1]])
         asset_var<-(asset_vol1)^2
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
         margins_r<-rep("norm",N)
         paramMargins_r <- list()
         extreme_stress_loss<-as.double(strsplit(input$extreme_stress_loss,",")[[1]])
         for(i in 1:N){
           paramMargins_r[[length(paramMargins_r)+1]] <- list(mean =asset_ret[i], sd =sqrt(asset_var[i]))
         }
         rand<-generate_N_rand(N,asset_corr,margins_r,paramMargins_r,M)
       }
      order<-as.integer(M*pro_dict_prob/sum(pro_dict_prob))
      pro_dict_matrix<-pro_dict_matrix
      for(i in 1:length(order)){
        if(i==1&&order[i]!=0){
          rand[1:order[i],]<-t(((1-pro_dict_matrix[,i])*t(rand[1:order[i],])))
          rand[1:order[i],]<-t(apply(rand[1:order[i],],1,function(x) x+as.vector(pro_dict_matrix[,i]*extreme_stress_loss)))
        }else if(order[i]!=0){
          rand[(sum(order[1:(i-1)])+1):(sum(order[1:i])),]<-t(((1-pro_dict_matrix[,i])*t(rand[(sum(order[1:(i-1)])+1):(sum(order[1:i])),])))
          if((sum(order[1:(i-1)])+1)==(sum(order[1:i]))){
            rand[(sum(order[1:(i-1)])+1):(sum(order[1:i])),]<-rand[(sum(order[1:(i-1)])+1):(sum(order[1:i])),]+as.vector(pro_dict_matrix[,i]*extreme_stress_loss)
          }else{
            rand[(sum(order[1:(i-1)])+1):(sum(order[1:i])),]<-t(apply(rand[(sum(order[1:(i-1)])+1):(sum(order[1:i])),],1,function(x) x+as.vector(pro_dict_matrix[,i]*extreme_stress_loss)))
          }
        }
      }
      pro_dict_matrix<-rand
      cor<-cor(pro_dict_matrix)
      for(i in 1:ncol(cor)){
        for(j in 1:nrow(cor)){
          if(is.na(cor[i,j]))  cor[i,j]=0
        }
      }
      cor<-data.frame(cor)
      colnames(cor)<-name
      cor<-as.data.frame(lapply(cor, format_num9))
      cor<-cbind(name,cor)
      colnames(cor)[1]<-c("Rescaled correlation matrix")
      return(cor)
    })
###############################################normal correlation matrix
    generate_normal_cor<-eventReactive(input$go_BN_cond,{
      N<-as.double(input$assets_num1)
      asset_corr1<-as.double(strsplit(input$asset_corr1,",")[[1]])
      asset_cor<-matrix(1,N,N)
      kk=1
      for(i in 1:(N-1)){
        for(j in (i+1):(N)){
          asset_cor[i,j]=asset_corr1[kk]
          kk=kk+1
        }
      }
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          asset_cor[j,i]=asset_cor[i,j]
        }
      }
      name<-strsplit(input$asset_name1,",")[[1]]
      asset_cor<-data.frame(asset_cor)
      colnames(asset_cor)<-name
      asset_cor[,1:ncol(asset_cor)]<-as.data.frame(lapply(asset_cor[,1:ncol(asset_cor)], format_num9))
      asset_cor<-cbind(name,asset_cor)
      colnames(asset_cor)[1]<-c("Normal part correlation matrix")
      return(asset_cor)
    })
############################################downside correlation matrix for extreme part
    generate_bayesian_cor_downside<-eventReactive(input$go_BN_cond,{
      pro_dict3<-as.matrix(t(bayesian_matrix1))
      pro_dict3<-as.matrix(pro_dict)
      pro_dict3<-pro_dict3[,-1]
      pro_dict_prob<-pro_dict3[1,]
      pro_dict_prob<-pro_dict_prob/sum(pro_dict_prob)
      pro_dict_matrix<-pro_dict3[2:nrow(pro_dict3),]
      extreme_stress_loss<-as.double(strsplit(input$extreme_stress_loss,",")[[1]])
      pro_dict_matrix<-(pro_dict_matrix)*extreme_stress_loss
      name<-rownames(pro_dict_matrix)
      N<-length(name)
      cor<-diag(N)
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          D1<-sqrt(pro_dict_prob%*%(pro_dict_matrix[i,]*pro_dict_matrix[i,]))
          D2<-sqrt(pro_dict_prob%*%(pro_dict_matrix[j,]*pro_dict_matrix[j,]))
          cor[i,j]<-(pro_dict_prob%*%(pro_dict_matrix[i,]*pro_dict_matrix[j,]))/(D1*D2)
          cor[i,j]=ifelse(is.na(cor[i,j]),0,cor[i,j])
          cor[j,i]=cor[i,j]
        }
        
      }
      cor<-data.frame(cor)
      colnames(cor)<-name
      cor<-as.data.frame(lapply(cor, format_num9))
      cor<-cbind(name,cor)
      colnames(cor)[1]<-c("downside correlation matrix for extreme part")
      return(cor)
    })
#################################format functions
    format_num <- function(col) {
      if (is.numeric(col))
        sprintf('%1.0f', col)
      else
        col
    }
    format_num2 <- function(col) {
      if (is.numeric(col))
        sprintf('%.3f', col)
      else
        col
    }
    format_num8 <- function(col) {
      if (is.numeric(col))
        sprintf('%.2f', 10000*col)
      else
        col
    }
    format_num9 <- function(col) {
      if (is.numeric(col))
        sprintf('%.2f%%', 100*col)
      else
        col
    }
    format_num10 <- function(col) {
      if (is.numeric(col))
        sprintf('%.2f', col)
      else
        col
    }
#######################################out put data tables
    output$hotable2 <- DT::renderDataTable(
      DT::datatable(con_tabe(), options = list(
        lengthMenu = list(c(5, 15, -1), c('5', '15', 'All')),
        pageLength = 15
    )))
    output$bayesian_matrix_cor_normal <- DT::renderDataTable(
      DT::datatable(generate_normal_cor(), options = list(
        lengthMenu = list(c(5, 15, -1), c('5', '15', 'All')),
        pageLength = 15
      )))
    output$bayesian_matrix_cor_extreme <- DT::renderDataTable(
      DT::datatable(generate_extreme_cor(), options = list(
        lengthMenu = list(c(5, 15, -1), c('5', '15', 'All')),
        pageLength = 15
      )))
    output$bayesian_matrix_cor <- DT::renderDataTable(
      DT::datatable(generate_bayesian_cor(), options = list(
        lengthMenu = list(c(5, 15, -1), c('5', '15', 'All')),
        pageLength = 15
      )))
    output$bayesian_matrix_cor_downside <- DT::renderDataTable(
      DT::datatable(generate_bayesian_cor_downside(), options = list(
        lengthMenu = list(c(5, 15, -1), c('5', '15', 'All')),
        pageLength = 15
      )))
    output$bayesian_matrix_cor2 <- DT::renderDataTable(
      DT::datatable(generate_bayesian_cor2(), options = list(
        lengthMenu = list(c(5, 15, -1), c('5', '15', 'All')),
        pageLength = 15
      )))
#################################plot baeysian network
    plot_bay<-eventReactive(input$go_plot_bayesian_net,{
      cond_tab=build_cond(affect_relation=input$affect_relation,prob_table=input$prob_table)   ###build condition probability table
      plist<-compileCPT(cond_tab)
      net<-grain(plist) 
      bn_net<-as.bn.fit(net)
      asset_name1=strsplit(input$asset_name1,",")[[1]]
      bayesian_matrix1<<-bayesian_matrix(cond_tab,asset_name=asset_name1)
      graphviz.plot(bn_net,layout='dot',shape='ellipse')
    })
    output$plot_bayesian<-renderPlot({return(plot_bay())})
############################generate finance cost matrix
    finan_cost2<-eventReactive(input$go_generate_finan,{
      if((input$assets_num1==4)&(trans_generaed==0)){
        condtion_name2='!=0'
        con_number2=1
        finan_cost1<-data.frame((matrix(0,nrow=input$assets_num1,ncol=con_number2+1)))
        col_names<-strsplit(condtion_name2,",")[[1]]
        colnames(finan_cost1)<-c("Asset_name",col_names)
        asset_name<-strsplit(input$asset_name1,",")[[1]]
        string2<-rep(asset_name,1)
        string2[1:input$assets_num1]=paste0(string2[1:input$assets_num1],"_finance_cost")
        finan_cost1[,1]<-string2
        finan_cost<<-finan_cost1
        return(finan_cost)
      }else{
      finan_cost<<-(finan_cost1)
      return(finan_cost)}
      }) 
    output$hotable5 <- DT::renderDataTable(
      DT::datatable(finan_cost2(), options = list(
        lengthMenu = list(c(5, 15, -1), c('5', '15', 'All')),
        pageLength = 15
      )))
    trans_cost2<-eventReactive(input$go_generate_trans,{
     if((input$assets_num1==4)&(trans_generaed==0)){
       con_number=1
       condtion_name="!=0"
       trans_cost1<-data.frame((matrix(0,nrow=input$assets_num1,ncol=con_number+1)))
       col_names<-strsplit(condtion_name,",")[[1]]
       colnames(trans_cost1)<-c("Asset_name",col_names)
       asset_name<-strsplit(input$asset_name1,",")[[1]]
       string1<-rep(asset_name,1)
       string1[1:input$assets_num1]=paste0(string1[1:input$assets_num1],"_trans_cost")
       trans_cost1[,1]<-string1
       trans_cost<<-(trans_cost1)
       return(trans_cost)
     }else{
     trans_cost<<-trans_cost1
     return(trans_cost)
     }
    })
    output$hotable1 <- DT::renderDataTable(
      DT::datatable(trans_cost2(),options = list(
        lengthMenu = list(c(5, 15, -1), c('5', '15', 'All')),
        pageLength = 15)))
########################################################################################generate tentative cost matrix
tentative_trans<-eventReactive(input$go_tentative_cost,{
      w_now<-as.double(strsplit(input$w_now,",")[[1]])
      w1<-as.double(strsplit(input$w1,",")[[1]])
      principal1<-input$principal1
      cost2(w_now,w1,trans_cost,finan_cost,haircut,real_finance_weight,principal1)
    })
output$tentative_transcost<-renderTable({ 
      weig<-tentative_trans()
      colnames(weig)<-c("trans_cost(bps)","finance_cost(bps)","finance_weights","haircut","allocated_weights","initial_Weights","weights_changed")
      weig<-data.frame(weig)
      name<-strsplit(input$asset_name1,",")[[1]]
      weig<-cbind(name,weig)
      colnames(weig)[1]="name"
      weig[,2:(ncol(weig)-5)]<-as.data.frame(lapply(weig[,2:(ncol(weig)-5)], format_num8))
      weig[,(ncol(weig)-4):(ncol(weig))]<-as.data.frame(lapply(weig[,(ncol(weig)-4):(ncol(weig))], format_num9))
      return(weig)
    }, readOnly = TRUE)
#########################################################the left plot when choosing lambda
 output$plot_ut_zero_equal<-renderPlot({
      ret<-seq(from=0,to=0.15,by=0.01)
      lam<-c(input$lambda1)
      gene<-data.frame(matrix(0,nrow=(length(lam)*length(ret)),ncol=3))
      colnames(gene)<-c("lambda","ret","probability")
      if(input$uti=='expo'){
        for(i in 1:length(lam)){
          for(j in 1:length(ret)){
            if(j==1){
              gene[((i-1)*length(ret)+j),]=c(lam[i],ret[j],0.5)
            }else{
              b=(1-exp(-lam[i]*ret[j]))/(exp(lam[i]*ret[j])-exp(-lam[i]*ret[j]))
              gene[((i-1)*length(ret)+j),]=c(lam[i],ret[j],b)
            }
          }
        }
      }else if(input$uti=='power'){
        for(i in 1:length(lam)){
          for(j in 1:length(ret)){
            if(j==1){
              gene[((i-1)*length(ret)+j),]=c(lam[i],ret[j],0.5)
            }else{
              b=((1+ret[j])^(1-lam[i])-1)/((1+ret[j])^(1-lam[i])-(1-ret[j])^(1-lam[i]))
              gene[((i-1)*length(ret)+j),]=c(lam[i],ret[j],b)
            }
          }
        }
        gene=na.omit(gene)
      }else{
        for(i in 1:length(lam)){
          for(j in 1:length(ret)){
            if(j==1){
              gene[((i-1)*length(ret)+j),]=c(lam[i],ret[j],0.5)
            }else{
              b=log(1+lam[i]*ret[j])/(log(1+lam[i]*ret[j])-log(1-lam[i]*ret[j]))
              gene[((i-1)*length(ret)+j),]=c(lam[i],ret[j],b)
            }
          }
        }
        gene=na.omit(gene)
      }
      gene[,3]<-format_num2(gene[,3])
      gene$probability<-1-as.double(gene$probability)
      ggplot(data=gene,aes(x=ret,y=probability,group=lambda,colour=lambda))+geom_line()+geom_text(aes(label=probability), size = 4,check_overlap = TRUE)+scale_y_continuous(breaks = round(seq(min(gene$probability), max(gene$probability), by = 0.03),5))+scale_x_continuous(breaks = round(seq(min(gene$ret), max(gene$ret), by = 0.01),5)) +xlab("Three-Month Return +/-")+ylab("probability of Win")
    })
 #########################################################the right plot when choosing lambda
    output$plot_ut_zero_equal_2<-renderPlot({
      neg_ret<-seq(from=0,to=input$neg_ret_ran,by=-0.01)
      lam<-c(input$lambda1)
      gene<-data.frame(matrix(0,nrow=(length(lam)*length(neg_ret)),ncol=3))
      colnames(gene)<-c("lambda","neg_ret","pos_ret")
      if(input$uti=='expo'){
        for(i in 1:length(lam)){
          for(j in 1:length(neg_ret)){
            if(j==1){
              gene[((i-1)*length(neg_ret)+j),]=c(lam[i],neg_ret[j],0)
            }else{
              b=log(2-exp(-lam[i]*neg_ret[j]))/(-lam[i])
              gene[((i-1)*length(neg_ret)+j),]=c(lam[i],neg_ret[j],b)
            }
          }
        }
        gene=na.omit(gene)
      }else if(input$uti=='power'){
        for(i in 1:length(lam)){
          for(j in 1:length(neg_ret)){
            if(j==1){
              gene[((i-1)*length(neg_ret)+j),]=c(lam[i],neg_ret[j],0)
            }else{
              b=(2-(1+neg_ret[j])^(1-lam[i]))^(1/(1-lam[i]))-1
              gene[((i-1)*length(neg_ret)+j),]=c(lam[i],neg_ret[j],b)
            }
          }
        }
        gene=na.omit(gene)
      }else{
        for(i in 1:length(lam)){
          for(j in 1:length(neg_ret)){
            if(j==1){
              gene[((i-1)*length(neg_ret)+j),]=c(lam[i],neg_ret[j],0)
            }else{
              b=(1/(1+lam[i]*neg_ret[j])-1)/lam[i]
              gene[((i-1)*length(neg_ret)+j),]=c(lam[i],neg_ret[j],b)
            }
          }
        }
        gene=na.omit(gene)
      }
      gene[,3]<-format_num2(gene[,3])
      gene$pos_ret<-as.double(gene$pos_ret)
      gene$neg_ret<-as.double(gene$neg_ret)
      ggplot(data=gene,aes(x=neg_ret,y=pos_ret,group=lambda,colour=lambda))+geom_line()+geom_text(aes(label=pos_ret), size = 4,check_overlap = TRUE)+xlab("Return in Bad Scenario")+ylab("Return in good Scenario")+scale_x_reverse(breaks = round(seq(from=min(as.double(gene$neg_ret)), to=max(as.double(gene$neg_ret)), by = 0.01),100))+scale_y_continuous(breaks = round(seq(from=min(as.double(gene$pos_ret)), to=max(as.double(gene$pos_ret)), by = 0.03),100))
    })
    ##############################out plot combo utility
    output$utility_out2<-renderPlot({
      res11<-utility_solve(x_extreme =input$x_extreme,x_downturning =input$x_downturning,0,x_upturning = input$x_upturning,prob2 = input$prob2)
      ggplot(data.frame(x=input$x_range2),aes(x))+stat_function(fun=vutility,args=list(x_extreme1=x_extreme,x_downturning1=x_downturning,k_21=k_2,k_11=k_1,x_11=x_1),geom='line',aes(colour="combination utility"))+xlab("return")+ylab("utility value")
    })
    ##############################out plot other three utility
    output$utility_out<-renderPlot({
      if(input$uti=='log'){
        ggplot(data.frame(x=input$x_range),aes(x))+stat_function(fun=function(x)log(1+(input$lambda1*x)),geom='line',aes(colour="log utility"))+xlab("return")+ylab("utility value")
      }else if(input$uti=='expo'){
        ggplot(data.frame(x=input$x_range),aes(x))+stat_function(fun=function(x)((1-exp(-input$lambda1*x))/input$lambda1),geom='line',aes(colour="exponential utility"))+xlab("return")+ylab("utility value")
      }else{
        ggplot(data.frame(x=input$x_range),aes(x))+stat_function(fun=function(x)(1/(1-input$lambda1)*((1+x)^(1-input$lambda1)-1)),geom='line',aes(colour="power utility"))+xlab("return")+ylab("utility value")
      }
    })
    weights2<-eventReactive(input$go_start_compute,{
      recalculate<<-0
      if(input$uti!="combo"){
        call_scenario_cost_BN(input$uti,input$assets_num1,input$lambda1,input$asset_ret1,input$asset_vol1,
                            input$asset_corr1,input$sample_number1,input$extreme_stress_loss,pro_dict,
                            input$principal1,trans_cost,finan_cost,haircut,real_finance_weight,input$w_now,input$lower_bound,input$upper_bound,input$subjective_k,
                            input$inequality_constraints,input$equality_constraints,input$maxeval_global,input$maxeval_local,input$asset_name1,input$moment,smallprob)}else{
        call_scenario_cost_BN2(input$assets_num1,input$x_extreme,input$x_downturning,0,input$x_upturning,input$prob2,input$asset_ret1,input$asset_vol1,
                                                    input$asset_corr1,input$sample_number1,input$extreme_stress_loss,pro_dict,
                                                    input$principal1,trans_cost,finan_cost,haircut,real_finance_weight,input$w_now,input$lower_bound,input$upper_bound,input$subjective_k,
                                                    input$inequality_constraints,input$equality_constraints,input$maxeval_global,input$maxeval_local,input$asset_name1,input$moment,smallprob)
      }
      })
    output$opt_weights<-renderPlot({
      weights3<<-weights2()
      if(input$uti=='expo'){
        u<--expo_find_w(weights3$weights[-length(weights3$weights)])
      }else if(input$uti=='power'){
        u<--power_find_w(weights3$weights[-length(weights3$weights)])
      }else if(input$uti=='log'){
        u<--log_find_w(weights3$weights[-length(weights3$weights)])
      }else{
        u<--combo_find_w(weights3$weights[-length(weights3$weights)],w_now,x_downturning1,AA1,ll1,k_21,k_11,x_11,trans_cost,finan_cost,haircut,real_finance_weight,principal1,rand2,rand_sub,rand2_moment,loss1,pro_dict,k,mu)
      }
      weights3$asset_name<-as.character(weights3$asset_name)
      ggplot(data=weights3, aes(x=asset_name, y=weights, fill=method)) +
        geom_bar(stat="identity", position=position_dodge())+
        geom_text(aes(label=round(weights,digits = 3),x=asset_name), vjust=1.6, size=5)+annotate(geom = "text", x = length(weights3$weights), y = max(weights3$weights),
                                                                                                 label = paste0("utility value:",sprintf("%.4f",u)),
                                                                                                 hjust = 1)
    })
#######how to caluculate stats for matrix
#####this is based on simulation
######basically devide the samples into 3 parts:assets,optimized total assets and cash
##########simulation 50000 paths for normal scenario and based on k we add some more paths as extreme paths and this extreme path is not relied on simualtion
################################then calculate the mean,vol,sharpe based on total paths
####################same idea when we plot CDF
calculate_stat<-function(weights3){ 
      w_1<-weights3$weights[-length(weights3$weights)] 
      name<-as.character(weights3$asset_name[-length(weights3$weights)])
      tcost<-cost(w_now,w_1,trans_cost,finan_cost,haircut,real_finance_weight,principal1)
      k=input$subjective_k
      mu1<-mu1
      loss1<-loss1
      pro_dict1<-pro_dict_all
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
      df1<-matrix(0,nrow=length(w_1)+1,ncol=3)
      for(i in 1:nrow(df1)){
        if(i==nrow(df1)){
        df1[nrow(df1),1]=mu_adjust%*%w_1+tcost
        df1[nrow(df1),2]=sqrt(w_1%*%cov_adjust%*%w_1)
        df1[nrow(df1),3]=(mu_adjust%*%w_1+tcost-RISK_FREE_RATE)/sqrt(w_1%*%cov_adjust%*%w_1)
        }else{
        w_temp<-rep(0,length(w_1))
        w_temp[i]<-1
        df1[i,1]=mu_adjust%*%w_temp
        df1[i,2]=sqrt(w_temp%*%cov_adjust%*%w_temp)
        df1[i,3]=(mu_adjust%*%w_temp-RISK_FREE_RATE)/sqrt(w_temp%*%cov_adjust%*%w_temp)
        }
      }
      df1<-data.frame(df1)
      return(df1)
    }
    output$tablle<-renderTable({
      weights3<<-weights2()
       ww_1<-weights3$weights[-length(weights3$weights)]
       mat=cost2(w_now,ww_1,trans_cost,finan_cost,haircut,real_finance_weight,principal1)
       colnames(mat)<-c("trans_cost(bps)","finance_cost(bps)","finance_weights","haircut","allocated_weights","initial_Weights","weights_changed")
       optimized_total<-c(apply(mat[,c(1:4)],2,function(x) sum(x)),calculate_total_weights(ww_1,finan_cost,haircut,real_finance_weight),calculate_total_weights(w_now,finan_cost,haircut,real_finance_weight),calculate_total_weights(ww_1,finan_cost,haircut,real_finance_weight)-calculate_total_weights(w_now,finan_cost,haircut,real_finance_weight))
       mat<-rbind(mat,optimized_total)
       df1<-calculate_stat(weights3)
       colnames(df1)<-c('mean_return','volatility','sharpe_ratio')
       mat<-cbind(df1,mat)
       mat<-data.frame(mat)
       name<-as.character(weights3$asset_name)
       name[length(name)]='optimized_overal'
       mat<-cbind(name,mat)
       colnames(mat)[1]="name"
       mat[,2:3]<-as.data.frame(lapply(mat[,2:3], format_num9))
       mat[,4]<-sprintf("%.2f",mat[,4])
       mat[,5:6]<-as.data.frame(lapply(mat[,5:6], format_num8))
       mat[,7:(ncol(mat))]<-as.data.frame(lapply(mat[,7:(ncol(mat))], format_num9))
       mat[nrow(mat),8]=NA
       mat<<-mat
       return(mat)
     }, readOnly = TRUE)

    get_weights_matrix<-eventReactive(input$go9,{
      if(length(which(input$uti2=='combo'))==0){
                weights_matrix<<-call_scenario_cost_BN_matrix(input$uti2,input$assets_num1,input$lambdas,input$asset_ret1,input$asset_vol1,
                              input$asset_corr1,input$sample_number2,input$extreme_stress_loss,pro_dict,
                              input$principal1,trans_cost,finan_cost,haircut,real_finance_weight,input$w_now,input$lower_bound,input$upper_bound,input$ks,
                              input$inequality_constraints,input$equality_constraints,input$maxeval_global,input$maxeval_local,input$asset_name1,input$moment2,smallprob2) 
                weights_matrix<<-data.frame(weights_matrix)
                  return(weights_matrix)
      }else if(length(which(input$uti2!='combo'))==0){
        weights_matrix5<<-call_scenario_cost_BN_matrix2(input$assets_num1,input$x_extreme2,input$x_downturning2,0,input$x_upturning2,input$prob22,input$asset_ret1,input$asset_vol1,
                                                        input$asset_corr1,input$sample_number1,input$extreme_stress_loss,pro_dict,
                                                        input$principal1,trans_cost,finan_cost,haircut,real_finance_weight,input$w_now,input$lower_bound,input$upper_bound,input$ks,
                                                        input$inequality_constraints,input$equality_constraints,input$maxeval_global,input$maxeval_local,input$asset_name1,input$moment2,smallprob2)
        weights_matrix5<<-data.frame(weights_matrix5)
        return(weights_matrix5)
      }else{
        uti2<-input$uti2[-which(input$uti2=='combo')]
        weights_matrix<<-call_scenario_cost_BN_matrix(uti2,input$assets_num1,input$lambdas,input$asset_ret1,input$asset_vol1,
                                                      input$asset_corr1,input$sample_number2,input$extreme_stress_loss,pro_dict,
                                                      input$principal1,trans_cost,finan_cost,haircut,real_finance_weight,input$w_now,input$lower_bound,input$upper_bound,input$ks,
                                                      input$inequality_constraints,input$equality_constraints,input$maxeval_global,input$maxeval_local,input$asset_name1,input$moment2,smallprob2) 
        weights_matrix5<<-call_scenario_cost_BN_matrix2(input$assets_num1,input$x_extreme2,input$x_downturning2,0,input$x_upturning2,input$prob22,input$asset_ret1,input$asset_vol1,
                                                        input$asset_corr1,input$sample_number1,input$extreme_stress_loss,pro_dict,
                                                        input$principal1,trans_cost,finan_cost,haircut,real_finance_weight,input$w_now,input$lower_bound,input$upper_bound,input$ks,
                                                        input$inequality_constraints,input$equality_constraints,input$maxeval_global,input$maxeval_local,input$asset_name1,input$moment2,smallprob2)
      
        }

      weights_matrix<<-data.frame(weights_matrix)
      weights_matrix5<<-data.frame(weights_matrix5)
      return(weights_matrix)
    })

    output$hotable3 <- renderTable({
      weights_matrix2<-get_weights_matrix()
      weights_matrix2[,4:ncol(weights_matrix2)]<-as.data.frame(lapply(weights_matrix2[,4:ncol(weights_matrix2)], format_num9))
    return(weights_matrix2)
      }, readOnly = TRUE)
    output$downloadData <- downloadHandler(
      filename = function() {
        paste('output', input$filetype, sep = ".")
      },
      content <- function(file) {
        sep <- switch(input$filetype, "csv" = ",", "tsv" = "\t")
        write.table(weights_matrix, file, sep = sep,
                    row.names = FALSE)
      }
    )
    output$hotable4 <- renderTable({
      weights_matrix2<-get_weights_matrix()
      if(length(which(input$uti2=='combo'))!=0){
        weights_matrix2<-weights_matrix5
      }
      weights_matrix2[,7:ncol(weights_matrix2)]<-as.data.frame(lapply(weights_matrix2[,7:ncol(weights_matrix2)], format_num9))
      return(weights_matrix2)
    }, readOnly = TRUE)
    Outvar<-reactive({
      vars<-strsplit(input$asset_name1,",")[[1]]
      vars<-c(vars,"optimized_total","user_input")
      vars<-as.list(vars)
      return(vars)
    })
    output$variables=renderUI({
      selectInput('see_return','Asset Name',Outvar())
    })
    output$variables2=renderUI({
      checkboxGroupInput('see_return2','Multiple assets comparison',Outvar(),selected = 'optimized_total')
    })

 out_plot2<-eventReactive(input$go_plot2,{ 
   ###get weights 
   weights3<-weights2()
   w_1<-weights3$weights[-length(weights3$weights)] 
   ###get name of each assets
   name<-as.character(weights3$asset_name[-length(weights3$weights)])
   #####get transaction cost and finance cost for all of these assets
   tcost<-cost(w_now,w_1,trans_cost,finan_cost,haircut,real_finance_weight,principal1)
   ####user input assets initialize
   if(exists("w_user_input")&&length(w_user_input)==nrow(weights3)-1&&warningmessage==0){
     w_user_input<-w_user_input
     tcost2<-tcost2
   }else{
     w_user_input<-w_1
     tcost2<-tcost
   }
   #######calculate nromal scnario return(considering cost) for optimized weights
   rand3<-rand2%*%w_1+tcost
   #######calculate nromal scnario return(considering cost) for user input weights
   rand32<-rand2%*%w_user_input+tcost2
   #########################a big matrix  contain the normal scenario simulation paths
   rand3<-cbind(rand3,rand32,rand2)
   rand3<-data.frame(rand3)
   colnames(rand3)<-c("optimized_total","user_input",name)
   k=input$subjective_k
   ##################calculate utility value for these two assts' weights
   if(input$uti=='expo'){
     u_uti<--expo_find_w(weights3$weights[-length(weights3$weights)],w_now,input$lambda1,trans_cost,finan_cost,haircut,real_finance_weight,principal1,rand2,rand_sub,rand2_moment,loss1,pro_dict,input$subjective_k,mu)
   }else if(input$uti=='power'){
     u_uti<--power_find_w(weights3$weights[-length(weights3$weights)],w_now,input$lambda1,trans_cost,finan_cost,haircut,real_finance_weight,principal1,rand2,rand_sub,rand2_moment,loss1,pro_dict,input$subjective_k,mu)
   }else if(input$uti=='log'){
     u_uti<--log_find_w(weights3$weights[-length(weights3$weights)],w_now,input$lambda1,trans_cost,finan_cost,haircut,real_finance_weight,principal1,rand2,rand_sub,rand2_moment,loss1,pro_dict,input$subjective_k,mu)
   }else{
     u_uti<--combo_find_w(weights3$weights[-length(weights3$weights)],w_now,x_downturning1,AA1,ll1,k_21,k_11,x_11,trans_cost,finan_cost,haircut,real_finance_weight,principal1,rand2,rand_sub,rand2_moment,loss1,pro_dict,input$subjective_k,mu)
   }
   if(input$uti=='expo'){
     u_uti2<-(1-k)*mean(na.omit(get_expo_ut(rand2 %*% w_user_input,input$lambda1,w_now,w_user_input,tcost2)))+k/(1-pro_dict$p0[1])*get_extreme_expo_uti(pro_dict,loss1,input$lambda1,w_now,w_user_input,tcost2,mu,rand_sub)
   }else if(input$uti=='power'){
     u_uti2<-(1-k)*mean(na.omit(get_power_ut(rand2 %*% w_user_input,input$lambda1,w_now,w_user_input,tcost2)))+(k/(1-pro_dict$p0[1]))*get_extreme_power_uti(pro_dict,loss1,input$lambda1,w_now,w_user_input,tcost2,mu,rand_sub)
   }else if(input$uti=='log'){
     u_uti2<-(1-k)*mean(na.omit(get_log_ut(rand2 %*% w_user_input,input$lambda1,w_now,w_user_input,tcost2)))+k/(1-pro_dict$p0[1])*get_extreme_log_uti(pro_dict,loss1,input$lambda1,w_now,w_user_input,tcost2,mu,rand_sub)
   }else{
     u_uti2<-(1-k)*mean(na.omit(get_combo_ut(rand2 %*% w_user_input,x_downturning,AA,ll1,k_2,k_1,x_1,w_now,w_user_input,tcost2)))+(k/(1-pro_dict$p0[1]))*get_extreme_combo_uti(pro_dict,loss1,x_downturning,AA,ll1,k_2,k_1,x_1,w_now,w_user_input,tcost2,mu,rand_sub)
   }
   ########################initialize M M records which assets we want to see its CDF
   M<-0
   for(i in 1:length(input$see_return2)){
     M=c(M,which(colnames(rand3)==input$see_return2[i]))
   }
   M<-M[-1]
   #########################################################if the user input weights is infeasible and we want to see the CDF of this weights,dont display it
   if(warningmessage!=0&&length(which(M==2))!=0){
     M<-M[-which(M==2)]
   }
   # pro_dict2<-as.matrix(pro_dict2)
   ######################pro_dict is bayesian network matrix
   pro_dict2_all<-pro_dict_all
   pro_dict2<-pro_dict
   u<-rep(1,nrow(pro_dict2_all)-1)
   u2<-rep(1,nrow(pro_dict2_all)-1)
   uu<-diag(length(w_1))
   ################pro_dict3 records the extreme evnt loss for assts use optimized weights
   pro_dict3<-t(t(pro_dict2_all)[2:ncol(pro_dict2_all),2:nrow(pro_dict2_all)]%*%(loss1*w_1))+tcost+loss1[length(loss1)]*w_1[length(w_1)]
   ################pro_dict3 records the extreme evnt loss for assts use user input weights   
   pro_dict32<-t(t(pro_dict2_all)[2:ncol(pro_dict2_all),2:nrow(pro_dict2_all)]%*%(loss1*w_user_input))+tcost2+loss1[length(loss1)]*w_user_input[length(w_user_input)]
   ################pro_dict4 records extreme event for each assts
   pro_dict4<-matrix(1,nrow=(ncol(pro_dict2_all)-1),ncol=(nrow(pro_dict2_all)-1))-t(pro_dict2_all)[2:ncol(pro_dict2_all),2:nrow(pro_dict2_all)]
   pro_dict4[,length(w_1)]<-rep(0,nrow(pro_dict4))
   pro_dict_loss_matrix<-data.frame(t(pro_dict2_all)[2:ncol(pro_dict2_all),2:nrow(pro_dict2_all)]%*%(loss1*uu))
#    if(is.null(rownames(pro_dict_loss_matrix))){
#      rownames(pro_dict_loss_matrix)=colnames(pro_dict2_all)[2]
#    }
   rand5<-c()
   ##############################################################################################################################
   #####start to iterate from M
   #############devide 2 parts:k=1 and k!=1
   ############for k=1,we don't need normal scenarios and k!=1 we need to combine normal scenario and extreme scenario
   ##########after these 2 parts we also need devide assets into 4 groups:optimized weights,user input weights,each assets except cash, cash
   ##########for the first two weights,we also need to add back other assets return when assets don't happen extreme events
   ##############################################################################################################################
   for(jk in 1:length(M)){
     rand6=rand3[,M[jk]]
     if(input$subjective_k!=1){  ##k!=1
       if(M[jk]==1){ #optimized total
         if(is.null(colnames(pro_dict3))){
          Numberr<-as.double(str_replace_all(colnames(unique(pro_dict2)),"p",""))[2]
         }else{
           Numberr<-as.double(str_replace_all(colnames(unique(pro_dict3)),"p",""))  ##record scenario's order         
         }
         ###record each scenario's number
         Number<-as.integer((input$subjective_k)/((1-input$subjective_k)*(1-pro_dict2$p0[1]))*nrow(rand3)*pro_dict2_all[1,1+Numberr])
         Numberr<-Numberr[which(Number!=0)]
         Number<-Number[which(Number!=0)]
         if(length(Number)==0){
           rand6<-data.frame(rand6)
         }else{
         ecdf_gene<-rep(pro_dict3[Numberr],Number)
         ecdf_gene2<-c()
         ######for each scenario generate return for assets don't happen extreme events
         for(i in 1:length(Number)){
           if(Number[i]%%2==1){
             temp_matrix<-generate_N_rand(length(w_1),asset_corr,margins_r,paramMargins_r,Number[i]+1)[1:Number[i],] 
             if(Number[i]==1){
               ecdf_gene2<-c(ecdf_gene2,sum((pro_dict4[Numberr[i],]*w_1)*t(temp_matrix))) 
             }else{
               ecdf_gene2<-c(ecdf_gene2,(pro_dict4[Numberr[i],]*w_1)%*%t(temp_matrix))
             }
           }else if(Number[i]!=0){
             temp_matrix<- generate_N_rand(length(w_1),asset_corr,margins_r,paramMargins_r,Number[i])
             ecdf_gene2<-c(ecdf_gene2,(pro_dict4[Numberr[i],]*w_1)%*%t(temp_matrix))
           }else{
             next
           }
         }
         ###add together,this is the total return for the extreme parts
         ecdf_gene<-ecdf_gene2+ecdf_gene
         ###combine normal parts and extreme parts
         rand6<-data.frame(c(rand6,ecdf_gene))
         }
         ####add a name for this
         colnames(rand6)<-colnames(rand3)[M[jk]]
       }else if(M[jk]==2){   ##user input weights same idea
         if(is.null(colnames(pro_dict32))){
           Numberr<-as.double(str_replace_all(colnames(unique(pro_dict2)),"p",""))[2]
         }else{
           Numberr<-as.double(str_replace_all(colnames(unique(pro_dict32)),"p",""))  ##record scenario's order         
         }
         Number<-as.integer((input$subjective_k)/((1-input$subjective_k)*(1-pro_dict2$p0[1]))*nrow(rand3)*pro_dict2_all[1,1+Numberr])
         Numberr<-Numberr[which(Number!=0)]
         Number<-Number[which(Number!=0)]
         if(length(Number)==0){
           rand6<-data.frame(rand6)
         }else{
           ecdf_gene<-rep(pro_dict32[Numberr],Number)
           ecdf_gene2<-c()
           for(i in 1:length(Number)){
             if(Number[i]%%2==1){
               temp_matrix<-generate_N_rand(length(w_user_input),asset_corr,margins_r,paramMargins_r,Number[i]+1)[1:Number[i],] 
               if(Number[i]==1){
                 ecdf_gene2<-c(ecdf_gene2,sum((pro_dict4[Numberr[i],]*w_user_input)*t(temp_matrix))) 
               }else{
                 ecdf_gene2<-c(ecdf_gene2,(pro_dict4[Numberr[i],]*w_user_input)%*%t(temp_matrix))
               }
             }else if(Number[i]!=0){
               temp_matrix<- generate_N_rand(length(w_user_input),asset_corr,margins_r,paramMargins_r,Number[i])
               ecdf_gene2<-c(ecdf_gene2,(pro_dict4[Numberr[i],]*w_user_input)%*%t(temp_matrix))
             }else{
               next
             }
           }
           ecdf_gene<-ecdf_gene2+ecdf_gene
           rand6<-data.frame(c(rand6,ecdf_gene))
         }
         colnames(rand6)<-colnames(rand3)[M[jk]]
       }else if(M[jk]!=(length(name)+2)){  ###each single assets
         set<-unique(pro_dict_loss_matrix[,M[jk]-2])
         number<-set[which(set!=0)]
         if(length(which((pro_dict_loss_matrix[,M[jk]-2])!=0))==0){
           Number<-0
         }else{
           Number<-as.integer((input$subjective_k)/((1-input$subjective_k)*(1-pro_dict2$p0[1]))*nrow(rand3)*sum(pro_dict2_all[1,1+as.double(str_replace_all(rownames(pro_dict_loss_matrix[which((pro_dict_loss_matrix[,M[jk]-2])!=0),]),"p",""))]))
         }
          if(length(which((pro_dict_loss_matrix[,M[jk]-2])==0))==0){
           Number1<-0
         }else{
         Number1<-as.integer((input$subjective_k)/((1-input$subjective_k)*(1-pro_dict2$p0[1]))*nrow(rand3)*sum(pro_dict2_all[1,1+as.double(str_replace_all(rownames(pro_dict_loss_matrix[which((pro_dict_loss_matrix[,M[jk]-2])==0),]),"p",""))]))
         }
         Number<-as.integer(nrow(rand3)/(nrow(rand3)+Number1)*Number)
         ecdf_gene<-rep(number,Number)
         rand6<-data.frame(c(rand6,ecdf_gene))
         colnames(rand6)<-colnames(rand3)[M[jk]]}else{  ####cash           
           rand6<-data.frame(rand6)
           colnames(rand6)<-c('cash')}
       if(jk==1){
         rand5<-as.vector(as.matrix(rand6))
         name_level_list<-c(colnames(rand6))
         length_level_list<-c(length(rand5))
       }else{
         rand5<-c(as.vector(rand5),as.vector(as.matrix(rand6)))
         name_level_list<-c(name_level_list,colnames(rand6))
         length_level_list<-c(length_level_list,length(as.vector(as.matrix(rand6))))
       }}else{  ###when k=1
         if(M[jk]==1){
           if(is.null(colnames(pro_dict3))){
             Numberr<-as.double(str_replace_all(colnames(unique(pro_dict2)),"p",""))[2]
           }else{
             Numberr<-as.double(str_replace_all(colnames(unique(pro_dict3)),"p",""))  ##record scenario's order         
           }
           ###record each scenario's number
            Number<-as.integer(nrow(rand3)*pro_dict2_all[1,1+Numberr])
            Numberr<-Numberr[which(Number!=0)]
            Number<-Number[which(Number!=0)]
            ecdf_gene<-rep(pro_dict3[Numberr],Number)
            Number_plot<<-Number
            ecdf_gene2<-c()
            ######for each scenario generate return for assets don't happen extreme events
            for(i in 1:length(Number)){
              if(Number[i]%%2==1){
                temp_matrix<-generate_N_rand(length(w_1),asset_corr,margins_r,paramMargins_r,Number[i]+1)[1:Number[i],] 
                if(Number[i]==1){
                  ecdf_gene2<-c(ecdf_gene2,sum((pro_dict4[Numberr[i],]*w_1)*t(temp_matrix))) 
                }else{
                  ecdf_gene2<-c(ecdf_gene2,(pro_dict4[Numberr[i],]*w_1)%*%t(temp_matrix))
                }
              }else if(Number[i]!=0){
                temp_matrix<- generate_N_rand(length(w_1),asset_corr,margins_r,paramMargins_r,Number[i])
                ecdf_gene2<-c(ecdf_gene2,(pro_dict4[Numberr[i],]*w_1)%*%t(temp_matrix))
              }else{
                next
              }
            }
            ecdf_gene<-ecdf_gene+ecdf_gene2
            ecdf_gene_plot<<-ecdf_gene
           ###don't need normal simulation just extreme
           rand6<-data.frame(ecdf_gene)
           colnames(rand6)<-colnames(rand3)[M[jk]]
         }else if(M[jk]==2){
           if(is.null(colnames(pro_dict32))){
             Numberr<-as.double(str_replace_all(colnames(unique(pro_dict2)),"p",""))[2]
           }else{
             Numberr<-as.double(str_replace_all(colnames(unique(pro_dict32)),"p",""))  ##record scenario's order         
           }
#            Number<-as.integer((input$subjective_k)/((1-input$subjective_k)*(1-pro_dict2$p0[1]))*nrow(rand3)*pro_dict2_all[1,1+Numberr])
           Number<-as.integer(nrow(rand3)*pro_dict2_all[1,1+Numberr])
           Numberr<-Numberr[which(Number!=0)]
           Number<-Number[which(Number!=0)]
           ecdf_gene<-rep(pro_dict32[Numberr],Number)
           ecdf_gene2<-c()
           for(i in 1:length(Number)){
             if(Number[i]%%2==1){
               temp_matrix<-generate_N_rand(length(w_user_input),asset_corr,margins_r,paramMargins_r,Number[i]+1)[1:Number[i],] 
               if(Number[i]==1){
                 ecdf_gene2<-c(ecdf_gene2,sum((pro_dict4[Numberr[i],]*w_user_input)*t(temp_matrix))) 
               }else{
                 ecdf_gene2<-c(ecdf_gene2,(pro_dict4[Numberr[i],]*w_user_input)%*%t(temp_matrix))
               }
             }else if(Number[i]!=0){
               temp_matrix<- generate_N_rand(length(w_user_input),asset_corr,margins_r,paramMargins_r,Number[i])
               ecdf_gene2<-c(ecdf_gene2,(pro_dict4[Numberr[i],]*w_user_input)%*%t(temp_matrix))
             }else{
               next
             }
           }
           ecdf_gene<-ecdf_gene2+ecdf_gene
           rand6<-data.frame(ecdf_gene)
           colnames(rand6)<-colnames(rand3)[M[jk]]
         }else if(M[jk]!=(length(name)+2)){
           set<-unique(pro_dict_loss_matrix[,M[jk]-2])
           number<-set[which(set!=0)]
           if(length(which((pro_dict_loss_matrix[,M[jk]-2])!=0))==0){
             Number<-0
           }else{
             Number<-as.integer(nrow(rand3)*sum(pro_dict2_all[1,1+as.double(str_replace_all(rownames(pro_dict_loss_matrix[which((pro_dict_loss_matrix[,M[jk]-2])!=0),]),"p",""))]))
           }
           if(length(which((pro_dict_loss_matrix[,M[jk]-2])==0))==0){
             Number2<-0
           }else{
             Number2<-as.integer(nrow(rand3)*sum(pro_dict2_all[1,1+as.double(str_replace_all(rownames(pro_dict_loss_matrix[which((pro_dict_loss_matrix[,M[jk]-2])==0),]),"p",""))]))
           }
           Number<-as.integer(nrow(rand3)/Number2*Number)
           ecdf_gene<-rep(number,Number)
           rand6<-data.frame(c(rand6,ecdf_gene))
           colnames(rand6)<-colnames(rand3)[M[jk]]
         }else{
           rand6<-data.frame(rand6)
           colnames(rand6)<-c('cash')}
         if(jk==1){
           rand5<-as.vector(as.matrix(rand6))
           name_level_list<-c(colnames(rand6))
           length_level_list<-c(length(rand5))
         }else{
           rand5<-c(as.vector(rand5),as.vector(as.matrix(rand6)))
           name_level_list<-c(name_level_list,colnames(rand6))
           length_level_list<-c(length_level_list,length(as.vector(as.matrix(rand6))))
         }
       }}
   gl=as.factor(rep(name_level_list,length_level_list))
   #####combine weights and assets'  name
   df<-data.frame(x=rand5,assets=gl)
   if(length(M)==1&M==(length(name)+1)){
     p2<-ggplot(df, aes(x,colour=assets)) +stat_ecdf(geom="step")+xlab(colnames(rand6))+ylab("Return Cumulative Density Function")+scale_y_continuous(breaks = round(seq(0,1, by = 0.05),5))+geom_vline(xintercept = 0,colour=1)+scale_x_continuous(labels=scales::percent,breaks = seq(from=round(100*input$x_range_return2[1])/100,to=input$x_range_return2[2], by = 0.05))+coord_cartesian(xlim=input$x_range_return2)
     for(i in seq(from=input$x_range_return2[1],to=input$x_range_return2[2], by = 0.01)){
       p2=p2+geom_vline(xintercept = i,colour=1,linetype="dotted")
     }
     for(i in seq(from=input$x_range_return2[1],to=input$x_range_return2[2], by = 0.05)){
       p2=p2+geom_vline(xintercept = i,colour=1)
     }
   }else{
     p2<-ggplot(df, aes(x,colour=assets)) +stat_ecdf(geom="step")+xlab(colnames(rand6))+ylab("Return Cumulative Density Function")+scale_y_continuous(breaks = round(seq(0,1, by = 0.05),5))+scale_x_continuous(labels=scales::percent,breaks = seq(from=round(100*input$x_range_return2[1])/100,to=input$x_range_return2[2], by = 0.05))+geom_vline(xintercept = 0,colour=1)+coord_cartesian(xlim=input$x_range_return2)
     for(i in seq(from=input$x_range_return2[1],to=input$x_range_return2[2], by = 0.01)){
       p2=p2+geom_vline(xintercept = i,colour=1,linetype="dotted")
     }
     for(i in seq(from=input$x_range_return2[1],to=input$x_range_return2[2], by = 0.05)){
       p2=p2+geom_vline(xintercept = i,colour=1)
     }
   }
   if(length(which(M==1))==0)
   {
     if(length(which(M==2))==0){
       return(p2)
     }else{
      p2<-p2+annotate(geom = "text", x = input$x_range_return2[1], y = 0.95,
                   label = paste0("user input utility value:",sprintf("%.4f",u_uti2)),
                   hjust = 0,colour="red",size=7) 
      return(p2)
     }}else{
       if(length(which(M==2))==0){
         p2<-p2+annotate(geom = "text", x = input$x_range_return2[1], y = 0.95,
                         label = paste0("optimized utility value:",sprintf("%.4f",u_uti)),
                         hjust = 0,colour="red",size=7)
         return(p2)
       }else{
         p2<-p2+annotate(geom = "text", x = input$x_range_return2[1], y = 0.95,
                         label = paste(paste0("optimized utility value:",sprintf("%.4f",u_uti),"\n"),paste0("user input utility value:",sprintf("%.4f",u_uti2))),
                         hjust = 0,colour="red",size=7)
         return(p2)
       }
     }
 })
 output$returnplot2<-renderPlot({
   out_plot2()
   })
 warnstring<<-reactiveValues(data=NULL)
 ret_owntable<<-function(own_check,subjective_k){
   if(recalculate==0){
     if(exists('mat')){
       mat<<-data.frame(lapply(mat, as.character), stringsAsFactors=FALSE)
       ownvalues<<-reactiveValues(data=as.data.frame(mat))
       mat2<<-mat
     }else{
       ownvalues<<-NULL
       return();
     }
   }
   if(recalculate!=0) mat2<<-as.data.frame(ownvalues$data)
   ######if check(which means automatically calculate finance weights we need to calculate optimized finance strategy)
   #########if fail just return -1
   if(own_check==TRUE){
     if(!exists('mat2'))  return();
     w_1<-as.double(str_replace_all(mat2[,9],"%",""))/100
     w_1<-w_1[-length(w_1)]
     w<-assign_finance_weights(w_1,finan_cost,haircut,real_finance_weight,principal1)
     if(w[1]<0) warningmessage<<--2
     w<-c(w,sum(w))
     w<-sprintf("%.2f",w*100)
     mat2[,7]<-paste0(w,"%")
   } 
   recalculate<<-1
   #####start to compute
   mat2[nrow(mat2),9]<-paste0(calculate_total_weights(as.double(str_replace_all(mat2[,9],"%",""))[-nrow(mat2)]/100,finan_cost,haircut,real_finance_weight)*100,"%")
   mat2[nrow(mat2),7]<-paste0(sum(as.double(str_replace_all(mat2[,7],"%",""))[-nrow(mat2)]),"%")
   mat2[,11]<-as.double(str_replace_all(mat2[,9],"%",""))-as.double(str_replace_all(mat2[,10],"%",""))
   mat2[,11]<-paste0(sprintf("%.2f",mat2[,11]),"%")
   w_now<-as.double(str_replace_all(mat2[,10],"%",""))/100
   w_now<-w_now[-length(w_now)]
   w_1<-as.double(str_replace_all(mat2[,9],"%",""))/100
   w_1<-w_1[-length(w_1)]
   # w_1<-weights3$weights[-length(weights3$weights)] 
   mat2[,5]<-c(trans_cost_vector(w_now,w_1,trans_cost1,principal1)*10000,sum(trans_cost_vector(w_now,w_1,trans_cost1,principal1)*10000))
   mat2[,5]<-sprintf("%.2f",mat2[,5])
   w_finan<-as.double(str_replace_all(mat2[,7],"%",""))/100
   w_finan<-w_finan[-length(w_finan)]
   mat2[,6]<-c(finance_cost_vector(w_finan,finan_cost1,principal1)*10000,sum(finance_cost_vector(w_finan,finan_cost1,principal1)*10000))
   mat2[,6]<-sprintf("%.2f",mat2[,6])
   rett<-get_mean_vol_sharpe(pro_dict_all,w_1,(as.double(mat2[nrow(mat2),5])+as.double(mat2[nrow(mat2),6]))/10000,subjective_k,mu1,loss1,asset_cov)
   rett<-c(paste0(sprintf("%.2f",rett[1:2]*100),"%"),sprintf("%.2f",rett[3]))
   mat2[nrow(mat2),2:4]<-rett
   mat2[1:(nrow(mat2)-1),2]=paste0(as.double(str_replace_all(mat[1:(nrow(mat)-1),2],"%","")),"%")
   mat2[1:(nrow(mat2)-1),3]=paste0(as.double(str_replace_all(mat[1:(nrow(mat)-1),3],"%","")),"%")
   mat2[1:(nrow(mat2)-1),4]=paste0(as.double(str_replace_all(mat[1:(nrow(mat)-1),4],"%","")),"%")
   warningmessage<<-getwarningmessage(w_1,w_now,w_finan,haircut,real_finance_weight,finan_cost1,lower_bound,upper_bound)
   warnstring$data<<-getwarnstring();
   ownvalues<<-reactiveValues(data=as.data.frame(mat2))
   w_user_input<<-w_1
   tcost2<<-(as.double(mat2[nrow(mat2),5])+as.double(mat2[nrow(mat2),6]))/10000
   return(ownvalues$data)
 }
#  observe({
#    # generate_owntable$data<<-ret_owntable(FALSE,input$subjective_k)
#    generate_owntable$data<<-ret_owntable(input$own_check,input$subjective_k)
#  })
# observeEvent(input$go_own,{
#   generate_owntable$data<<-ret_owntable(input$own_check,input$subjective_k)
# })

observeEvent(input$go_start_compute,{
  generate_owntable$data<<-ret_owntable(input$own_check,input$subjective_k)
  recalculate<<-0
})
generate_owntable<-reactiveValues(data = NULL)
observe({
if(!is.null(input$owntable))
ownvalues$data <- hot_to_r(input$owntable)
})

 output$owntable<- renderRHandsontable({
   generate_owntable$data<<-ret_owntable(input$own_check,input$subjective_k)
   if(is.null(generate_owntable$data)) return()
   rhandsontable(generate_owntable$data,readOnly =TRUE) %>%
     hot_col(col="finance_weights",readOnly = FALSE ) %>%
     hot_col(col="allocated_weights",readOnly = FALSE )
 })
 #######warnmessage output for user input weights
 getwarnstring<-function(){
   if(warningmessage>1){
     s<- paste0("WARNING:check haircut condtion for your input! code=",warningmessage)
   }else if(warningmessage==1){
     s<- paste0("WARNING:cash position for short margin doesn't meet! code=",warningmessage)
   }else if(warningmessage==0){
     s<-"ALL CHECK!"
   }else if(warningmessage==-1){
     s<-paste0("WARNING:TOTAL AMOUNT Of finance Weight is not right! code=",warningmessage)
   }else if(warningmessage<(-1)&warningmessage>(-2)){
     s<-paste0("WARNING:Your input breaks the lower bound! The first asset number is ",round(-10000*(warningmessage+1)))
   }else if(warningmessage==-2){
     s<-paste0("WARNING:Your input weights is infeasible for getting a possible solution to do finance! code is",warningmessage)
   }else if(warningmessage<(-2)&warningmessage>(-3)){ 
     s<-paste0("WARNING:Your input breaks the upper bound! The first asset number is ",round(-10000*(warningmessage+2)))
   }else if(warningmessage<(-3)&warningmessage>(-4)){
     if(round(-10000*(warningmessage+3))==1){
       s<-paste0("WARNING:Your input weights don't meet inequality constraints! Break constraints: cash constraints")
     }else if(round(-10000*(warningmessage+3))<length(w_now)+2){
       s<-paste0("WARNING:Your input weights don't meet inequality constraints! Break constraints: haircut constraints,asset number",round(-10000*(warningmessage+3))-1)
     }else{
       s<-paste0("WARNING:Your input weights don't meet inequality constraints! Break constraints: inequality constraints NO.",round(-10000*(warningmessage+3))-length(w_now)+1)
     }
   }else if(warningmessage==-4){
     s<-paste0("WARNING:Your input weights don't meet equality constraints! code=",warningmessage)
   }
   return(s)
 }

#  observeEvent(input$go_own,{
#    warnstring$data<<-getwarnstring();
#  })
 observeEvent(input$go_start_compute,{
   warnstring$data<<-getwarnstring();
 })
 output$warnmsg<-renderText({
   warnstring$data
 })
 globals<-eventReactive(input$go_start_compute,{
   if(globalsearch==0){
     s<-"globalsearch complete!"
   }else if(globalsearch==-1){
     s<-"Use Moment matching method to skip global search"
   }else{
     s<-"globalsearch failed,Only local optimal available at this time!"
   }
 })
 
 output$gsmsg<-renderText({
   globals()
 })
 #  out_plot<-eventReactive(input$go_plot,{
 #    weights3<-weights2() 
 #    w_1<-weights3$weights[-length(weights3$weights)] 
 #    name<-as.character(weights3$asset_name[-length(weights3$weights)])
 #    tcost<-cost(w_now,w_1,trans_cost,finan_cost,haircut,real_finance_weight,principal1)
 #    if(exists("w_user_input")&&length(w_user_input)==nrow(weights3)-1){
 #      w_user_input<-w_user_input
 #      tcost2<-tcost2
 #    }else{
 #      w_user_input<-w_1
 #      tcost2<-tcost
 #    }
 #    if(input$uti=='expo'){
 #      u_uti<--expo_find_w(weights3$weights[-length(weights3$weights)],w_now,input$lambda1,trans_cost,finan_cost,haircut,real_finance_weight,principal1,rand2,loss1,pro_dict,k,mu)
 #    }else if(input$uti=='power'){
 #      u_uti<--power_find_w(weights3$weights[-length(weights3$weights)],w_now,input$lambda1,trans_cost,finan_cost,haircut,real_finance_weight,principal1,rand2,loss1,pro_dict,k,mu)
 #    }else if(input$uti=='log'){
 #      u_uti<--log_find_w(weights3$weights[-length(weights3$weights)],w_now,input$lambda1,trans_cost,finan_cost,haircut,real_finance_weight,principal1,rand2,loss1,pro_dict,k,mu)
 #    }else{
 #      u_uti<--combo_find_w(weights3$weights[-length(weights3$weights)],w_now,x_downturning1,AA1,ll1,k_21,k_11,x_11,trans_cost,finan_cost,haircut,real_finance_weight,principal1,rand2,loss1,pro_dict,k,mu)
 #    }
 #    if(input$uti=='expo'){
 #      u_uti2<-(1-k)*mean(na.omit(get_expo_ut(rand2 %*% w_user_input,input$lambda1,w_now,w_user_input,tcost2)))+k/(1-pro_dict$p0[1])*get_extreme_expo_uti(pro_dict,loss1,input$lambda1,w_now,w_user_input,tcost2,mu)
 #    }else if(input$uti=='power'){
 #      u_uti2<-(1-k)*mean(na.omit(get_power_ut(rand2 %*% w_user_input,input$lambda1,w_now,w_user_input,tcost2)))+(k/(1-pro_dict$p0[1]))*get_extreme_power_uti(pro_dict,loss1,input$lambda1,w_now,w_user_input,tcost2,mu)
 #    }else if(input$uti=='log'){
 #      u_uti2<-(1-k)*mean(na.omit(get_log_ut(rand2 %*% w_user_input,input$lambda1,w_now,w_user_input,tcost2)))+k/(1-pro_dict$p0[1])*get_extreme_log_uti(pro_dict,loss1,input$lambda1,w_now,w_user_input,tcost2,mu)
 #    }else{
 #      u_uti2<-(1-k)*mean(na.omit(get_combo_ut(rand2 %*% w_user_input,x_downturning,AA,ll,k_2,k_1,x_1,w_now,w_user_input,tcost2)))+(k/(1-pro_dict$p0[1]))*get_extreme_combo_uti(pro_dict,loss1,x_downturning,AA,ll,k_2,k_1,x_1,w_now,w_user_input,tcost2,mu)
 #    }
 #    rand3<-rand2%*%w_1+tcost
 #    rand32<-rand2%*%w_user_input+tcost2
 #    rand3<-cbind(rand3,rand32,rand2)
 #    rand3<-data.frame(rand3)
 #    colnames(rand3)<-c("optimized_total","user_input",name)
 #    rand5<-c()
 #    M=which(input$see_return==colnames(rand3))
 #      u<-rep(1,nrow(pro_dict2)-1)
 #      u2<-rep(1,nrow(pro_dict2)-1)
 #      uu<-diag(length(w_1))
 #      pro_dict3<-t(t(pro_dict2)[2:ncol(pro_dict2),2:nrow(pro_dict2)]%*%(loss1*w_1))+tcost+loss1[length(loss1)]*w_1[length(w_1)]
 #      pro_dict32<-t(t(pro_dict2)[2:ncol(pro_dict2),2:nrow(pro_dict2)]%*%(loss1*w_user_input))+tcost2+loss1[length(loss1)]*w_user_input[length(w_user_input)]
 #      pro_dict4<-matrix(1,nrow=(ncol(pro_dict2)-1),ncol=(nrow(pro_dict2)-1))-t(pro_dict2)[2:ncol(pro_dict2),2:nrow(pro_dict2)]
 #      pro_dict4[,length(w_1)]<-rep(0,nrow(pro_dict4))
 #      pro_dict_loss_matrix<-data.frame(t(pro_dict2)[2:ncol(pro_dict2),2:nrow(pro_dict2)]%*%(loss1*uu))
 #      rand4=rand3[,M]
 #      if(M==1){
 #        if(input$subjective_k!=1){
 #       Numberr<-as.double(str_replace_all(colnames(unique(pro_dict3)),"p",""))
 #       Number<-as.integer((input$subjective_k)/((1-input$subjective_k)*(1-pro_dict2$p0[1]))*nrow(rand3)*pro_dict2[1,1+as.double(str_replace_all(colnames(unique(pro_dict3)),"p",""))])
 #       ecdf_gene<-rep(unique(pro_dict3),Number)
 #       ecdf_gene2<-c()
 #       for(i in 1:length(Number)){
 #         if(Number[i]%%2==1){
 #           temp_matrix<-generate_N_rand(length(w_1),asset_corr,margins_r,paramMargins_r,Number[i]+1)[1:Number[i],] 
 #           ecdf_gene2<-c(ecdf_gene2,pro_dict4[Numberr[i],]%*%t(temp_matrix))
 #         }else if(Number[i]!=0){
 #           temp_matrix<- generate_N_rand(length(w_1),asset_corr,margins_r,paramMargins_r,Number[i])
 #           ecdf_gene2<-c(ecdf_gene2,pro_dict4[Numberr[i],]%*%t(temp_matrix))
 #         }else{
 #           next
 #         }
 #       }
 #       ecdf_gene<-ecdf_gene2+ecdf_gene
 #       rand4<-data.frame(c(rand4,ecdf_gene))
 #       colnames(rand4)<-colnames(rand3)[M]
 #       p2<-ggplot(rand4, aes(x=rand4)) +
 #         stat_ecdf(geom="step")+xlab(colnames(rand4))+ylab("Return Cumulative Density Function")+scale_y_continuous(breaks = round(seq(0,1, by = 0.05),5))+scale_x_continuous(labels=scales::percent,breaks =  seq(from=round(100*input$x_range_return[1])/100,to=input$x_range_return[2], by = 0.05))+geom_vline(xintercept = 0,colour=1)+coord_cartesian(xlim=input$x_range_return)
 #       for(i in seq(from=input$x_range_return[1],to=input$x_range_return[2], by = 0.01)){
 #         p2=p2+geom_vline(xintercept = i,colour=1,linetype="dotted")
 #       }
 #       for(i in seq(from=input$x_range_return[1],to=input$x_range_return[2], by = 0.05)){
 #         p2=p2+geom_vline(xintercept = i,colour=1)
 #       }
 #       p2<-p2+annotate(geom = "text", x = input$x_range_return[1], y = 0.95,
 #                 label = paste0("optimized utility value:",sprintf("%.4f",u_uti)),
 #                 hjust = 0,colour="red",size=7)
 #       return(p2)}else{
 #         Numberr<-as.double(str_replace_all(colnames(unique(pro_dict3)),"p",""))
 #         Number<-as.integer(nrow(rand3)*pro_dict2[1,1+as.double(str_replace_all(colnames(unique(pro_dict3)),"p",""))])
 #         ecdf_gene<-rep(unique(pro_dict3),Number)
 #         ecdf_gene2<-c()
 #         for(i in 1:length(Number)){
 #           if(Number[i]%%2==1){
 #             temp_matrix<-generate_N_rand(length(w_1),asset_corr,margins_r,paramMargins_r,Number[i]+1)[1:Number[i],] 
 #             ecdf_gene2<-c(ecdf_gene2,pro_dict4[Numberr[i],]%*%t(temp_matrix))
 #           }else if(Number[i]!=0){
 #             temp_matrix<- generate_N_rand(length(w_1),asset_corr,margins_r,paramMargins_r,Number[i])
 #             ecdf_gene2<-c(ecdf_gene2,pro_dict4[Numberr[i],]%*%t(temp_matrix))
 #           }else{
 #             next
 #           }
 #         }
 #         ecdf_gene<-ecdf_gene2+ecdf_gene
 #         rand4<-data.frame(ecdf_gene)
 #         colnames(rand4)<-colnames(rand3)[M]
 #         p2<-ggplot(rand4, aes(x=rand4)) +
 #           stat_ecdf(geom="step")+xlab(colnames(rand4))+ylab("Return Cumulative Density Function")+scale_y_continuous(breaks = round(seq(0,1, by = 0.05),5))+scale_x_continuous(labels=scales::percent,breaks =  seq(from=round(100*input$x_range_return[1])/100,to=input$x_range_return[2], by = 0.05))+geom_vline(xintercept = 0,colour=1)+coord_cartesian(xlim=input$x_range_return)
 #         for(i in seq(from=input$x_range_return[1],to=input$x_range_return[2], by = 0.01)){
 #           p2=p2+geom_vline(xintercept = i,colour=1,linetype="dotted")
 #         }
 #         for(i in seq(from=input$x_range_return[1],to=input$x_range_return[2], by = 0.05)){
 #           p2=p2+geom_vline(xintercept = i,colour=1)
 #         }
 #         p2<-p2+annotate(geom = "text", x = input$x_range_return[1], y = 0.95,
 #                         label = paste0("optimized utility value:",sprintf("%.4f",u_uti)),
 #                         hjust = 0,colour="red",size=7)
 #         return(p2)
 #       }}else if(M==2){
 #         if(input$subjective_k!=1){
 #           Numberr<-as.double(str_replace_all(colnames(unique(pro_dict32)),"p",""))
 #           Number<-as.integer((input$subjective_k)/((1-input$subjective_k)*(1-pro_dict2$p0[1]))*nrow(rand3)*pro_dict2[1,1+as.double(str_replace_all(colnames(unique(pro_dict32)),"p",""))])
 #           ecdf_gene<-rep(unique(pro_dict32),Number)
 #           ecdf_gene2<-c()
 #           for(i in 1:length(Number)){
 #             if(Number[i]%%2==1){
 #               temp_matrix<-generate_N_rand(length(w_user_input),asset_corr,margins_r,paramMargins_r,Number[i]+1)[1:Number[i],] 
 #               ecdf_gene2<-c(ecdf_gene2,pro_dict4[Numberr[i],]%*%t(temp_matrix))
 #             }else if(Number[i]!=0){
 #               temp_matrix<- generate_N_rand(length(w_user_input),asset_corr,margins_r,paramMargins_r,Number[i])
 #               ecdf_gene2<-c(ecdf_gene2,pro_dict4[Numberr[i],]%*%t(temp_matrix))
 #             }else{
 #               next
 #             }
 #           }
 #           ecdf_gene<-ecdf_gene2+ecdf_gene
 #           rand4<-data.frame(c(rand4,ecdf_gene))
 #           colnames(rand4)<-colnames(rand3)[M]
 #           p2<-ggplot(rand4, aes(x=rand4)) +
 #             stat_ecdf(geom="step")+xlab(colnames(rand4))+ylab("Return Cumulative Density Function")+scale_y_continuous(breaks = round(seq(0,1, by = 0.05),5))+scale_x_continuous(labels=scales::percent,breaks =  seq(from=round(100*input$x_range_return[1])/100,to=input$x_range_return[2], by = 0.05))+geom_vline(xintercept = 0,colour=1)+coord_cartesian(xlim=input$x_range_return)
 #           for(i in seq(from=input$x_range_return[1],to=input$x_range_return[2], by = 0.01)){
 #             p2=p2+geom_vline(xintercept = i,colour=1,linetype="dotted")
 #           }
 #           for(i in seq(from=input$x_range_return[1],to=input$x_range_return[2], by = 0.05)){
 #             p2=p2+geom_vline(xintercept = i,colour=1)
 #           }
 #           p2<-p2+annotate(geom = "text", x = input$x_range_return[1], y = 0.95,
 #                           label = paste0("User input utility value:",sprintf("%.4f",u_uti2)),
 #                           hjust = 0,colour="red",size=7)
 #           return(p2)}else{
 #             Numberr<-as.double(str_replace_all(colnames(unique(pro_dict32)),"p",""))
 #             Number<-as.integer(nrow(rand3)*pro_dict2[1,1+as.double(str_replace_all(colnames(unique(pro_dict32)),"p",""))])
 #             ecdf_gene<-rep(unique(pro_dict32),Number)
 #             ecdf_gene2<-c()
 #             for(i in 1:length(Number)){
 #               if(Number[i]%%2==1){
 #                 temp_matrix<-generate_N_rand(length(w_user_input),asset_corr,margins_r,paramMargins_r,Number[i]+1)[1:Number[i],] 
 #                 ecdf_gene2<-c(ecdf_gene2,pro_dict4[Numberr[i],]%*%t(temp_matrix))
 #               }else if(Number[i]!=0){
 #                 temp_matrix<- generate_N_rand(length(w_user_input),asset_corr,margins_r,paramMargins_r,Number[i])
 #                 ecdf_gene2<-c(ecdf_gene2,pro_dict4[Numberr[i],]%*%t(temp_matrix))
 #               }else{
 #                 next
 #               }
 #             }
 #             ecdf_gene<-ecdf_gene2+ecdf_gene
 #             rand4<-data.frame(ecdf_gene)
 #             colnames(rand4)<-colnames(rand3)[M]
 #             p2<-ggplot(rand4, aes(x=rand4)) +
 #               stat_ecdf(geom="step")+xlab(colnames(rand4))+ylab("Return Cumulative Density Function")+scale_y_continuous(breaks = round(seq(0,1, by = 0.05),5))+scale_x_continuous(labels=scales::percent,breaks =  seq(from=round(100*input$x_range_return[1])/100,to=input$x_range_return[2], by = 0.05))+geom_vline(xintercept = 0,colour=1)+coord_cartesian(xlim=input$x_range_return)
 #             for(i in seq(from=input$x_range_return[1],to=input$x_range_return[2], by = 0.01)){
 #               p2=p2+geom_vline(xintercept = i,colour=1,linetype="dotted")
 #             }
 #             for(i in seq(from=input$x_range_return[1],to=input$x_range_return[2], by = 0.05)){
 #               p2=p2+geom_vline(xintercept = i,colour=1)
 #             }
 #             p2<-p2+annotate(geom = "text", x = input$x_range_return[1], y = 0.95,
 #                             label = paste0("User input utility value:",sprintf("%.4f",u_uti2)),
 #                             hjust = 0,colour="red",size=7)
 #             return(p2)
 #           }
 #       }else if(M!=(length(name)+2)){
 #        if(input$subjective_k!=1){
 #        set<-unique(pro_dict_loss_matrix[,M-1])
 #        number<-set[which(set!=0)]
 #        Number<-as.integer((input$subjective_k)/((1-input$subjective_k)*(1-pro_dict2$p0[1]))*nrow(rand3)*sum(pro_dict2[1,1+as.double(str_replace_all(rownames(pro_dict_loss_matrix[which((pro_dict_loss_matrix[,M-1])!=0),]),"p",""))]))
 #        Number1<-as.integer((input$subjective_k)/((1-input$subjective_k)*(1-pro_dict2$p0[1]))*nrow(rand3)*sum(pro_dict2[1,1+as.double(str_replace_all(rownames(pro_dict_loss_matrix[which((pro_dict_loss_matrix[,M-1])==0),]),"p",""))]))
 #        Number<-as.integer(nrow(rand3)/(nrow(rand3)+Number1)*Number)
 #        ecdf_gene<-rep(number,Number)
 #       rand4<-data.frame(c(rand4,ecdf_gene))
 #       colnames(rand4)<-colnames(rand3)[M]
 #        p2<-ggplot(rand4, aes(x=rand4)) +stat_ecdf(geom="step")+xlab(colnames(rand4))+ylab("Return Cumulative Density Function")+scale_y_continuous(breaks = round(seq(0,1, by = 0.05),5))+scale_x_continuous(labels=scales::percent,breaks =  seq(from=round(100*input$x_range_return[1])/100,to=input$x_range_return[2], by = 0.05))+geom_vline(xintercept = 0,colour=1)+coord_cartesian(xlim=input$x_range_return)
 #        for(i in seq(from=input$x_range_return[1],to=input$x_range_return[2], by = 0.01)){
 #          p2=p2+geom_vline(xintercept = i,colour=1,linetype="dotted")
 #        }
 #        for(i in seq(from=input$x_range_return[1],to=input$x_range_return[2], by = 0.05)){
 #          p2=p2+geom_vline(xintercept = i,colour=1)
 #        }
 #        return(p2)}else{
 #          set<-unique(pro_dict_loss_matrix[,M-1])
 #          number<-set[which(set!=0)]
 #          Number<-as.integer(nrow(rand3)*sum(pro_dict2[1,1+as.double(str_replace_all(rownames(pro_dict_loss_matrix[which((pro_dict_loss_matrix[,M-1])!=0),]),"p",""))]))
 #          Number2<-as.integer(nrow(rand3)*sum(pro_dict2[1,1+as.double(str_replace_all(rownames(pro_dict_loss_matrix[which((pro_dict_loss_matrix[,M-1])==0),]),"p",""))]))
 #          Number<-as.integer(nrow(rand3)/Number2*Number)
 #          ecdf_gene<-rep(number,Number)
 #          rand4<-data.frame(c(rand4,ecdf_gene))
 #          colnames(rand4)<-colnames(rand3)[M]
 #          p2<-ggplot(rand4, aes(x=rand4)) +stat_ecdf(geom="step")+xlab(colnames(rand4))+ylab("Return Cumulative Density Function")+scale_y_continuous(breaks = round(seq(0,1, by = 0.05),5))+scale_x_continuous(labels=scales::percent,breaks =  seq(from=round(100*input$x_range_return[1])/100,to=input$x_range_return[2], by = 0.05))+geom_vline(xintercept = 0,colour=1)+coord_cartesian(xlim=input$x_range_return)
 #          for(i in seq(from=input$x_range_return[1],to=input$x_range_return[2], by = 0.01)){
 #            p2=p2+geom_vline(xintercept = i,colour=1,linetype="dotted")
 #          }
 #          for(i in seq(from=input$x_range_return[1],to=input$x_range_return[2], by = 0.05)){
 #            p2=p2+geom_vline(xintercept = i,colour=1)
 #          }
 #          return(p2)}
 #      }else{
 #        rand4<-rand3[,M]
 #         p2<-ggplot(rand3, aes(x=rand4)) +
 #          stat_ecdf(geom="step")+xlab(colnames(rand3)[M])+ylab("Return Cumulative Density Function")+scale_y_continuous(breaks = round(seq(0,1, by = 0.05),5))+geom_vline(xintercept = 0,colour=1)+scale_x_continuous(labels=scales::percent,breaks = seq(from=round(100*input$x_range_return[1])/100,to=input$x_range_return[2], by = 0.05))+coord_cartesian(xlim=input$x_range_return)
 #         for(i in seq(from=input$x_range_return[1],to=input$x_range_return[2], by = 0.01)){
 #           p2=p2+geom_vline(xintercept = i,colour=1,linetype="dotted")
 #         }
 #         for(i in seq(from=input$x_range_return[1],to=input$x_range_return[2], by = 0.05)){
 #           p2=p2+geom_vline(xintercept = i,colour=1)
 #         }
 #        return(p2)
 #      }
 #    })
 #  
 #  output$returnplot<-renderPlot({
 #   out_plot()
 #  })

  }
)


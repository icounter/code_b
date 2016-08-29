#####include packages
library(shiny)
library(shinysky)
library(shinyAce)
library(rhandsontable)
###warning close
options(warn=-1)
modes <- getAceModes()
themes <- getAceThemes()
###########################################################################################################################
############################################global variable initialze######################################################
###########################################################################################################################
MAXIMUM_LOSS<<--99999999999999      ##this is used to stand for a very negative numble
MINIMUM_RET<<-0.01                  ##this is used to make log utility and exponential utility function smooth
trans_generaed<<-0                  ##a flag record whether we have generated trans&finance cost or not 
NN<<-7                              ##asset numble initialize  
RISK_FREE_RATE<<-0.0065             ##risk free rate
haircut<<-c(0,0,0,1)                ##haircut vector initialize
real_finance_weight<<-c(1,1,1,1)    ##real finance weights vector initialize
globalsearch<<-0                    ##warning message records global search succeed or not 
warningmessage<<-0                  ##warning message records the user input varibales feasible or not
tolerance<<-0.0002                  ##1 bps tolerance for user input
smallprob<<-1*1e-5;smallprob2<<-1*1e-5  ##negelect scenarios happen with less than this probability
# sum(pro_dict[1,which(pro_dict[1,]>smallprob)])
# length(which(pro_dict[1,]>smallprob))
sub_sample<<-1000                    ##how much to simulate 
###########################################################################################################################
############################################optimization control###########################################################
###########################################################################################################################
ftol_abs<<-1e-8;ftol_abs_global<<-1e-2   ###control objective functions value change_abs
xtol_abs<<-1e-4;xtol_abs_global<<-2e-1    ###control every varibales values' smallest value change 
########here we allow a very slack change for global optimizer and comparatively rigid change for local optimizer
print_level<<-3    ########control print level this controls the output information of our optimizer 
## first 2 tol vector controls the global optimization if the value is within the range then we think the solution meets the constraints     
tol_constraints_eq<<-1e-2;tol_constraints_ineq<<-1e-2;recalculate<<-0 
####max evaluation rouns this controls the speed of search process
maxeval<<-300
displayList<-c("expo","power","log","combo")
##############################################################################################################################
########################################APP UI starts here####################################################################
##############################################################################################################################
shinyUI(fluidPage(
  img(src="logo.jpg", height = 50, width = 200),hr(),  ##load image
  sidebarLayout(
      sidebarPanel(
                  fileInput('file1', 'Choose CSV File',accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
                  helpText("Upload your input csv file here and it will automatically update your assumptions"),
                  helpText("You can also update your assumtions by hand"),
                  tags$hr(),
                  checkboxInput('header', 'Header', TRUE),
                  radioButtons('sep', 'Separator',c(Comma=',',Semicolon=';',Tab='\t'),','),
                  radioButtons('quote', 'Quote',c(None='','Double Quote'='"','Single Quote'="'"),selected = '"')),  ##load scv file),
       mainPanel(tableOutput('contents'))
               ),hr(),
  titlePanel("Mean-Variance Framework and Exponential Utility Function"),    
  helpText("Note:This is only a samll example to show the consistency of exponential utility function and mean-variance framework"),
  helpText("Set the assumption step by step,number of samples means how many simulation paths you would use when using exponential utility"),
  fluidRow(
          column(1,numericInput("assets_num",label = h3("Assets Number"), value = 14)),
          column(2,numericInput("lambda", label = h3("Lambda"), value = 5)),
          column(2,textInput("asset_ret", label = h3("Return of each asset"), value = "0.0635,0.0747,0.1835,0.097,0.0492,0.0769,0.0749,0.1035,0.199,0.2338,0.0619,0.0888,0.1221,0.0065")),
          column(2,textInput("asset_vol", label = h3("vol of each asset"), value="0.068,0.146,0.189,0.004,0.068,0.146,0.06,0.094,0.201,0.274,0.021,0.055,0.202,0")),
          column(2,textInput("asset_corr", label = h3("correlation of each asset"), value = "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0")),
          column(2, textInput("lower_bound1",label = h3("Lower bounds for assets"), value = "0,0,0,0,0,0,0,0,0,0,0,0,0,0")),
          column(2,textInput("upper_bound1", label = h3("Upper bounds for assets"), value = "1,1,1,1,1,1,1,1,1,1,1,1,1,1")),
          column(2, sliderInput("sample_number", label = h3("number of samples"),min = 10000, max = 1000000,step=10000, value =50000)),
          column(2,actionButton("go1","start to compute"))
          ),
  hr(),
  fluidRow(plotOutput("ggplot")),### plot based on results
  ###first part ends here
  hr(),
  titlePanel("Real scenario considering transaction cost and set up Bayesian Net"),
  hr(),
  titlePanel("1.Setup normal return part of assets"),
  hr(),
  fluidRow(
           column(1,numericInput("assets_num1",label = h3("Assets Number"),value = 4)),
           helpText("Note:All the below inputs are comma delimited"),
           column(2,textInput("asset_name1", label = h3("Names of each Asset"), value = "CLO,senior_loan,NPL1stpay,cash")),
           column(2,textInput("asset_ret1", label = h3("Return of each asset"), value = "0.05,0.035,0.04,0.0045")),
           column(2,textInput("asset_vol1", label = h3("Volatility of each asset"), value="0.085,0.022,0.025,0")),
           column(2,textInput("asset_corr1", label = h3("correlation of each asset"), value = "0.6,0.6,0,0.6,0,0"))),
  hr(),
  titlePanel("2.Setup Extreme part of assets"),
           helpText("Note:Define it before you want you want to do assets optimization"),
  sidebarLayout(
  sidebarPanel(
               textInput("extreme_stress_loss",label = h3("Extreme Loss of each Asset"), value ="-.2,-0.03,-0.03,0.0065"),
               helpText("Note:This is the return of each assests when extreme event happens"),
               textInput("affect_relation", label = h3("Affect Relationship between assets"), value ="Hidden1;CLO+Hidden1;senior_loan+Hidden1;NPL1stpay+Hidden1;cash"),
               helpText("Define Structure of your BN,semicomma delimited.and A+B means A is affected by B"),
               actionButton("go_plot_bayesian_net","Generate Bayesian Net plot"),
               helpText("outplot the structure of BN at the right panel"),
               textInput("prob_table", label = h3("Probability table"), value ="10,90;80,20,10,90;80,20,10,90;80,20,10,90;0,100"),
               helpText("Define conditional probability table for each node,semicomma delimited")
               ),
  mainPanel(
    plotOutput("plot_bayesian"),
    helpText("In the plot,A->B means A will affect  B"),
    helpText("Joint probability table shows the probability value in each scenario,1 means  extreme events will happen and 0 means extreme events won't happen")
           )
               ),
  hr(),
  titlePanel("Joint probability table"),
  sidebarLayout(
  sidebarPanel(
  helpText("Condtional probability matrix, Line i and row j mean if i happens the probability of j happens"),
  sliderInput("tentative_k", label = h3("Rescale the probabbility of Bayesian Network"),min = 0.0, max =1.0,step=0.01, value =0.01),actionButton("go_BN_cond","Show BN conditional Table"),
  selectInput("show_part_or_all", label = h3("show correlation matrix"),choices = list("part" = 'part1', "all" = 'all1'), selected = 'part1')),
  mainPanel(DT::dataTableOutput("hotable2"),hr(),titlePanel("Normal Scenario Correlation Matrix"),
                 DT::dataTableOutput("bayesian_matrix_cor_normal"),  ##normal correaltion matrix
                 hr(),
                 conditionalPanel("input.show_part_or_all=='all1'",             ##conditional panel display when we select to show all
                 titlePanel("Extreme Scenario Correlation Matrix"),
                 DT::dataTableOutput("bayesian_matrix_cor_extreme"),
                 hr()),
                 titlePanel("Conditional probability Matrix in Extreme Scenario"),
                 DT::dataTableOutput("bayesian_matrix_cor"),
                 hr(),
                 conditionalPanel("input.show_part_or_all=='all1'",
                                 titlePanel("Downside correlation in Extreme Scenario"),
                                 DT::dataTableOutput("bayesian_matrix_cor_downside"),hr()),
                 titlePanel("Rescaled Correlation Matrix in Bayesian Network"),
                 DT::dataTableOutput("bayesian_matrix_cor2"))
          ),
  hr(),
  titlePanel("3.Set Transaction Cost & Finance cost"),
  sidebarLayout(
  sidebarPanel(
    numericInput("principal1", label = h3("Principal($M)"), value = 1720),
    helpText("How much money we hold right now"),
    textInput("w_now", label = h3("Present Portofolio Position"), value = "0,0,0,0"),
    helpText("present portofolio position,comma delimited"),
    textInput("w1", label = h3("Tentative Portofolio Position"), value = "0,0.25,0.25,0"),
    helpText("The portofolio position you would like to change to,click compute transaction cost to see how much you should pay for position change,comma delimited"),
    actionButton("go_generate_trans","generate transcation cost Matrix"),
    actionButton("go_generate_finan","generate finance cost Matrix"),
    helpText("You can make your piecewise linear assumption here,input the number of your intervals")
    ),
  mainPanel(
    tabPanel("trans_cost",h4("Transaction Cost matrix"),div(class="well container-fluid",DT::dataTableOutput("hotable1"))),  
    tabPanel("finan_cost" ,h4("Finance Cost matrix"),div(class="well container-fluid",DT::dataTableOutput("hotable5"))))),
  hr(),
  titlePanel("How much you will spend on transaction cost and finance cost"), 
  actionButton("go_tentative_cost","start to compute transaction cost and finance cost"),
  tableOutput("tentative_transcost"),
  hr(),
  titlePanel("4.Select Utility Function"),
  sidebarLayout(
    sidebarPanel(
    selectInput("uti", label = h3("Utility Function"),choices = list("Exponential Utility" = 'expo', "Power Utility" = 'power', "Log Utility" = 'log',"Risk seeking combination"='combo'), selected = 'expo'),   #
    helpText("Four kinds of utility functions to choose from:power,log,expo and comb."),
    conditionalPanel("input.uti!='combo'",
    numericInput("lambda1",label = h3("Input Lambda"),  value = 5),
    helpText("This is a two-scenario selection question:(1) you can get the mid return with 100%;(2) with a probability of p to get downside lose and (1-p) to win upside return."),
    helpText("Make these two scenarios indifference to you")),
    conditionalPanel("input.uti=='combo'",
    numericInput("x_extreme",label = h3("extreme loss tolerance"),step=0.01,value = -0.3),
    numericInput("x_downturning",label = h3("downside turning return"),step=0.01,value = -0.05),
    numericInput("x_upturning",label = h3("upside turning return"),step=0.01,value = 0.05),
    numericInput("prob2",label = h3("prob to get downside return"),step=0.01,value = 0.8),
    helpText("This assumes when a PM has a small negative return or a small positive return, he will be risk seeking to get more positive return."),
    helpText("downside turning return means before he loses this much he is risk seeking, after that he is risk averse"),
    helpText("upside turning return means before he wins this much he is risk seeking, after that he is risk averse"),
    helpText("Mid return and probability to fail is the same idea as other three functions,to make yourself indifference between these two scenarios"))
    ),  
  mainPanel(
    conditionalPanel("input.uti=='combo'",hr(),sliderInput("x_range2", label = h3("Show utility in the range"), min = -1, max = 1, step=0.01,value = c(-0.3, 0.2)),
                     plotOutput("utility_out2")),
    conditionalPanel("input.uti!='combo'",titlePanel("Reference plot for lambda:"), 
                     helpText("By changing the range of negative return you can change the display range of right plot"),
                     helpText("Left plot is the plot for (1) zero return,(2) equal return and loss in the other scenario and the probability you set to lose when makes this an indifference choice"),
                     helpText("Right plot is the plot for (1) zero return,(2) equal probability of win and lose,if knowing how much you will lose in the bad scenario how much you will need to win to compensate for it"),
                     sliderInput("neg_ret_ran", label = h3("select the range of negative return"), min = -0.2, max = 0, step=0.01,value = -0.11),
                     splitLayout(cellWidths = c("50%", "50%"), plotOutput("plot_ut_zero_equal"),plotOutput("plot_ut_zero_equal_2")),
                     hr(),
    sliderInput("x_range", label = h3("Show utility in the range"), min = -1, max = 1, step=0.01,value = c(-0.3, 0.2)),
      plotOutput("utility_out"))
    )),
    hr(),
  titlePanel("5.Set bounds & subjective value"),
  fluidRow(
    helpText("The lower bound of your weights,can't short then make lower bound to be a zero vector"),
    column(2, textInput("lower_bound",label = h3("Lower bounds for assets"),value = "0,0,0,0")),
    helpText("The upper bound of your weights,if you can't do leverage then make upper bound to be a vector of one"),
    column(2,textInput("upper_bound", label = h3("Upper bounds for assets"), value = "1,1,1,1")),
    column(2, sliderInput("sample_number1", label = h3("number of samples"),min = 10000, max = 1000000,step=10000, value =10000)),
    helpText("Subjective Value parameter is the probability you assigned that in the next interval at least one extreme event will happen"),
    column(2, sliderInput("subjective_k", label = h3("Subjective Value"),min = 0.0, max =1.0,step=0.01, value =0.01)),
    column(2, sliderInput("maxeval_global", label = h3("Global or Moment Match Search rounds"),min = 50, max =1000,step=50, value =200)),
    column(2, sliderInput("maxeval_local", label = h3("Final Step Search rounds"),min = 50, max =1000,step=50, value =200)),
    column(2, checkboxInput("moment", label = "Use Moment matching method to replace global search part", value = TRUE))
    # column(2, sliderInput("smallprob", label = h3("negelect scenarios with small probabilities"),min = 0, max = 1e-4,step=1e-6, value =1e-5))
    ),hr(),
  fluidRow(
##################################################################################################################################
#################################input equality and inequality constraints########################################################
##################################################################################################################################
  column(6, aceEditor("inequality_constraints", value="##input inequality bounds here(hint format is g(x)>=0)
eval_g0 <- function(w1,w_now=w_now1,beta1=beta11,rand2=rand21,rand_sub=rand_sub1,rand2_moment=rand21_moment,loss1=loss11,pro_dict=pro_dict1,k=k1,trans_cost=trans_cost1,finan_cost=finan_cost1,haircut=haircut2,real_finance_weight=real_finance_weight1,principal1=principal11,mu=mu1)
{
list1=which(real_finance_weight!=1)
haircut1=(1-haircut)*real_finance_weight*w1
h=ifelse(length(list1)!=0,w1[length(w1)]-sum(abs(((1-haircut2)*real_finance_weight*w1)[list1]))-sum(abs((haircut*w1)[list1])),1)
h1=ifelse(length(list1)!=0,sum(w1[-list1])-1,sum(w1)-1) ###borrow part should be equal to 1
h<--c(h,cost_vector(w1,finan_cost,haircut,real_finance_weight,principal1),h1)
return(h)
}")),
  column(6,aceEditor("equality_constraints", value="##input equality bounds here
#     eval_h0 <- function(w1,w_now=w_now1,beta1=beta11,rand2=rand21,rand_sub=rand_sub1,rand2_moment=rand21_moment,loss1=loss11,pro_dict=pro_dict1,k=k1,trans_cost=trans_cost1,finan_cost=finan_cost1,haircut=haircut2,real_finance_weight=real_finance_weight1,principal1=principal11,mu=mu1)
#     {
#      return(sum(w1)-1)
#     }")),
  column(2,actionButton("go_start_compute","start to compute"))),
  hr(),
  fluidRow(plotOutput("opt_weights"),
  textOutput('gsmsg'),tags$head(tags$style("#gsmsg{color: red;font-size: 20px;font-style: italic;}")),
  hr(),
  tableOutput("tablle"),hr(),
  #actionButton("go_own","Submit your input weights"),
  checkboxInput("own_check", label = "Automatically calculate finance weights for given weights", value = TRUE), rHandsontableOutput('owntable'),
        textOutput('warnmsg'),tags$head(tags$style("#warnmsg{color: red;font-size: 20px;font-style: italic;}"))),
     
  hr(),
######################################################################
#################################single CDF part commented############
######################################################################
# 
# titlePanel("Cumulative Density Function based on Optimized Weights"),
#  sidebarLayout(
#  sidebarPanel(uiOutput('variables'),
#               sliderInput("x_range_return", label = h3("Show Return CDF in the range"), min = -1, max = 1, step=0.05,value = c(-0.3, 0.3)),
#               actionButton("go_plot","start to generate plot")
#              ),mainPanel(plotOutput("returnplot"))
#  ),
# hr(),
######################################################################
#################################multi CDF part ######################
######################################################################
sidebarLayout(
  sidebarPanel(uiOutput('variables2'),sliderInput("x_range_return2", label = h3("Show Return CDF in the range"), min = -1, max = 1, step=0.05,value = c(-0.3, 0.3)),
              actionButton("go_plot2","start to generate plot")),
  mainPanel(plotOutput("returnplot2",height = 500))),
hr(),
######################################################################
#################################PART 6 ##############################
######################################################################
titlePanel("6.Optimize weights table"),
sidebarLayout(
  sidebarPanel(textInput("lambdas", label = h3("lambdas you are interested in"), value = "5,10"),
  textInput("ks", label = h3("probabilities of at least one extreme events to happen"), value = "0.2,0.5"),
  helpText("Above are the lambdas and subjective_ks you are interested in,you may input multiple lambdas and ks,comma delimited"),
  checkboxInput("moment2", label = "Use Moment matching method to replace global search part", value = TRUE),
  # sliderInput("smallprob2", label = h3("negelect scenarios with small probabilities"),min = 0, max = 1e-4, step=1e-6,value =1e-5),
  checkboxGroupInput("uti2", label = h3("Utility function typies you are interested in"), 
                               choices = list("Exponential Utility" = 'expo', "Power Utility" = 'power', "Log Utility" = 'log',"Risk seeking combination"='combo'), selected = 'expo'),
  conditionalPanel(condition="input.uti2.indexOf('combo')>-1",hr(),
                   helpText("Below is just for COMBO function"),
                   numericInput("x_extreme2",label = h3("extreme loss tolerance"),step=0.01,value = -0.3),
                   numericInput("x_downturning2",label = h3("downside turning return"),step=0.01,value = -0.05),
                   # numericInput("x_mid22",label = h3("mid return between downside and upside return"),step=0.01,value = 0),
                   numericInput("x_upturning2",label = h3("upside turning return"),step=0.01,value = 0.05),
                   numericInput("prob22",label = h3("prob to get downside return"),step=0.01,value = 0.8)
                   ),hr(),
  sliderInput("sample_number2", label = h3("number of samples to generate"),min = 10000, max = 1000000,step=10000, value =50000),
  actionButton("go9","click here to show tables"),
  radioButtons("filetype", "File type:",choices = c("csv", "tsv")),
  downloadButton('downloadData','Download'),
  helpText("This download function is only available when open in browser")
  ),mainPanel(
    titlePanel("The optimized weights of the assets are:"),
    conditionalPanel(condition="input.uti2.indexOf('expo')>-1||input.uti2.indexOf('log')>-1||input.uti2.indexOf('power')>-1",tableOutput("hotable3")),
    conditionalPanel(condition="input.uti2.indexOf('combo')>-1",hr(),titlePanel("Using Combination Utility function:"),tableOutput("hotable4"))
    )
)))
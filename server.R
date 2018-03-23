library(shiny)
library(mc2d)
library(ggplot2)
library(scales)
library(hexbin)
library(plyr)
source("multiplot.R")
source("ggGuage.R") 


shinyServer(function(input,output){
  
  simulate <- eventReactive(input$Analyse, {
  #simulate <- reactive( {
    N <- input$noOfSimulations

    PLEFestimate <- data.frame(L = input$PLEFMin,  ML = input$PLEFML,  H = input$PLEFMax, CONF = as.numeric(input$PLEFConf))
    TEFestimate <- data.frame(L = input$TEFMin,  ML = input$TEFML,  H = input$TEFMax, CONF = as.numeric(input$TEFConf))
    CFestimate <- data.frame(L = input$CFMin,  ML = input$CFML,  H = input$CFMax, CONF = as.numeric(input$CFConf))
    POAestimate <- data.frame(L = input$POAMin,  ML = input$POAML,  H = input$POAMax, CONF = as.numeric(input$POAConf))
    Vulnestimate <- data.frame(L = input$VulnMin,  ML = input$VulnML,  H = input$VulnMax, CONF = as.numeric(input$VulnConf))
    TCestimate <- data.frame(L = input$TCMin,  ML = input$TCML,  H = input$TCMax, CONF = as.numeric(input$TCConf))
    RSestimate <- data.frame(L = input$RSMin,  ML = input$RSML,  H = input$RSMax, CONF = as.numeric(input$RSConf))
    SLEFestimate <- data.frame(L = input$SLEFMin,  ML = input$SLEFML,  H = input$SLEFMax, CONF = as.numeric(input$SLEFConf))
    PrimaryProductivityestimate <- data.frame(L = input$PrimaryProductivityMin,  ML = input$PrimaryProductivityML,  H = input$PrimaryProductivityMax, CONF = as.numeric(input$PrimaryProductivityConf))
    PrimaryResponseestimate <- data.frame(L = input$PrimaryResponseMin,  ML = input$PrimaryResponseML,  H = input$PrimaryResponseMax, CONF = as.numeric(input$PrimaryResponseConf))
    PrimaryReplacementestimate <- data.frame(L = input$PrimaryReplacementMin,  ML = input$PrimaryReplacementML,  H = input$PrimaryReplacementMax, CONF = as.numeric(input$PrimaryReplacementConf))
    PrimaryFandJestimate <- data.frame(L = input$PrimaryFandJMin,  ML = input$PrimaryFandJML,  H = input$PrimaryFandJMax, CONF = as.numeric(input$PrimaryFandJConf))
    PrimaryCompAdvestimate <- data.frame(L = input$PrimaryCompAdvMin,  ML = input$PrimaryCompAdvML,  H = input$PrimaryCompAdvMax, CONF = as.numeric(input$PrimaryCompAdvConf))
    PrimaryReputationestimate <- data.frame(L = input$PrimaryReputationMin,  ML = input$PrimaryReputationML,  H = input$PrimaryReputationMax, CONF = as.numeric(input$PrimaryReputationConf))
    SecondaryProductivityestimate <- data.frame(L = input$SecondaryProductivityMin,  ML = input$SecondaryProductivityML,  H = input$SecondaryProductivityMax, CONF = as.numeric(input$SecondaryProductivityConf))
    SecondaryResponseestimate <- data.frame(L = input$SecondaryResponseMin,  ML = input$SecondaryResponseML,  H = input$SecondaryResponseMax, CONF = as.numeric(input$SecondaryResponseConf))
    SecondaryReplacementestimate <- data.frame(L = input$SecondaryReplacementMin,  ML = input$SecondaryReplacementML,  H = input$SecondaryReplacementMax, CONF = as.numeric(input$SecondaryReplacementConf))
    SecondaryFandJestimate <- data.frame(L = input$SecondaryFandJMin,  ML = input$SecondaryFandJML,  H = input$SecondaryFandJMax, CONF = as.numeric(input$SecondaryFandJConf))
    SecondaryCompAdvestimate <- data.frame(L = input$SecondaryCompAdvMin,  ML = input$SecondaryCompAdvML,  H = input$SecondaryCompAdvMax, CONF = as.numeric(input$SecondaryCompAdvConf))
    SecondaryReputationestimate <- data.frame(L = input$SecondaryReputationMin,  ML = input$SecondaryReputationML,  H = input$SecondaryReputationMax, CONF = as.numeric(input$SecondaryReputationConf))
    ServiceDowntimeestimate <- data.frame(L = input$ServiceDowntimeMin,  ML = input$ServiceDowntimeML,  H = input$ServiceDowntimeMax, CONF = as.numeric(input$ServiceDowntimeConf))
    DataLossestimate <- data.frame(L = input$DataLossMin,  ML = input$DataLossML,  H = input$DataLossMax, CONF = as.numeric(input$DataLossConf))
    DataLossScale <- data.frame(SVL = input$SVLow, HL = input$HLow, HH = input$HHigh, SGL = input$SGLow, SGH = input$SGHigh, ML = input$MLow, MH = input$MHigh, LL = input$LLow, LH = input$LHigh, VLH = input$VLHigh, LP = input$lossPercentile)


    # Change after Code Review - 2015
    # Run rpert calculation under different seed for vulnerability simulation
    if(input$chkRS && input$chkTC){
      set.seed(input$vulnSeed)
      TCsamples <- sapply(rep(NA,N),function(x){rpert(N, TCestimate$L/100, TCestimate$ML/100, TCestimate$H/100, shape = TCestimate$CONF)})
      RSsamples <- sapply(rep(NA,N),function(x){rpert(N, RSestimate$L/100, RSestimate$ML/100, RSestimate$H/100, shape = RSestimate$CONF)})
    }
    # set seed for rest of model
    set.seed(input$seed)
    
    if(input$chkPOA && input$chkCF){
      CFsamples <- rpert(N, CFestimate$L, CFestimate$ML, CFestimate$H, shape = CFestimate$CONF)
      POAsamples <- rpert(N, POAestimate$L, POAestimate$ML, POAestimate$H, shape = POAestimate$CONF)
      TEFsamples <- POAsamples * CFsamples
    } 

    if(input$chkTEF){
      TEFsamples <- rpert(N, TEFestimate$L, TEFestimate$ML, TEFestimate$H, shape = TEFestimate$CONF)
    }

    # Change after Code Review - 2015
    # Note that Vulnerability Monte Carlo simulation has been simplified to a vector calculation.
    # The colSums() function operates down the colum to collapse an NxN matrix to an Nx1 matrix
    if(input$chkRS && input$chkTC){
      vuln <- TCsamples > RSsamples
      Vulnsamples = (colSums(vuln)/N) #* 100
      PLEFsamples <- TEFsamples * Vulnsamples
    }

    if(input$chkVuln){
      Vulnsamples <- rpert(N, Vulnestimate$L/100, Vulnestimate$ML/100, Vulnestimate$H/100, shape = Vulnestimate$CONF)
      PLEFsamples <- TEFsamples * Vulnsamples
    }

    if(input$chkPLEF){
      PLEFsamples <- rpert(N, PLEFestimate$L, PLEFestimate$ML, PLEFestimate$H, shape = PLEFestimate$CONF)
    }

    PrimaryProductivitysamples <- rpert(N, PrimaryProductivityestimate$L, PrimaryProductivityestimate$ML, PrimaryProductivityestimate$H, shape = PrimaryProductivityestimate$CONF)
    PrimaryResponsesamples <- rpert(N, PrimaryResponseestimate$L, PrimaryResponseestimate$ML, PrimaryResponseestimate$H, shape = PrimaryResponseestimate$CONF)
    PrimaryReplacementsamples <- rpert(N, PrimaryReplacementestimate$L, PrimaryReplacementestimate$ML, PrimaryReplacementestimate$H, shape = PrimaryReplacementestimate$CONF)
    PrimaryFandJsamples <- rpert(N, PrimaryFandJestimate$L, PrimaryFandJestimate$ML, PrimaryFandJestimate$H, shape = PrimaryFandJestimate$CONF)
    PrimaryCompAdvsamples <- rpert(N, PrimaryCompAdvestimate$L, PrimaryCompAdvestimate$ML, PrimaryCompAdvestimate$H, shape = PrimaryCompAdvestimate$CONF)
    PrimaryReputationsamples <- rpert(N, PrimaryReputationestimate$L, PrimaryReputationestimate$ML, PrimaryReputationestimate$H, shape = PrimaryReputationestimate$CONF)
    PLMsamples <- PrimaryProductivitysamples + PrimaryResponsesamples + PrimaryReplacementsamples + PrimaryFandJsamples + PrimaryCompAdvsamples + PrimaryReputationsamples
    #PLMsamples <- sum(PrimaryProductivitysamples, PrimaryResponsesamples, PrimaryReplacementsamples, PrimaryFandJsamples, PrimaryCompAdvsamples, PrimaryReputationsamples)

    SecondaryProductivitysamples <- rpert(N, SecondaryProductivityestimate$L, SecondaryProductivityestimate$ML, SecondaryProductivityestimate$H, shape = SecondaryProductivityestimate$CONF)
    SecondaryResponsesamples <- rpert(N, SecondaryResponseestimate$L, SecondaryResponseestimate$ML, SecondaryResponseestimate$H, shape = SecondaryResponseestimate$CONF)
    SecondaryReplacementsamples <- rpert(N, SecondaryReplacementestimate$L, SecondaryReplacementestimate$ML, SecondaryReplacementestimate$H, shape = SecondaryReplacementestimate$CONF)
    SecondaryFandJsamples <- rpert(N, SecondaryFandJestimate$L, SecondaryFandJestimate$ML, SecondaryFandJestimate$H, shape = SecondaryFandJestimate$CONF)
    SecondaryCompAdvsamples <- rpert(N, SecondaryCompAdvestimate$L, SecondaryCompAdvestimate$ML, SecondaryCompAdvestimate$H, shape = SecondaryCompAdvestimate$CONF)
    SecondaryReputationsamples <- rpert(N, SecondaryReputationestimate$L, SecondaryReputationestimate$ML, SecondaryReputationestimate$H, shape = SecondaryReputationestimate$CONF)
    SLMsamples <- SecondaryProductivitysamples + SecondaryResponsesamples + SecondaryReplacementsamples + SecondaryFandJsamples + SecondaryCompAdvsamples + SecondaryReputationsamples
    #SLMsamples <- sum(SecondaryProductivitysamples, SecondaryResponsesamples, SecondaryReplacementsamples, SecondaryFandJsamples, SecondaryCompAdvsamples, SecondaryReputationsamples)

    SLEFsamples <- (rpert(N, SLEFestimate$L, SLEFestimate$ML, SLEFestimate$H, shape = SLEFestimate$CONF) / 100) * PLEFsamples

    PrimaryALE <- PLEFsamples * PLMsamples
    SecondaryALE <- SLEFsamples * SLMsamples
    TotalALE <- PrimaryALE + SecondaryALE

#Added Service Downtime and Data Loss
    ServiceDowntimesamples <- rpert(N, ServiceDowntimeestimate$L, ServiceDowntimeestimate$ML, ServiceDowntimeestimate$H, shape = ServiceDowntimeestimate$CONF)
    DataLosssamples <- rpert(N, DataLossestimate$L, DataLossestimate$ML, DataLossestimate$H, shape = DataLossestimate$CONF)

    ALEServiceDowntime <- PLEFsamples * ServiceDowntimesamples
    ALEDataLoss <- PLEFsamples * DataLosssamples

##############################################################################

#########################################################################################################
    
    output$ProbPlot <- renderPlot({
  
      ALEsamples <- simulate()
      
      ninety <- quantile(ALEsamples, probs=(0.9))
      maximum <- quantile(ALEsamples, probs = (1))
      ten <- quantile(ALEsamples, probs = (0.1))
      average <- mean(ALEsamples)
      minimum <- quantile(ALEsamples, probs = (0))
      LossPercentile <- quantile(ALEsamples, probs=(DataLossScale$LP/100))
      
      h <- hist(ALEsamples, breaks=100, col="steelblue", main="Annualised Risk Exposure", xlab="Loss", axes = FALSE, ylab="Simulation Distribution", panel.first=grid(NULL,NULL,lty=1))
      
      par(new = T)
      
      ec <- ecdf(ALEsamples)
      #plot(x = h$mids, y=ec(h$mids)*max(h$counts), col = rgb(0,0,0,alpha=0), axes=FALSE, xlab=NA, ylab=NA)
      axis(1, at=axTicks(1), labels=sprintf("£%s", axTicks(1)))
      
      
      abline(v = LossPercentile, col ="red", lwd = 3, lty = 2)
      
      abline(v = ninety, col ="black", lwd = 3, lty = 2)
      text(ninety, max(h$counts)-2, "90th", col = "black", pos = 4, srt = 45)
      
      abline(v = maximum, col ="black", lwd = 3, lty = 2)
      text(maximum, max(h$counts)-2, "Max", col = "black", pos = 4, srt = 45)
      
      abline(v = average, col ="black", lwd = 3, lty = 2)
      text(average, max(h$counts)-2, "Avg", col = "black", pos = 4, srt = 45)
      
      abline(v = ten, col ="black", lwd = 3, lty = 2)
      text(ten, max(h$counts)-2, "10th", col = "black", pos = 4, srt = 45)
      
      abline(v = minimum, col ="black", lwd = 3, lty = 2)
      text(minimum, max(h$counts)-2, "Min", col = "black", pos = 4, srt = 45)
      
      
      lines(x = h$mids, y=ec(h$mids)*max(h$counts), col ='firebrick', lwd=2, axes=F)
      ticks <- c(0,10,20,30,40,50,60,70,80,90,100)
      axis(4, at=seq(from = 0, to = max(h$counts), length.out = 11), labels=seq(0, 100, 10), col = 'firebrick', col.axis = 'firebrick', las = 1)
      mtext(side = 4, line = 2, "Cumulative Frequency", col = 'firebrick')
      box()
      
      
      
    })
#######################################################################################################
#Using losses at the loss percentile    
    output$GuagePlot <- renderPlot({
      ALEsamples <- simulate()
      VARMoney <- quantile(ALEsamples, probs=(DataLossScale$LP/100))
      
      if(VARMoney <= DataLossScale$VLH){
        rating <- 8.4
      }
      else if(VARMoney > DataLossScale$VLH && VARMoney <= DataLossScale$LH){
        rating <- 25.1
      }
      else if(VARMoney > DataLossScale$LH && VARMoney <= DataLossScale$MH){
        rating <- 41.8
      }
      else if(VARMoney > DataLossScale$MH && VARMoney <= DataLossScale$SGH){
        rating <- 58.5
      }
      else if(VARMoney > DataLossScale$SGH && VARMoney <= DataLossScale$HH){
        rating <- 75.2
      }
      else {
        rating <- 91.9
      }
      
      
      gg.gauge(rating,breaks=c(0,16.7,33.4,50.1,66.8,83.5,100))
      
    })
    
#######################################################################################################    
  
  output$plot <- renderPlot({
    
    ALEsamples <- simulate()

    ALEDotted <- data.frame(
      PLMS = PLMsamples,
      SLMS = SLMsamples,
      PLEF = PLEFsamples,
      SLEF = SLEFsamples,
      chkTC = input$chkTC,
      chkRS = input$chkRS
      )
    

    #gg <- ggplot(data.frame(ALEsamples), aes(x = ALEsamples))
    #gg <- gg + ggtitle("Annualized Risk Exposure")
    #gg <- gg + ylab("Frequency")
    #gg <- gg + xlab("Loss")
    #gg <- gg + geom_histogram(color = "black", fill = "steelblue", binwidth = diff(range(ALEsamples)/50))
    #gg <- gg + scale_x_continuous(labels = comma)
    #gg <- gg + theme_bw()
    #scenario1p1 <- gg

    #gg1 <- ggplot(data.frame(ALEsamples), aes(x = ALEsamples))
    #gg1 <- gg1 + ggtitle("Annualized Risk Exposure (Financial)")
    #gg1 <- gg1 + ylab("Cumulative Frequency")
    #gg1 <- gg1 + xlab("Loss")
    #gg1 <- gg1 + stat_ecdf(color="firebrick")
    #gg1 <- gg1 + theme_bw()
    #scenario1p2 <- gg1

    gg2 <- ggplot(ALEDotted, aes(x = PLEF))
    gg2 <- gg2 + ggtitle("Risk")
    gg2 <- gg2 + ylab("Loss Magnitude")
    gg2 <- gg2 + xlab("Loss Event Frequency (Loss Events/year)")
    gg2 <- gg2 + geom_point(aes(x = SLEF, y = SLMS, colour="Secondary"), na.rm=TRUE, position = "jitter", color = "firebrick")
    gg2 <- gg2 + geom_point(aes(x = PLEF, y = PLMS, colour="Primary"), na.rm=TRUE, position = "jitter", color = "steelblue")
    gg2 <- gg2 + scale_x_log10(limits = c(0.01, 10000), labels = comma)
    gg2 <- gg2 + coord_cartesian(xlim=c(1,10E3))
    gg2 <- gg2 + scale_y_log10(labels = comma)
    gg2 <- gg2 + coord_cartesian(ylim=c(1,10E7))
    gg2 <- gg2 + theme_bw()
    scenario1p3 <- gg2

    if(ALEDotted$chkTC && ALEDotted$chkRS){
      # Change after Code Review - 2015
      # Note that the simulation has been removed and replaced by a call to the simulations run earlier.
      # This allows the plot to remain consistent with the overall model output but caveats should still be presented that this is the results of one simulation of many.
      TCPlotsamples <- TCsamples[,1]
      RSPlotsamples <- RSsamples[,1]
      TCPlot <- data.frame(Vulnerability = "Threat Capability", Continuum = TCPlotsamples)
      RSPlot <- data.frame(Vulnerability = "Control Strength", Continuum = RSPlotsamples)
      datt <- rbind(TCPlot, RSPlot)
      #datt <- data.frame(Vulnerability = factor(rep(c("Threat Capability", "Control Strength"), each=10000)), TCRSValues = c(TCPlotsamples, RSPlotsamples))
      gg3 <- ggplot(datt,aes(x = Vulnerability, y = Continuum, fill=Vulnerability))
      gg3 <- gg3 + ggtitle("Vulnerability")
      gg3 <- gg3 + xlab("")
      gg3 <- gg3 + geom_boxplot()
      gg3 <- gg3 + coord_cartesian(ylim=c(1,100))
      gg3 <- gg3 + coord_flip()
      gg3 <- gg3 + theme_bw() +  theme(legend.position = "none")
      scenario1p4 <- gg3
      #multiplot(scenario1p1, scenario1p3, scenario1p2, scenario1p4, cols=2)
      multiplot(scenario1p3, scenario1p4, cols=2)
    }    else{
      #multiplot(scenario1p1, scenario1p3, scenario1p2, cols=2)
      multiplot(scenario1p3, cols=2)
    }

   #multiplot(scenario1p1, scenario1p3, scenario1p2, cols=2)
    
  })

  output$detail2 <- renderPrint({
    ALEsamples <- simulate()
    LossPercentile <- DataLossScale$LP/100
    VARMoney <- quantile(ALEsamples, probs=(LossPercentile))
    VARServiceDownTime <- quantile(ALEServiceDowntime, probs=(LossPercentile))
    VARDataLoss <- quantile(ALEDataLoss, probs=(LossPercentile))
    print(paste0("Financial Losses at the ", DataLossScale$LP, "th percentile are estimated at £", format(round_any(VARMoney, 10), nsmall = 0, big.mark = ",")));
    print(paste0("Service Downtime Losses at the ", DataLossScale$LP, "th percentile are estimated at ", format(round_any(VARServiceDownTime, 1), nsmall = 0, big.mark = ","), " hours"));
    print(paste0("Data Record Losses at the ", DataLossScale$LP, "th percentile are estimated at ", format(round_any(VARDataLoss, 1), nsmall = 0, big.mark = ","), " KB"));
  })
#############################################################################################################
    return(TotalALE)
  })

##############################################################################

simulate1 <- eventReactive(input$Analyse, {
#simulate1 <- reactive( {
    
    N <- input$noOfSimulations
    set.seed(input$seed)


    PLEFestimate1 <- data.frame(L = input$PLEFMin1,  ML = input$PLEFML1,  H = input$PLEFMax1, CONF = as.numeric(input$PLEFConf1))
    TEFestimate1 <- data.frame(L = input$TEFMin1,  ML = input$TEFML1,  H = input$TEFMax1, CONF = as.numeric(input$TEFConf1))
    CFestimate1 <- data.frame(L = input$CFMin1,  ML = input$CFML1,  H = input$CFMax1, CONF = as.numeric(input$CFConf1))
    POAestimate1 <- data.frame(L = input$POAMin1,  ML = input$POAML1,  H = input$POAMax1, CONF = as.numeric(input$POAConf1))
    Vulnestimate1 <- data.frame(L = input$VulnMin1,  ML = input$VulnML1,  H = input$VulnMax1, CONF = as.numeric(input$VulnConf1))
    TCestimate1 <- data.frame(L = input$TCMin1,  ML = input$TCML1,  H = input$TCMax1, CONF = as.numeric(input$TCConf1))
    RSestimate1 <- data.frame(L = input$RSMin1,  ML = input$RSML1,  H = input$RSMax1, CONF = as.numeric(input$RSConf1))
    SLEFestimate1 <- data.frame(L = input$SLEFMin1,  ML = input$SLEFML1,  H = input$SLEFMax1, CONF = as.numeric(input$SLEFConf1))
    PrimaryProductivityestimate1 <- data.frame(L = input$PrimaryProductivityMin1,  ML = input$PrimaryProductivityML1,  H = input$PrimaryProductivityMax1, CONF = as.numeric(input$PrimaryProductivityConf1))
    PrimaryResponseestimate1 <- data.frame(L = input$PrimaryResponseMin1,  ML = input$PrimaryResponseML1,  H = input$PrimaryResponseMax1, CONF = as.numeric(input$PrimaryResponseConf1))
    PrimaryReplacementestimate1 <- data.frame(L = input$PrimaryReplacementMin1,  ML = input$PrimaryReplacementML1,  H = input$PrimaryReplacementMax1, CONF = as.numeric(input$PrimaryReplacementConf1))
    PrimaryFandJestimate1 <- data.frame(L = input$PrimaryFandJMin1,  ML = input$PrimaryFandJML1,  H = input$PrimaryFandJMax1, CONF = as.numeric(input$PrimaryFandJConf1))
    PrimaryCompAdvestimate1 <- data.frame(L = input$PrimaryCompAdvMin1,  ML = input$PrimaryCompAdvML1,  H = input$PrimaryCompAdvMax1, CONF = as.numeric(input$PrimaryCompAdvConf1))
    PrimaryReputationestimate1 <- data.frame(L = input$PrimaryReputationMin1,  ML = input$PrimaryReputationML1,  H = input$PrimaryReputationMax1, CONF = as.numeric(input$PrimaryReputationConf1))
    SecondaryProductivityestimate1 <- data.frame(L = input$SecondaryProductivityMin1,  ML = input$SecondaryProductivityML1,  H = input$SecondaryProductivityMax1, CONF = as.numeric(input$SecondaryProductivityConf1))
    SecondaryResponseestimate1 <- data.frame(L = input$SecondaryResponseMin1,  ML = input$SecondaryResponseML1,  H = input$SecondaryResponseMax1, CONF = as.numeric(input$SecondaryResponseConf1))
    SecondaryReplacementestimate1 <- data.frame(L = input$SecondaryReplacementMin1,  ML = input$SecondaryReplacementML1,  H = input$SecondaryReplacementMax1, CONF = as.numeric(input$SecondaryReplacementConf1))
    SecondaryFandJestimate1 <- data.frame(L = input$SecondaryFandJMin1,  ML = input$SecondaryFandJML1,  H = input$SecondaryFandJMax1, CONF = as.numeric(input$SecondaryFandJConf1))
    SecondaryCompAdvestimate1 <- data.frame(L = input$SecondaryCompAdvMin1,  ML = input$SecondaryCompAdvML1,  H = input$SecondaryCompAdvMax1, CONF = as.numeric(input$SecondaryCompAdvConf1))
    SecondaryReputationestimate1 <- data.frame(L = input$SecondaryReputationMin1,  ML = input$SecondaryReputationML1,  H = input$SecondaryReputationMax1, CONF = as.numeric(input$SecondaryReputationConf1))
    ServiceDowntimeestimate1 <- data.frame(L = input$ServiceDowntimeMin1,  ML = input$ServiceDowntimeML1,  H = input$ServiceDowntimeMax1, CONF = as.numeric(input$ServiceDowntimeConf1))
    DataLossestimate1 <- data.frame(L = input$DataLossMin1,  ML = input$DataLossML1,  H = input$DataLossMax1, CONF = as.numeric(input$DataLossConf1))
    DataLossScale1 <- data.frame(SVL = input$SVLow, HL = input$HLow, HH = input$HHigh, SGL = input$SGLow, SGH = input$SGHigh, ML = input$MLow, MH = input$MHigh, LL = input$LLow, LH = input$LHigh, VLH = input$VLHigh, LP = input$lossPercentile)

    #Run rpert calculation under different seed for vulnerability simulation
    if(input$chkRS1 && input$chkTC1){
      set.seed(input$vulnSeed)
      TCsamples1 <- sapply(rep(NA,N),function(x){rpert(N, TCestimate1$L/100, TCestimate1$ML/100, TCestimate1$H/100, shape = TCestimate1$CONF)})
      RSsamples1 <- sapply(rep(NA,N),function(x){rpert(N, RSestimate1$L/100, RSestimate1$ML/100, RSestimate1$H/100, shape = RSestimate1$CONF)})
    }
    #set seed for rest of model
    set.seed(input$seed)
    
    if(input$chkPOA1 && input$chkCF1){
      CFsamples1 <- rpert(N, CFestimate1$L, CFestimate1$ML, CFestimate1$H, shape = CFestimate1$CONF)
      POAsamples1 <- rpert(N, POAestimate1$L, POAestimate1$ML, POAestimate1$H, shape = POAestimate1$CONF)
      TEFsamples1 <- POAsamples1 * CFsamples1
    } 

    if(input$chkTEF1){
      TEFsamples1 <- rpert(N, TEFestimate1$L, TEFestimate1$ML, TEFestimate1$H, shape = TEFestimate1$CONF)
    }

    # Change after Code Review - 2015
    # Note that Vulnerability Monte Carlo simulation has been simplified to a vector calculation.
    # The colSums() function operates down the colum to collapse an NxN matrix to an Nx1 matrix
    if(input$chkRS1 && input$chkTC1){
      vuln1 <- TCsamples1 > RSsamples1
      Vulnsamples1 = (colSums(vuln1)/N) #* 100
      PLEFsamples1 <- TEFsamples1 * Vulnsamples1
    }

    if(input$chkVuln1){
      Vulnsamples1 <- rpert(N, Vulnestimate1$L/100, Vulnestimate1$ML/100, Vulnestimate1$H/100, shape = Vulnestimate1$CONF)
      PLEFsamples1 <- TEFsamples1 * Vulnsamples1
    }

    if(input$chkPLEF1){
      PLEFsamples1 <- rpert(N, PLEFestimate1$L, PLEFestimate1$ML, PLEFestimate1$H, shape = PLEFestimate1$CONF)
    }

    PrimaryProductivitysamples1 <- rpert(N, PrimaryProductivityestimate1$L, PrimaryProductivityestimate1$ML, PrimaryProductivityestimate1$H, shape = PrimaryProductivityestimate1$CONF)
    PrimaryResponsesamples1 <- rpert(N, PrimaryResponseestimate1$L, PrimaryResponseestimate1$ML, PrimaryResponseestimate1$H, shape = PrimaryResponseestimate1$CONF)
    PrimaryReplacementsamples1 <- rpert(N, PrimaryReplacementestimate1$L, PrimaryReplacementestimate1$ML, PrimaryReplacementestimate1$H, shape = PrimaryReplacementestimate1$CONF)
    PrimaryFandJsamples1 <- rpert(N, PrimaryFandJestimate1$L, PrimaryFandJestimate1$ML, PrimaryFandJestimate1$H, shape = PrimaryFandJestimate1$CONF)
    PrimaryCompAdvsamples1 <- rpert(N, PrimaryCompAdvestimate1$L, PrimaryCompAdvestimate1$ML, PrimaryCompAdvestimate1$H, shape = PrimaryCompAdvestimate1$CONF)
    PrimaryReputationsamples1 <- rpert(N, PrimaryReputationestimate1$L, PrimaryReputationestimate1$ML, PrimaryReputationestimate1$H, shape = PrimaryReputationestimate1$CONF)
    PLMsamples1 <- PrimaryProductivitysamples1 + PrimaryResponsesamples1 + PrimaryReplacementsamples1 + PrimaryFandJsamples1 + PrimaryCompAdvsamples1 + PrimaryReputationsamples1


    SecondaryProductivitysamples1 <- rpert(N, SecondaryProductivityestimate1$L, SecondaryProductivityestimate1$ML, SecondaryProductivityestimate1$H, shape = SecondaryProductivityestimate1$CONF)
    SecondaryResponsesamples1 <- rpert(N, SecondaryResponseestimate1$L, SecondaryResponseestimate1$ML, SecondaryResponseestimate1$H, shape = SecondaryResponseestimate1$CONF)
    SecondaryReplacementsamples1 <- rpert(N, SecondaryReplacementestimate1$L, SecondaryReplacementestimate1$ML, SecondaryReplacementestimate1$H, shape = SecondaryReplacementestimate1$CONF)
    SecondaryFandJsamples1 <- rpert(N, SecondaryFandJestimate1$L, SecondaryFandJestimate1$ML, SecondaryFandJestimate1$H, shape = SecondaryFandJestimate1$CONF)
    SecondaryCompAdvsamples1 <- rpert(N, SecondaryCompAdvestimate1$L, SecondaryCompAdvestimate1$ML, SecondaryCompAdvestimate1$H, shape = SecondaryCompAdvestimate1$CONF)
    SecondaryReputationsamples1 <- rpert(N, SecondaryReputationestimate1$L, SecondaryReputationestimate1$ML, SecondaryReputationestimate1$H, shape = SecondaryReputationestimate1$CONF)
    SLMsamples1 <- SecondaryProductivitysamples1 + SecondaryResponsesamples1 + SecondaryReplacementsamples1 + SecondaryFandJsamples1 + SecondaryCompAdvsamples1 + SecondaryReputationsamples1


    SLEFsamples1 <- (rpert(N, SLEFestimate1$L, SLEFestimate1$ML, SLEFestimate1$H, shape = SLEFestimate1$CONF) / 100) * PLEFsamples1

    PrimaryALE1 <- PLEFsamples1 * PLMsamples1
    SecondaryALE1 <- SLEFsamples1 * SLMsamples1
    TotalALE1 <- PrimaryALE1 + SecondaryALE1

#Added Service Downtime and Data Loss
    ServiceDowntimesamples1 <- rpert(N, ServiceDowntimeestimate1$L, ServiceDowntimeestimate1$ML, ServiceDowntimeestimate1$H, shape = ServiceDowntimeestimate1$CONF)
    DataLosssamples1 <- rpert(N, DataLossestimate1$L, DataLossestimate1$ML, DataLossestimate1$H, shape = DataLossestimate1$CONF)
    ALEServiceDowntime1 <- PLEFsamples1 * ServiceDowntimesamples1
    ALEDataLoss1 <- PLEFsamples1 * DataLosssamples1


##############################################################################

#########################################################################################################

    
    output$ProbPlot1 <- renderPlot({
  
      ALEsamples1 <- simulate1()
      
      ninety <- quantile(ALEsamples1, probs=(0.9))
      maximum <- quantile(ALEsamples1, probs = (1))
      ten <- quantile(ALEsamples1, probs = (0.1))
      average <- mean(ALEsamples1)
      minimum <- quantile(ALEsamples1, probs = (0))
      LossPercentile <- quantile(ALEsamples1, probs=(DataLossScale1$LP/100))
      
      h <- hist(ALEsamples1, breaks=100, col="steelblue", main="Annualised Risk Exposure", xlab="Loss", axes = FALSE, ylab="Simulation Distribution", panel.first=grid(NULL,NULL,lty=1))
      
      par(new = T)
      
      ec <- ecdf(ALEsamples1)
      #plot(x = h$mids, y=ec(h$mids)*max(h$counts), col = rgb(0,0,0,alpha=0), axes=F, xlab=NA, ylab=NA)
      axis(1, at=axTicks(1), labels=sprintf("£%s", axTicks(1)))
      
      abline(v = LossPercentile, col ="red", lwd = 3, lty = 2)
      
      abline(v = ninety, col ="black", lwd = 3, lty = 2)
      text(ninety, max(h$counts)-2, "90th", col = "black", pos = 4, srt = 45)
      
      abline(v = maximum, col ="black", lwd = 3, lty = 2)
      text(maximum, max(h$counts)-2, "Max", col = "black", pos = 4, srt = 45)
      
      abline(v = average, col ="black", lwd = 3, lty = 2)
      text(average, max(h$counts)-2, "Avg", col = "black", pos = 4, srt = 45)
      
      abline(v = ten, col ="black", lwd = 3, lty = 2)
      text(ten, max(h$counts)-2, "10th", col = "black", pos = 4, srt = 45)
      
      abline(v = minimum, col ="black", lwd = 3, lty = 2)
      text(minimum, max(h$counts)-2, "Min", col = "black", pos = 4, srt = 45)
      
      
      
      
      lines(x = h$mids, y=ec(h$mids)*max(h$counts), col ='firebrick', lwd=2)
      ticks <- c(0,10,20,30,40,50,60,70,80,90,100)
      axis(4, at=seq(from = 0, to = max(h$counts), length.out = 11), labels=seq(0, 100, 10), col = 'firebrick', col.axis = 'firebrick', las = 1)
      mtext(side = 4, line = 2, "Cumulative Frequency", col = 'firebrick')
      box()
      

    })
#######################################################################################################
#Using losses at the loss percentile    
    output$GuagePlot1 <- renderPlot({
      ALEsamples1 <- simulate1()
      VARMoney <- quantile(ALEsamples1, probs=(DataLossScale1$LP/100))
      
      if(VARMoney <= DataLossScale1$VLH){
        rating <- 8.4
      }
      else if(VARMoney > DataLossScale1$VLH && VARMoney <= DataLossScale1$LH){
        rating <- 25.1
      }
      else if(VARMoney > DataLossScale1$LH && VARMoney <= DataLossScale1$MH){
        rating <- 41.8
      }
      else if(VARMoney > DataLossScale1$MH && VARMoney <= DataLossScale1$SGH){
        rating <- 58.5
      }
      else if(VARMoney > DataLossScale1$SGH && VARMoney <= DataLossScale1$HH){
        rating <- 75.2
      }
      else {
        rating <- 91.9
      }
      
      
      gg.gauge(rating,breaks=c(0,16.7,33.4,50.1,66.8,83.5,100))
      
    })
    
#######################################################################################################    
  
  output$plot1 <- renderPlot({
    
    ALEsamples1 <- simulate1()

    ALEDotted1 <- data.frame(
      PLMS1 = PLMsamples1,
      SLMS1 = SLMsamples1,
      PLEF1 = PLEFsamples1,
      SLEF1 = SLEFsamples1,
      chkTC1 = input$chkTC1,
      chkRS1 = input$chkRS1
      )

    # gg <- ggplot(data.frame(ALEsamples1), aes(x = ALEsamples1))
    # gg <- gg + ggtitle("Annualized Risk Exposure")
    # gg <- gg + ylab("Frequency")
    # gg <- gg + xlab("Loss")
    # gg <- gg + geom_histogram(color = "black", fill = "steelblue", binwidth = diff(range(ALEsamples1)/50))
    # gg <- gg + scale_x_continuous(labels = comma)
    # gg <- gg + theme_bw()
    # scenario2p1 <- gg

    # gg1 <- ggplot(data.frame(ALEsamples1), aes(x = ALEsamples1))
    # gg1 <- gg1 + ggtitle("Annualized Risk Exposure (Financial)")
    # gg1 <- gg1 + ylab("Cumulative Frequency")
    # gg1 <- gg1 + xlab("Loss")
    # gg1 <- gg1 + stat_ecdf(color="firebrick")
    # gg1 <- gg1 + theme_bw()
    # scenario2p2 <- gg1

    gg2 <- ggplot(ALEDotted1, aes(x = PLEF1))
    gg2 <- gg2 + ggtitle("Risk")
    gg2 <- gg2 + ylab("Loss Magnitude")
    gg2 <- gg2 + xlab("Loss Event Frequency (Loss Events/year)")
    gg2 <- gg2 + geom_point(aes(x = SLEF1, y = SLMS1, colour="Secondary"), na.rm=TRUE, position = "jitter", color = "firebrick")
    gg2 <- gg2 + geom_point(aes(x = PLEF1, y = PLMS1, colour="Primary"), na.rm=TRUE, position = "jitter", color = "steelblue")
    gg2 <- gg2 + scale_x_log10(limits = c(0.01, 10000), labels = comma)
    gg2 <- gg2 + coord_cartesian(xlim=c(1,10E3))
    gg2 <- gg2 + scale_y_log10(labels = comma)
    gg2 <- gg2 + coord_cartesian(ylim=c(1,10E7))
    gg2 <- gg2 + theme_bw()
    scenario2p3 <- gg2

    if(ALEDotted1$chkTC && ALEDotted1$chkRS){
      # Change after Code Review - 2015
      # Note that the simulation has been removed and replaced by a call to the simulations run earlier.
      # This allows the plot to remain consistent with the overall model output but caveats should still be presented that this is the results of one simulation of many.
      TCPlotsamples1 <- TCsamples1[,1]
      RSPlotsamples1 <- RSsamples1[,1]
      TCPlot1 <- data.frame(Vulnerability1 = "Threat Capability", Continuum = TCPlotsamples1)
      RSPlot1 <- data.frame(Vulnerability1 = "Control Strength", Continuum = RSPlotsamples1)
      datt1 <- rbind(TCPlot1, RSPlot1)
      gg3 <- ggplot(datt1,aes(x = Vulnerability1, y = Continuum, fill=Vulnerability1))
      gg3 <- gg3 + ggtitle("Vulnerability")
      gg3 <- gg3 + xlab("")
      gg3 <- gg3 + geom_boxplot()
      gg3 <- gg3 + coord_cartesian(ylim=c(1,100))
      gg3 <- gg3 + coord_flip()
      gg3 <- gg3 + theme_bw() +  theme(legend.position = "none")
      scenario2p4 <- gg3
      # multiplot(scenario2p1, scenario2p3, scenario2p2, scenario2p4, cols=2)
      multiplot(scenario2p3, scenario2p4, cols=2)
    }    else{
      multiplot(scenario2p3, cols=2)
    }



  })
  

  output$detail4 <- renderPrint({
    ALEsamples1 <- simulate1()
    VARMoney1 <- quantile(ALEsamples1, probs=(0.95))
    VARServiceDownTime1 <- quantile(ALEServiceDowntime1, probs=(0.95))
    VARDataLoss1 <- quantile(ALEDataLoss1, probs=(0.95))
    print(paste0("Financial Losses at the 95th percentile are estimated at £", format(round_any(VARMoney1, 100), nsmall = 0, big.mark = ",")));
    print(paste0("Service Downtime Losses at the 95th percentile are estimated at ", format(round_any(VARServiceDownTime1, 1), nsmall = 0, big.mark = ","), " hours"));
    print(paste0("Data Record Losses at the 95th percentile are estimated at ", format(round_any(VARDataLoss1, 1), nsmall = 0, big.mark = ","), " KB"));
  })

  #############################################################################################################

    return(TotalALE1)

  })





#######################################################################################################


     output$detail <- renderPrint({          
     ALEsamples <- simulate()
     print(summary(ALEsamples));
   })
  


#######################################################################################################
  
   output$detail3 <- renderPrint({          
     ALEsamples1 <- simulate1()
     print(summary(ALEsamples1));
   })
  

##########################################################################################################


output$plot2 <- renderPlot({
    Scenario1 <- simulate()
    Scenario2 <- simulate1()
    Scenario1Plot <- data.frame(Scenario = "Scenario 1", ALE = Scenario1)
    Scenario2Plot <- data.frame(Scenario = "Scenario 2", ALE = Scenario2)
    dat <- rbind(Scenario1Plot, Scenario2Plot)
    
    
    Comparison <- data.frame(
      S1 = Scenario1,
      S2 = Scenario2,
      cont = seq(0, 100, length = 10)
    )
    
    
    
    gg1 <- ggplot(dat, aes(x=Scenario, y=ALE, fill=Scenario))
    gg1 <- gg1 + geom_boxplot()
    gg1 <- gg1 + coord_cartesian(ylim=c(1,100))
    gg1 <- gg1 + coord_flip()
    gg1 <- gg1 + theme_bw() +  theme(legend.position = "none")
    gg1 <- gg1 + xlab("")
    gg1 <- gg1 + ylab("")
    gg1 <- gg1 + scale_y_continuous(labels=dollar_format(prefix="£"))
    Compare1 <- gg1
    
    
    gg2 <- ggplot(Comparison, aes(x = cont))
    gg2 <- gg2 + geom_density(aes(x = S1, colour="Scenario 1"), alpha = 1/3, fill = "red")
    gg2 <- gg2 + geom_density(aes(x = S2, colour="Scenario 2"), alpha = 1/3, fill = "green")
    gg2 <- gg2 + theme_bw()
    gg2 <- gg2 + theme(axis.text.y=element_blank(),
                       axis.ticks.y=element_blank())
    gg2 <- gg2 + ylab("Simulation Distribution")
    gg2 <- gg2 + xlab("")
    gg2 <- gg2 + scale_x_continuous(labels=dollar_format(prefix="£"))
    Compare2 <- gg2
    
    
    multiplot(Compare2, Compare1, cols=1)
    
    #print(gg)

  })


})



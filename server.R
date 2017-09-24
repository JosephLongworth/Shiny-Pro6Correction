################################################################################################
################################################################################################
################################################################################################
# This is the ui script for the Proline 6 correction Shiny app developed by Joseph Longworth
################################################################################################
################################################################################################
################################################################################################
# ----- Load packages required
library(shiny)
library(stringr)
# ----- Increase the size of accepted upload
options(shiny.maxRequestSize=50*1024^2)

shinyServer(function(input, output) {
  output$plot <- renderPlot({
    
# ----- Define the files being uploaded the proteinGroups, evidence, and msmsScans
  inFile1 <- input$file1
  inFile2 <- input$file2

# ----- These statments control what happens until the three files have been uploaded 
  if (is.null(inFile1))
    return(NULL)
  if (is.null(inFile2))
    return(NULL)

# ----- Create the Progress bar and start
  progress <- shiny::Progress$new()
  progress$set(message = "1/6 Upload files", value = 0)  

# ----- Sets the incomming data as data.frames
  ProteinGroups=  read.delim(inFile1$datapath,sep='\t')
  evidence=  read.delim(inFile2$datapath,sep='\t')

  # rm (list=ls())
  # setwd("E:/Google Drive/Shiny_apps/P6_Correction_V5")
  # ProteinGroups=read.delim("test_data/proteinGroups.txt",sep='\t')
  # evidence=  read.delim("test_data/evidence.txt",sep='\t')
  # 
# -----  Clear Contaminants, Decoy and Only identified by site
  ProteinGroups=ProteinGroups[!ProteinGroups$Reverse=="+",]
  ProteinGroups=ProteinGroups[!ProteinGroups$Potential.contaminant=="+",]
  ProteinGroups=ProteinGroups[!ProteinGroups$Only.identified.by.site=="+",]
  evidence=evidence[!evidence$Reverse=="+",]
  evidence=evidence[!evidence$Potential.contaminant=="+",]
  
# -----  Use this line to remove the Pro6  from evidence. Requires the column head "Pro6"
  evidence=evidence[evidence$Pro6==0,]
  
  
# -----  Count the number of Proline in each evidence sequence and split into sub df
progress$set(message = "2/6 Count Prolines", value = 20)  
  p=c()
  for(i in 1:nrow(evidence)){p=c(p,str_count(as.character(evidence$Sequence[i]),"P"))}
  df=cbind(evidence,p)
  df_split=split(df,df$p)
  for (i in 1:length(df_split)){
    df_temp=as.data.frame(df_split[i])
    colnames(df_temp)=colnames(df)
    assign(df_temp,x = paste0("df",i))}
  
  
  # -----  Plot P=0 and P=1 to determine correction factors
  boxplot=boxplot(df1$Ratio.H.L,df2$Ratio.H.L,outline=F)
  correction_factor=(boxplot$stats[3,1]/(boxplot$stats[3,2]))-1
  boxplot$stats
  a=median(df1$Ratio.H.L,na.rm = T)
  b=median(df2$Ratio.H.L,na.rm = T)
  
  a/b-1
  
  # -----  Correct evidence file and split into sub df
progress$set(message = "3/6 Correct Ratios", value = 40)  
  corrected_evidence=evidence
 corrected_evidence$Ratio.H.L=evidence$Ratio.H.L*(1+(correction_factor)*p)
 df=cbind(corrected_evidence,p)
  df_split=split(df,df$p)
  for (i in 1:length(df_split)){
    df_temp=as.data.frame(df_split[i])
    colnames(df_temp)=colnames(df)
    assign(df_temp,x = paste0("df_corrected",i))}
  
  
  # -----  Plot 4 graphs with and without correction and with and without outliers
progress$set(message = "4/6 Plot boxplots", value = 60) 
par(mfrow=c(2,2))
  boxplot_pre=boxplot(log10(df1$Ratio.H.L),log10(df2$Ratio.H.L),log10(df3$Ratio.H.L), 
                      log10(df4$Ratio.H.L),log10(df5$Ratio.H.L),log10(df6$Ratio.H.L),
                      log10(df7$Ratio.H.L),outline=T,main="Uncorrected with outliers",
                      names = c("0","1","2","3","4","5","6"),xlab="Proline count",ylab="H/L ratio")
  
  boxplot_pre=boxplot(log10(df1$Ratio.H.L),log10(df2$Ratio.H.L),log10(df3$Ratio.H.L),
                      log10(df4$Ratio.H.L),log10(df5$Ratio.H.L),log10(df6$Ratio.H.L),
                      log10(df7$Ratio.H.L),outline=F,main="Uncorrected without outliers",
                      names = c("0","1","2","3","4","5","6"),xlab="Proline count",ylab="H/L ratio")
  
  boxplot_pre=boxplot(log10(df_corrected1$Ratio.H.L),log10(df_corrected2$Ratio.H.L),
                      log10(df_corrected3$Ratio.H.L),log10(df_corrected4$Ratio.H.L),
                      log10(df_corrected5$Ratio.H.L),log10(df_corrected6$Ratio.H.L),
                      log10(df_corrected7$Ratio.H.L),outline=T,main="Corrected with outliers",
                      names = c("0","1","2","3","4","5","6"),xlab="Proline count",ylab="H/L ratio",
                      sub=paste0("Correction factor = ",round(correction_factor,3)))
  
  boxplot_pre=boxplot(log10(df_corrected1$Ratio.H.L),log10(df_corrected2$Ratio.H.L),
                      log10(df_corrected3$Ratio.H.L),log10(df_corrected4$Ratio.H.L),
                      log10(df_corrected5$Ratio.H.L),log10(df_corrected6$Ratio.H.L),
                      log10(df_corrected7$Ratio.H.L),outline=F,main="Corrected without outliers",
                      names = c("0","1","2","3","4","5","6"),xlab="Proline count",ylab="H/L ratio",
                      sub=paste0("Correction factor = ",round(correction_factor,3)))

# -----  Split evidence by protein group to determin the corrected protein ratios
progress$set(message = "5/6 Compile Protein ratios", value = 80)   
IDs=corrected_evidence$Protein.group.IDs
IDs=IDs[!duplicated(IDs)]
Corrected_ProtienGroups_ratios=matrix(nrow = length(IDs),ncol = 3)
Ratio_col=which(colnames(corrected_evidence)=="Ratio.H.L")
for(i in 1:length(IDs)){
  Corrected_ProtienGroups_ratios[i,1]=as.character(IDs[i])
  Corrected_ProtienGroups_ratios[i,2]= median(corrected_evidence[corrected_evidence$Protein.group.IDs==IDs[i],Ratio_col],na.rm = T)}


# -----  Create normalised corrected ratios and merge with current Protein group file
as.numchar=function(x){as.numeric(as.character(x))}
Corrected_ProtienGroups_ratios[,3]=as.numchar(Corrected_ProtienGroups_ratios[,2])/median(as.numchar(Corrected_ProtienGroups_ratios[,2]),na.rm = T)
Corrected_ProtienGroups_ratios=as.data.frame(Corrected_ProtienGroups_ratios)
colnames(Corrected_ProtienGroups_ratios)=c("id","Corrected.Ratio.H.L","Corrected.Ratio.H.L.normalized" )
Corrected_ProteinGroups<<-merge.data.frame(ProteinGroups,Corrected_ProtienGroups_ratios,by="id",all.x = T)
  
      })

# ----- Build the UI download button preventing construction until Download is available.
output$download_button <- renderUI({
  
  # ----- Define the files being uploaded the proteinGroups, evidence, and msmsScans
  inFile1 <- input$file1
  inFile2 <- input$file2
  
  # ----- These statments control what happens until the three files have been uploaded 
  if (is.null(inFile1))
    return(NULL)
  if (is.null(inFile2))
    return(NULL)
downloadButton('downloadData', 'Download')
})

# ----- Creats the download as part of the 
output$downloadData <- downloadHandler( filename =  function() {paste("proteinGroups_ProlienCorrected","_",Sys.time(), '.csv', sep='') 
}  ,  content = function(file) {
  write.csv(Corrected_ProteinGroups,file,row.names =F)
})
})

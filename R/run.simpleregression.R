#' @title Linear modeling of drug to a binary(or categorical) covariate. With ability to control for one other covariate, output correlation of continuous value with drug sensitivity as log10 signed pvals.
#' @title writes out to a file and returns a data frame of the results.
#' @param drug.frame columns are samples drugs are rows. Data frame with 'drug' in first column header, and list of drug names underneath.Other columns have sample names as column headers and drug values in the frame itself
#' @param categorical.frame columns are samples continous values  are rows. data frame of expression or dependency etc data,  'gene' is first column header with list of gene names or feature names, and list of sample names. Other columns have sample names as column names and values (eg. expression) in the frame itself
#' @param drug.name name of the drug you will be comparing to i.e. column name to extract from the drug frame
#' @param type.frame Default null. If you want to correct for a covariate create a dataframe with headers 'sample' 'type' . sample names in first column , cancer type in second column (or whatever covariate you want to correct for)
#' @param output.file where to write the output file. default is current working directory with name categorical.with.type.signedlog10pvals.txt (nmuts will not make sense if you are not using mutation data (i.e. 1s and 0s))
#' @param percent.zeros remove rows that have more than this percent of zeros  0 to 1 scale.  0.98 (98 percent) is default
#' @param keep.na keep genes with NA in the output. default T
#' @param reverse.sign reverse sign of the output pval table, default F

#' @export
#' @importFrom mixOmics pls

run.simpleregression<-function(GeneNumber, GeneName, exp_matrix, Predictor, Transpose)
{
  output=as.data.frame(matrix(ncol=9, nrow=GeneNumber))
  colnames(output)=c("Gene","Intercept","Intercept Pval", "Predictor1 Pval", "Predictor1 Coefficient", "Predictor1 tval", "Predictor2 Pval","Predictor2 Coefficient", "Predictor 2 tval")
  output[,1]=GeneName
  
  for (i in 1:GeneNumber)
  {
    skip=F
    if (Transpose==F)
    {
      gene=as.numeric(exp_matrix[,i])
    }
    else
    {
      gene=as.numeric(exp_matrix[i,2:ncol(exp_matrix)])
      #First column of expression data is gene name
      gene=melt(gene, id.vars="gene")
    }
    merge=data.frame(gene,Predictor)
    colnames(merge)=c("Gene", "Predictor1")
    
    tryCatch({summary(lm(Gene ~ Predictor1, merge)) },
             
             error=function(e){
               
               print(e);
               
               skip<<-T;
               
             })
    
    if(skip==T)
    {
      output[i,2]=NA;
      next;
    }
    
    results= summary(lm(Gene ~ Predictor1, merge))$coefficients
    #intercept
    output[i,2]=results[1,1]
    #Intercept pvalue
    output[i,3]=results[1,4]
    #Predictor log sign pval
    output[i,4]=-log10(results[2,4])*sign(results[2,3])
    #Predictor coefficient
    output[i,5]=results[2,1]
    #Predictor t value
    output[i,6]=results[2,3]
    output[i,7]=NA
    output[i,8]=NA
    output[i,9]=NA
  }
  return(output)
}

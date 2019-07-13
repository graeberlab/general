#' @title Linear modeling of drug to a continuous covariate. With ability to control for one other covariate, output correlation of continuous value with drug sensitivity as log10 signed pvals.
#' @title writes out to a file and returns a data frame of the results.
#' @param drug.frame columns are samples drugs are rows. Data frame with 'drug' in first column header, and list of drug names underneath.Other columns have sample names as column headers and drug values in the frame itself
#' @param continuous.frame columns are samples continous values  are rows. data frame of expression or dependency etc data,  'gene' is first column header with list of gene names or feature names, and list of sample names. Other columns have sample names as column names and values (eg. expression) in the frame itself
#' @param drug.name name of the drug you will be comparing to i.e. column name to extract from the drug frame
#' @param type.frame Default null. If you want to correct for a covariate create a dataframe with headers 'sample' 'type' . sample names in first column , cancer type in second column (or whatever covariate you want to correct for)
#' @param output.file where to write the output file. default is current working directory with name continuous.with.type.signedlog10pvals.txt
#' @param percent.zeros remove rows that have more than this percent of zeros  0 to 1 scale.  0.98 (98 percent) is default
#' @export
#' @importFrom mixOmics pls

lm.one.drug.vs.continuous=function (continuous.frame,drug.frame,type.frame=NULL,drug.name="TRE515",
                                             output.file="./continuous.with.signedlog10pvals.txt",percent.zeros=.98) {
  common_samps=intersect(colnames(continuous.frame)[-1],colnames(drug.frame)[-1])
  colnames(continuous.frame)[1]="gene"
  colnames(drug.frame)[1]='drug'
  drug.frame=drug.frame[,c("drug",common_samps)]
  continuous.frame=continuous.frame[,c("gene",common_samps)]
  continuous.frame$gene=gsub(x =continuous.frame$gene,pattern = " ",replacement = "" )

  continuous.frame=continuous.frame[apply(continuous.frame[,-1],1,function(x) var(x,na.rm=T)) >0,]
  myidx=  !apply(continuous.frame[,-1],1, function (x) ((sum(na.omit(x) == 0))/length(na.omit(x)) >= percent.zeros))
  continuous.frame= continuous.frame[myidx,]
  pval = as.data.frame(matrix(nrow = nrow(continuous.frame), ncol = 2))
  colnames(pval)[1]="gene"
  colnames(pval)[2]=paste0(drug.name,"_log10pval")
  #colnames(pval)[3]="nmuts"
  # pval$nmuts=rowSums(continuous.frame[,-1])
  pval$gene = continuous.frame[,1]
  if(!is.null(type.frame)){
    colnames(type.frame)[1]="sample"
    colnames(type.frame)[2]="type"
  }

  for (i in 1:nrow(continuous.frame)){
    print(i)
    gene = continuous.frame[i,]
    tgene = reshape2::dcast(reshape2::melt(gene, id.vars = "gene"), variable ~ gene)
    colnames(tgene)[1] ="sample"

    #you could put another for loop here to go over all drugs
    drug.individual=drug.frame[drug.frame$drug==drug.name,]
    tdrug= reshape2::dcast(reshape2::melt(drug.individual, id.vars = "drug"), variable ~ drug)
    colnames(tdrug)[1] ="sample"

    cmerge = merge(tgene, tdrug, by =  c("sample"))
    colnames(cmerge)[2] ="gene"
    colnames(cmerge)[3]="drug"
    if(!is.null(type.frame)){
      cmerge=dplyr::inner_join(cmerge,type.frame,by= "sample")
    }

    cmerge=na.omit(cmerge)
    if( sum(colSums(is.na(cmerge))  %in% ncol(cmerge)) == T) { #if all NAs for any drug or gene skip
      pval[i,2]=NA
      next
    }
    tryCatch({ #if errors skip, usually because there arent enough lines with values for both measurements to run lm
       if(!is.null(type.frame)){
         tstat = sign(summary(lm(drug ~ gene + type , data = cmerge,na.action="na.omit"))$coefficients[2,3])
         pv = -1*log10(summary(lm(drug ~ gene + type , data = cmerge,na.action="na.omit"))$coefficients[2,4])
         pval[i,2] = tstat * pv
       } else {
         tstat = sign(summary(lm(drug ~ gene , data = cmerge,na.action="na.omit"))$coefficients[2,3])
         pv = -1*log10(summary(lm(drug ~ gene , data = cmerge,na.action="na.omit"))$coefficients[2,4])
         pval[i,2] = tstat * pv
       }

    }, error = function (e) {
      print (e);
      pval[i,2] = NA
    })






  }
  write.table(pval,output.file,quote=F,row.names=F,sep="\t")

  return(pval)

}




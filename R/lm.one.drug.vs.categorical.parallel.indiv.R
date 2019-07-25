#' @title Used only for running individual_cancers. Linear modeling of drug to a binary(or categorical) covariate. With ability to control for one other covariate, output correlation of continuous value with drug sensitivity as log10 signed pvals.
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

lm.one.drug.vs.categorical.parallel2=function (categorical.frame,drug.frame,type.frame=NULL,drug.name="TRE515",name.frame=NULL,
                                              output.file="./categorical.with.signedlog10pvals.txt",percent.zeros=1,keep.na=T,reverse.sign=F,covar.frame=NULL) {
  common_samps=intersect(colnames(categorical.frame)[-1],colnames(drug.frame)[-1])
  colnames(categorical.frame)[1]="gene"
  colnames(drug.frame)[1]='drug'
  drug.frame=drug.frame[,c("drug",common_samps)]
  categorical.frame=categorical.frame[,c("gene",common_samps)]
  categorical.frame$gene=gsub(x =categorical.frame$gene,pattern = " ",replacement = "" )
  categorical.frame=categorical.frame[apply(categorical.frame[,-1],1,function(x) var(x,na.rm=T)) >0,]
  myidx=  !apply(categorical.frame[,-1],1, function (x) ((sum(na.omit(x) == 0))/length(na.omit(x)) >= percent.zeros))
  categorical.frame= categorical.frame[myidx,]
  
  
  if(!is.null(type.frame)){
    colnames(type.frame)[1]="sample"
    colnames(type.frame)[2]="type"
    
    pval = as.data.frame(matrix(nrow = nrow(categorical.frame), ncol = 5),stringsAsFactors = F)
    colnames(pval)[1]="gene"
    colnames(pval)[2]="full.name"
    colnames(pval)[3]="beta"
    colnames(pval)[4]=paste0(drug.name,"_log10pval")
    colnames(pval)[5]="nmuts"
    #pval$nmuts=rowSums(categorical.frame[,-1])
    pval$gene = categorical.frame[,1]
    pval$full.name=name.frame$name[match(pval$gene, name.frame$symbol)]
    rownames(categorical.frame) = NULL #reset indexes
  } else {
    pval = as.data.frame(matrix(nrow = nrow(categorical.frame), ncol = 5),stringsAsFactors = F)
    colnames(pval)[1]="gene"
    colnames(pval)[2]="full.name"
    colnames(pval)[3]="beta"
    colnames(pval)[4]=paste0(drug.name,"_log10pval")
    colnames(pval)[5]="nmuts"
    # pval$nmuts=rowSums(categorical.frame[,-1])
    pval$gene = categorical.frame[,1]
    pval$full.name=name.frame$name[match(pval$gene, name.frame$symbol)]
  }
  
  
  do.categorical.lm= function(i,my.output=c(NA,NA,NA)){
    print(i)
    gene = categorical.frame[i,]
    tgene = data.frame(reshape2::dcast(reshape2::melt(gene, id.vars = "gene"), variable ~ gene),stringsAsFactors=F)
    colnames(tgene)[1] ="sample"
    
    #you could put another for loop here to go over all drugs
    drug.individual=drug.frame[drug.frame$drug==drug.name,]
    tdrug= data.frame(reshape2::dcast(reshape2::melt(drug.individual, id.vars = "drug"), variable ~ drug),stringsAsFactors=F)
    colnames(tdrug)[1] ="sample"
    
    cmerge = data.frame(merge(tgene, tdrug, by =  c("sample")),stringsAsFactors = F)
    colnames(cmerge)[2] ="gene"
    colnames(cmerge)[3]="drug"
    if(!is.null(type.frame)){
      cmerge=dplyr::inner_join(cmerge,type.frame,by= "sample")
    }
    if(!is.null(covar.frame)){
      cmerge=dplyr::inner_join(cmerge,covar.frame,by= "sample")
    }
    
    cmerge=na.omit(cmerge)
    cmerge$drug=as.numeric(cmerge$drug)
    cmerge$sample=as.character(cmerge$sample)
    cmerge$gene=as.numeric(cmerge$gene)
    if(!is.null(type.frame)){
      cmerge$type=as.character(cmerge$type)
    }
    my.muts=sum(cmerge$gene,na.rm=T)
    if(my.muts <=0) {return( c(NA,NA,my.muts))} # dont return 1 mutation versions.
    # if(  sum(((colSums(is.na(cmerge))  == nrow(cmerge)) > 0)) >=1) { #if all NAs for any drug or gene skip
    #   my.pval=NA
    #   return(my.pval)
    # }
    
    
    
    tryCatch({ #if errors skip, usually because there arent enough lines with values for both measurements to run lm
      if(!is.null(type.frame)){
        cmerge$type=as.character(cmerge$type)
        if(!is.null(covar.frame)){
          model.summary=summary(lm(drug ~ gene + type + covar1 , data = cmerge,na.action="na.omit"))$coefficients
          beta = model.summary[2,1]
          tstat = sign(model.summary[2,3])
          pv = -1*log10(model.summary[2,4])
          my.pval = tstat * pv
          my.output=c(beta,my.pval,my.muts)
          
          
          
        } else{
          
          
          model.summary=summary(lm(drug ~ gene  , data = cmerge,na.action="na.omit"))$coefficients
          beta = model.summary[2,1]
          tstat = sign(model.summary[2,3])
          pv = -1*log10(model.summary[2,4])
          my.pval = tstat * pv
          my.output=c(beta,my.pval,my.muts)
        }} else {
          #cmerge$gene=factor(cmerge$gene, levels= c(0,1))
          model.summary=summary(lm(drug ~ gene , data = cmerge,na.action="na.omit"))$coefficients
          beta = model.summary[2,1]
          tstat = sign(model.summary[2,3])
          pv = -1*log10(model.summary[2,4])
          my.pval = tstat * pv
          my.output=c(beta,my.pval,my.muts)
        }
      return(my.output)
    }, error = function (e) {
      print (e);
      
      my.output=c(NA,NA,NA)
      return(my.output)
    })
    
    
  }
  
  
  library(parallel)
  n_core <- detectCores()
  cl <- makeCluster(n_core-1,outfile="log.txt")
  clusterExport(cl,varlist=c('do.categorical.lm',ls()),envir=environment())
  output.vect=parSapply(cl = cl,X = 1:nrow(categorical.frame), FUN = function(x) do.categorical.lm(i=x))
  pval[,c(3:5)]=t(output.vect)
  stopCluster(cl)
  
  
  if(reverse.sign==T){
    #numeric_cols <- vapply(pval, is.numeric, logical(1))
    pval[, 3:4] <- -1* pval[,3:4]
  }
  
  if (keep.na==F){
    pval=na.omit(pval,cols=c(3,4))
  }
  
  
  # pval=pval[order(pval[,4],decreasing= F),]
  write.table(pval,output.file,quote=F,row.names=F,sep="\t")
  
  return(pval)
  
}



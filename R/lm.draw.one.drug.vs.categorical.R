
lm.draw.one.drug.vs.categorical=function (categorical.frame,drug.frame,type.frame=NULL,drug.name="TRE515",gene.name="KRAS",name.frame=NULL,add.facet=T,
                                              output.directory="./",percent.zeros=1,keep.na=T,reverse.sign=F,covar.frame=NULL,xlab="mutation status",ylab="TRE515 IC50") {
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
  
  

    print(gene.name)
    gene = categorical.frame[categorical.frame[,1] == gene.name,]
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
          pv = model.summary[2,4]
          my.pval = tstat * pv
          
          
          
          
        } else{
          
          
          model.summary=summary(lm(drug ~ gene + type , data = cmerge,na.action="na.omit"))$coefficients
          beta = model.summary[2,1]
          tstat = sign(model.summary[2,3])
          pv = model.summary[2,4]
          my.pval = tstat * pv
         
        }} else {
          #cmerge$gene=factor(cmerge$gene, levels= c(0,1))
          model.summary=summary(lm(drug ~ gene , data = cmerge,na.action="na.omit"))$coefficients
          beta = model.summary[2,1]
          tstat = sign(model.summary[2,3])
          pv = model.summary[2,4]
          my.pval = tstat * pv
         
        }
   
    }, error = function (e) {
      print (e);
      
    })
    
    
   h = ggpubr::ggscatter(data = cmerge,x = "gene", y = "drug",cor.coef = T, add='reg.line',title = gene.name,xlab=xlab,ylab=ylab) + 
     annotate("text", x=Inf, y=Inf, label= paste0("type controlled signed pval= ",round(my.pval,2)),vjust=1, hjust=1)
     ggsave(filename = paste0(output.directory, gene.name,".categorical.with.signedpvals.png"), plot = h,width = 4,height = 4)
  
     
   my.types = cmerge %>% group_by(type) %>%
          summarize(n_unique = n_distinct(gene))
   my.types=as.data.frame(my.types)
   my.types=my.types[my.types$n_unique >1,]$type
   cmerge2=cmerge[cmerge$type %in% my.types,]
   if(add.facet==T){
   i = ggpubr::ggscatter(data = cmerge2,x = "gene", y = "drug",cor.coef = T,cor.coef.size=2.8, add='reg.line',title = gene.name,xlab=xlab,ylab=ylab) + facet_wrap( ~ type,scales="free") 
   ggsave(filename = paste0(output.directory, gene.name,".facet.categorical.with.signedpvals.png"), plot = i,width = 12,height = 12)
   }
    # 
    
    
     
  
  # library(parallel)
  # n_core <- detectCores()
  # cl <- makeCluster(n_core-1,outfile="log.txt")
  # clusterExport(cl,varlist=c('do.categorical.lm',ls()),envir=environment())
  # output.vect=parSapply(cl = cl,X = 1:nrow(categorical.frame), FUN = function(x) do.categorical.lm(i=x))
  # pval[,c(3:5)]=t(output.vect)
  # stopCluster(cl)
  # 
  # 
  # if(reverse.sign==T){
  #   #numeric_cols <- vapply(pval, is.numeric, logical(1))
  #   pval[, 3:4] <- -1* pval[,3:4]
  # }
  # 
  # if (keep.na==F){
  #   pval=na.omit(pval,cols=c(3,4))
  # }
  # 
  # 
  # # pval=pval[order(pval[,4],decreasing= F),]
  # write.table(pval,output.file,quote=F,row.names=F,sep="\t")
  # 
  # return(pval)
  
}





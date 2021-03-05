plink2ccret = function(tag, gset, workdir, Rout=FALSE)
{
  #+++++This is tested with bedtools version 2.25.0 and R version 3.2.2
  #+++++https://github.com/arq5x/bedtools2/releases/bedtools-2.25.0.tar.gz
  #+++++Example: tag='mycnv1';  gset='mygeneset.bed'; workdir='./Data'
  
  setwd(workdir)
  f.fam = paste(tag, '.fam', sep='')
  f.cnv = paste(tag, '.cnv', sep='')
  
  #++++sort and index the input cnvs and write the files 
  x      = read.table(f.cnv, head=T, as.is=T)
  x      = x[order(x$CHR, x$BP1),]
  x1     = data.frame(idx=1:nrow(x), x)
  write.table(x1, file = paste(f.cnv, '.indxed.txt',sep=''), sep='\t', quote=F, row.name=F)
  x2     = x1
  x2$BP1 = x2$BP1-1   #convert 1-based to 0-based
  write.table(x2[,c(4:6,1)], file = paste(f.cnv, '.indxed.bed',sep=''), sep='\t', quote=F, row.name=F, col.name=F)
  rm(x); rm(x1); rm(x2)
  
  #+++++++++Run bedtools in order to construct CNVRs
  f.cnvindxed     = paste(f.cnv, '.indxed.bed',sep='')
  f.cnvmerged.n   = paste(f.cnv, '.indxed.merged.n.bed',sep='')
  f.cnvmerged.nms = paste(f.cnv, '.indxed.merged.nms.bed',sep='')
  f.cnv.gset=paste(f.cnv, '.indxed.', gset, sep='') 
  try(system(paste('bedtools merge -i ', f.cnvindxed, ' -c 1 -o count > ', f.cnvmerged.n)))
  try(system(paste('bedtools merge -i ', f.cnvindxed, ' -c 4 -o collapse > ', f.cnvmerged.nms)))
  try(system(paste('bedtools intersect -a ', f.cnvindxed, ' -b ', gset, ' -wa -wb > ', f.cnv.gset))) 
 
  #+++++++Note the rows in CCRET matrices are ordered by usam
  #+++++++Get yy - vector indicating affection status
  fam = read.table(f.fam, as.is=T)
  colnames(fam) = c('FID','IID','PID','MID','SEX','PHE')
  n   = nrow(fam)
  if( length(unique(fam[,1]))!=n ) {
      stop('samples are not unique')
  }    
  usam = fam[,1]
  yy   = data.frame(fam[,6])
  rownames(yy) = usam
  
  #++++++++Get the subjectbygene matrix
  allgene = read.table(gset, sep='\t', as.is=T)
  colnames(allgene) = c('gene.chr','gene.startPos','gene.endPos','gene.name')
  ngene   = nrow(allgene)  
  x       = read.table(paste(f.cnv, '.indxed.txt',sep=''), head=T, as.is=T)
  x1      = read.table(f.cnv.gset, sep='\t', as.is=T)
  colnames(x1) = c('chr','bp1','bp2','idx2','gene.chr','gene.startPos','gene.endPos','gene.name')
  xm      = merge(x, x1, by.x='idx', by.y='idx2')  #merging by CNV index attaches FIDs
  out     = matrix(NA, nrow=n, ncol=ngene)
  for (i in 1:n){
     for (j in 1:ngene){
           tmp = xm$TYPE[xm$FID==usam[i] & xm$gene.name==allgene$gene.name[j]]
           if (length(tmp) == 0){
                 out[i,j] = 0
           }
           if (length(tmp) > 0){
                 out[i,j] = 1+sum(tmp>2)   #Coding: normal:0, del:1; dup:2
           }
     } 
  }
  rownames(out) = usam
  colnames(out) = allgene$gene.name
  subjectbygene=out

  #+++++++++++
  # Verified "subjectbygene" by checking against indivual burden vector e.g.  
  # tt1=tapply(xm$gene.name, xm$FID, function(x) length(unique(x)))
  # Next construct the CNVRs and save the records for checking (cnvr.cnvs)
  #++++++++++
  
  x  = read.table(f.cnv, head=T, as.is=T)
  x  = x[order(x$CHR, x$BP1),]
  x  = data.frame(idx=1:nrow(x), x) 
  y1 = read.table(f.cnvmerged.n, sep='\t', as.is=T)
  colnames(y1) = c('cnvr.chr','cnvr.str','cnvr.stp', 'cnvr.ncnv')
  y2 = read.table(f.cnvmerged.nms, sep='\t', as.is=T)
  colnames(y2) = c('cnvr.chr','cnvr.str','cnvr.stp', 'cnv.idx')
  cnvr  = data.frame(y2, cnvr.ncnv=y1$cnvr.ncnv, cnvr.len = y1$cnvr.stp - y1$cnvr.str)
  ncnvr = nrow(cnvr)

  cnvr.cnvs = NULL
  geno      =   matrix(2, ncol=ncnvr, nrow=n)  #dosage matrix
  w.len     =   matrix(0, ncol=ncnvr, nrow=n)  #length matrix
  med.olap.cnv.cnvr = rep(NA, ncnvr)
  for (k in 1:ncnvr){
    thisidx = strsplit(y2$cnv.idx[k],",")
    thisidx = thisidx[[1]]
    ncnv    = length(thisidx)
    ulen    = cnvr$cnvr.stp[k] - cnvr$cnvr.str[k] 
    thiscnvs = x[match(thisidx, x$idx),]
    olap.cnv.cnvr = round((thiscnvs$BP2 - thiscnvs$BP1 +1)/ulen,3)
    thiscnvr = data.frame(cnvr.idx=rep(k,ncnv), cnvr[k,], thiscnvs, olap.cnv.cnvr)
    med.olap.cnv.cnvr[k] = median(olap.cnv.cnvr)
    cnvr.cnvs = rbind(cnvr.cnvs,thiscnvr)
    g      = match(thiscnvr$FID, usam)  #this gives unique iid by merging split-cnvs from the same subject
    tmpl   = tapply(thiscnvs$BP2 - thiscnvs$BP1 +1, thiscnvr$FID, sum)
    tmplg  = match(names(tmpl), usam) 
    geno[g,k] = thiscnvs$TYPE
    w.len[tmplg,k]   = tmpl
 }
 rownames(geno) = rownames(w.len) = usam
 colnames(geno) = colnames(w.len) = c(1:ncnvr)
 cnvr = data.frame(cnvr, med.olap.cnv.cnvr)

print(paste('writing output files to ', workdir))
write.table(yy,    file = paste(tag, '_yy.txt', sep=''), sep='\t', quote=F, row.name=T, col.name=T)  #yy
write.table(geno,  file = paste(tag, '_ds.txt', sep=''), sep='\t', quote=F, row.name=T, col.name=T)  #ds.mat
write.table(w.len, file = paste(tag, '_ln.txt', sep=''), sep='\t', quote=F, row.name=T, col.name=T)  #ln.mat
write.table(subjectbygene, file = paste(tag, '_gi.txt', sep=''), sep='\t', quote=F, row.name=T, col.name=T)  #gi.mat
#++optional
if (Rout == TRUE) {
   save(fam, cnvr, cnvr.cnvs, yy, geno, w.len, subjectbygene, file=paste(tag,'.rda',sep=''))
}   

} #end of the function






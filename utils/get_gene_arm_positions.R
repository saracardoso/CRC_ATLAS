chrs_genes = read.table(paste(project_dir, 'utils/gencode_v19_gene_pos.txt', sep='/'), row.names=1)
chrs_genes = chrs_genes[rownames(predictions_obs),]

chrs_arms = c()
for(chr in 1:22){
  message(chr)
  chr_url = RCurl::getURL(paste('http://atlasgeneticsoncology.org/Indexbychrom/idxg_', chr, '.html', sep=''), .opts=list(ssl.verifypeer=FALSE))
  chr_table = XML::readHTMLTable(chr_url)
  chr_table = rlist::list.clean(chr_table, fun=is.null, recursive=FALSE)
  
  if(length(chrs_arms)==0) chrs_arms = cbind(chr_table[[1]][[1]], chr_table[[1]][[3]])
  else chrs_arms = rbind(chrs_arms, cbind(chr_table[[1]][[1]], chr_table[[1]][[3]]))
}
chrs_arms = unique(chrs_arms)

location = rep(NA, dim(chrs_genes)[1])
chrs_genes = cbind(chrs_genes, location)
genes = chrs_arms[chrs_arms[,1]%in%rownames(chrs_genes),1] # Genes from chrs_genes we have info for
chrs_genes[genes,'location'] = chrs_arms[chrs_arms[,1]%in%rownames(chrs_genes), 2]
arm = gsub('[0-9]*[.]*[-]*', '', chrs_genes$location)
chrs_genes = cbind(chrs_genes, arm)
chrs_genes$arm = substr(chrs_genes$arm, 1, 1)

sum(is.na(chrs_genes$arm)) # Nº genes we do not have info for: 476
for(chr in unique(chrs_genes$V2)){
  message(chr)
  chr_se = range(which(chrs_genes$V2==chr))
  p_e = max(which(chrs_genes$arm=='p' & chrs_genes$V2==chr))
  q_se = range(which(chrs_genes$arm=='q' & chrs_genes$V2==chr))
  
  nas = which(is.na(chrs_genes$arm))
  if(!p_e%in%c(-Inf, Inf)){
    nas_p = nas[nas>chr_se[1] & nas<p_e]
    chrs_genes$arm[nas_p] = 'p'
  } 
  nas_q = nas[nas>q_se[1] & q_se[2]]
  chrs_genes$arm[nas_q] = 'q'
}
sum(is.na(chrs_genes$arm))
chrs_genes$arm[which(is.na(chrs_genes$arm))] = 'p'
write.table(chrs_genes, paste(project_dir, 'utils/gene_pos_wArms.txt', sep='/'))
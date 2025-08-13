df <- data.frame()
for ( i in list){
d <- read.table(i,sep="\t",header=T)

hla_file <- gsub("_rsem.genes.results","",i)
hla_file <- gsub("./","",hla_file)

hla <- read.table(paste("/BiO/Hyein/OptiType/",hla_file,"/",hla_file,"_expr_result.tsv",sep=""),sep="\t",header=T)

d$RPK <- ifelse(d$effective_length > 0,
                d$expected_count / (d$effective_length / 1000),
                0)
if (nrow(hla) != 6) {

  table_file <- paste("/BiO/Hyein/OptiType/", hla_file, "/", hla_file, "_table_result.tsv", sep = "")
  table <- read.table(table_file, sep = "\t", header = TRUE)
  table <- table %>%
  filter(!grepl("^[0-9.]+$", Original)) %>%  # Original이 숫자(소수 포함)인 행 제거
  filter(!grepl("^[0-9.]+$", HLA_4digit))  

  ## 여기서 table을 이용해 동형 locus 찾기 & allele ID 붙이기
  ## 예: table 컬럼 Original=HLA000xx, HLA_4digit=A*xx:xx
  hla <- hla %>%
    left_join(table, by = c("Allele" = "Original"))

hla <- hla %>%
  mutate(locus = substr(HLA_4digit, 1, 1))  

expand_homozygous <- function(df_locus) {
  if (nrow(df_locus) == 1) {
    out <- df_locus[rep(1, 2), , drop = FALSE]  # 행 복제
    out$Allele  <- paste0(out$Allele,  c("_1","_2"))
    # HLA_4digit, Length 등은 그대로 복제
    out$Perfect_paired_unique_read_count <- out$Perfect_paired_unique_read_count / 2
    out
  } else if (nrow(df_locus) == 2) {
    df_locus
  } else {
    stop(sprintf("Locus %s has %d alleles; expected 1 or 2.",
                 df_locus$locus[1], nrow(df_locus)))
  }
}

hla2 <- hla %>%
  group_by(locus) %>%
  group_modify(~expand_homozygous(.x)) %>%
  ungroup()

hla2 <- data.frame(hla2)
hla <- hla2
hla$RPK <- hla$Perfect_paired_unique_read_count / (hla$Length / 1000)

d$RPK <- ifelse(d$effective_length > 0,
                d$expected_count / (d$effective_length / 1000),
                0)

original_all_RPK <- sum(d$RPK)
hla_a_rpk <- d[d$gene_id == hla_a,]$RPK
hla_b_rpk <- d[d$gene_id == hla_b,]$RPK
hla_c_rpk <- d[d$gene_id == hla_c,]$RPK

original_all_RPK_except_HLA <- original_all_RPK - (hla_a_rpk + hla_b_rpk + hla_c_rpk)
hla_allele_all_RPK <- sum(hla$RPK)
new_all_RPK <- original_all_RPK_except_HLA + hla_allele_all_RPK

d$new_TPM <- d$RPK / new_all_RPK * 1000000
hla$new_TPM <- hla$RPK / new_all_RPK * 1000000
tmp <- data.frame(gene_id = c("HLA-A1","HLA-A2","HLA-B1","HLA-B2","HLA-C1","HLA-C2"),new_TPM = hla$new_TPM)
d2 <- d[,c("gene_id","new_TPM")]
df_t <- rbind(d2,tmp)

df_t <- df_t %>% filter(gene_id != hla_a & gene_id != hla_b & gene_id != hla_c)
names(df_t)[2] <- hla_file

}
else{
 table_file <- paste("/BiO/Hyein/OptiType/", hla_file, "/", hla_file, "_table_result.tsv", sep = "")
  table <- read.table(table_file, sep = "\t", header = TRUE)
  table <- table %>%
  filter(!grepl("^[0-9.]+$", Original)) %>%  # Original이 숫자(소수 포함)인 행 제거
  filter(!grepl("^[0-9.]+$", HLA_4digit))  

  ## 여기서 table을 이용해 동형 locus 찾기 & allele ID 붙이기
  ## 예: table 컬럼 Original=HLA000xx, HLA_4digit=A*xx:xx
  hla <- hla %>%
    left_join(table, by = c("Allele" = "Original"))

##Calculate original each RPK

hla$RPK <- hla$Perfect_paired_unique_read_count / (hla$Length / 1000)

d$RPK <- ifelse(d$effective_length > 0,
                d$expected_count / (d$effective_length / 1000),
                0)

original_all_RPK <- sum(d$RPK)
hla_a_rpk <- d[d$gene_id == hla_a,]$RPK
hla_b_rpk <- d[d$gene_id == hla_b,]$RPK
hla_c_rpk <- d[d$gene_id == hla_c,]$RPK

original_all_RPK_except_HLA <- original_all_RPK - (hla_a_rpk + hla_b_rpk + hla_c_rpk)
hla_allele_all_RPK <- sum(hla$RPK)
new_all_RPK <- original_all_RPK_except_HLA + hla_allele_all_RPK

d$new_TPM <- d$RPK / new_all_RPK * 1000000
hla$new_TPM <- hla$RPK / new_all_RPK * 1000000
tmp <- data.frame(gene_id = c("HLA-A1","HLA-A2","HLA-B1","HLA-B2","HLA-C1","HLA-C2"),new_TPM = hla$new_TPM)
d2 <- d[,c("gene_id","new_TPM")]
df_t <- rbind(d2,tmp)

df_t <- df_t %>% filter(gene_id != hla_a & gene_id != hla_b & gene_id != hla_c)
names(df_t)[2] <- hla_file

if( nrow(df) == 0){
df <- df_t
}
else{
df <- left_join(df,df_t,by="gene_id")
}
}
hla$locus <- c("HLA-A1","HLA-A2","HLA-B1","HLA-B2","HLA-C1","HLA-C2")

write.table(hla,paste("/BiO/Hyein/OptiType/",hla_file,"/",hla_file,"_all_table.txt",sep=""),sep="\t",quote=F,row.names=F)

}






##estimate molecules
dd2 <- data.frame()

 for ( i in 2:(ncol(dd)-1)){
 tmp <- dd[,c(1,i)]
 mol <- (tmp[[2]]/1000000) * 300000
tmp[[2]] <- mol
 if ( i ==2){

dd2 <- tmp
 }
 else{
 dd2 <- left_join(dd2,tmp,by="gene_id")
 }
 }


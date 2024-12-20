library(DrugTargetMR)
library(tidyverse)
devtools::install_version("parallel",version ="4.2.0")
BiocManager::install("BiocParallel")


#____________________________________________________________________________________________# load and transform to standard data
prepare_finngen(
  file_path = "finngen_R12_D3_SARCOIDOSIS.gz",
  out_path = "./",
  generate_mr = T,
  generate_smr = F,
  sample_size = 497722
)
a = data.table::fread("finngen_R12_D3_SARCOIDOSIS.txt",nrows = 3)


#transform_hg38ToHg19(file_path = "finngen_R12_D3_SARCOIDOSIS.txt")
files <- dir(pattern = ".gz", full.names = TRUE)[!dir(pattern = ".gz", full.names = TRUE) %in% "./finngen_R12_D3_SARCOIDOSIS.gz"]
prepare_finngen(
  file_path = c(files),
  out_path = "./",
  generate_mr = T,
  generate_smr = F,
  sample_size = 497722
)
#____________________________________________________________________________________________# MR analysis

mr_clump(
  local_exp_path =c("./finngen_R12_D3_SARCOIDOSIS.txt") ,
  exp_p = 5e-8,
  maf = 0.01,
  Fvalue = 10,
  cis =F,
  #cis_trans_col = "cis_trans_hg19",
  clump =T,
  clump_kb = 10000,
  clump_r2 = 0.001,
  bfile = "d:/GWAS/data_ref/1kg.v3/EUR",
  plink_exe = get_plink_exe(),
  out_path = "./MR")

files <- dir(pattern = ".txt", full.names = TRUE)[dir(pattern = ".txt", full.names = TRUE) != "./finngen_R12_D3_SARCOIDOSIS.txt"]
mr_local_clumped(
  local_exp_path = "./MR/mr_clump_method1-5e-08-clumped-10000kb-0.001r2.csv",
  out_p = 0,
  local_out_path=c(files),
  bfile = "d:/GWAS/data_ref/1kg.v3/EUR",
  out_path = "./MR")


#____________________________________________________________________________________________# LDSC analysis

files <- dir(pattern = ".txt", full.names = TRUE)[dir(pattern = ".txt", full.names = TRUE) != "./finngen_R12_D3_SARCOIDOSIS.txt"]
for (file in files) {
  postgwas_ldsc_rg(file_path =  c("finngen_R12_D3_SARCOIDOSIS.txt",file),out_prefix = paste0("Sarcoidosis"," with ",sub(".*/(.*?)\\.txt", "\\1", file)),
                  out_path = "./LDSC/",ancestry = "EUR")}



#____________________________________________________________________________________________# SCCA跨组织关联分析

## 提前导入用于fusion分析的sumstat_data数据，类型为data.frame，可减少循环过程中的运算量
sumstat_data <- data.table::fread("./Sarcoidosis.sumstats",data.table=F)

## 使用循环，挨个组织执行FUSION分析
files_pos<- dir(path = "d:/GWAS/data_fusion/weight_file/sCCA构建的GTEx跨组织基因权重/sCCA_weights_v8_2/sCCA_weights_v8/",".pos",full.names = T)

df_weights <- data.frame(files = files_pos,
                         weight_dir = dirname(files_pos),
                         sCCA = basename(files_pos))


for (i in 1:nrow(df_weights)) {
  print(i)
  fusion_assoc(sumstat_data = sumstat_data,
               weights = df_weights$files[i],
               weights_dir =df_weights$weight_dir[i],
               ref_ld_chr_prefix = "d:/GWAS/data_fusion/LDREF/EUR/1000G.EUR.",
               ref_ld_chr_num = 1:22,
               out_path ="./results/results_fusion_scca/",
               out_prefix = df_weights$sCCA[i])
}

# acat-o汇总p值
fusion_acatP(file_path = dir("./results/results_fusion_scca/",pattern = ".csv",full.names = T),
             out_path = "./results/results_fusion_scca/",
             out_prefix = "scca")
actP = read.csv("results/results_fusion_scca/scca_compute_acatP.csv")
actP$P.adj = p.adjust(actP$acatP,method = "fdr")
actP$gene = sub("\\..*", "", actP$ID)
actP = actP %>% left_join(Ensembl_GRCh38[,c("gene_id","gene_name")],by = c("gene" = "gene_id"))
actP$Gene.name = actP$gene_name
SCCA_gene = actP %>% filter(P.adj<0.05) %>% pull(Gene.name)
table(actP$P.adj<0.05)
write.csv(actP,file = "results/results_fusion_scca/actP.csv")
write.csv(SCCA_gene,file = "results/results_fusion_scca/SCCA.fdr.gene.csv")

SCCA_gene = read.csv("results/results_fusion_scca/SCCA.fdr.gene.csv")%>% pull(x)

#____________________________________________________________________________________________________# UTMOST
utmost <- utmost_single_tissue(gwas_file = "./finngen_R12_D3_SARCOIDOSIS.txt",
                               cov_folder  = "d:/GWAS/data_utmost/GTEX_v8_20230609/covariance_GTEx8_normalized_pruned-20240619T071616Z-001/covariance_GTEx8_normalized_pruned/",
                               weight_db_folder  = "d:/GWAS/data_utmost/GTEX_v8_20230609/database_normalized_pruned-20240619T072828Z-001/database_normalized_pruned/",
                               # single_tissue = c("Adipose_Subcutaneous","Artery_Tibial"),  # 默认使用指定文件夹下的所有.db文件
                               out_path = "./results/results_utmost",workers = 25,
                               out_prefix = "utmost")



utmost_joint_GBJ(gene_info_file = "./results/results_utmost/utmost_summary_gene_info.txt",
                 input_file = "./results/results_utmost/utmost_summary_results.csv",
                 
                 # cov_rds，设置预先加载的协方差数据，covariance_joint_GTEx8_normalized_pruned.rds
                 cov_rds = "d:/GWAS/data_utmost/GTEX_v8_20230609/covariance_joint_GTEx8_normalized_pruned.rds",

                 # weight_db_rds，设置预先加载的权重数据，database_normalized_pruned.rds
                 weight_db_rds = "d:/GWAS/data_utmost/GTEX_v8_20230609/database_normalized_pruned.rds",

                 # 设置并行线程
                 workers = 1, 
                 
                 # 设置输出
                 out_path = "./results/results_utmost",
                 out_prefix ="utmost_sarcoidosis")


abc = read.csv("results/results_utmost/joint_GBJ_joint_GBJ.csv")
abc$p.adj = p.adjust(abc$p_value,method = "fdr")
abc$gene = sub("\\..*", "", abc$gene)
abc = abc %>% left_join(Ensembl_GRCh38[,c("gene_id","gene_name")],by = c("gene" = "gene_id"))
abc$Gene.name = abc$gene_name
utmost_gene = abc %>% filter(p.adj<0.05) %>% pull(Gene.name)
write.csv(abc,file = "./results/results_utmost/utmost.summary.result.fdr.csv")
write.csv(utmost_gene,file = "./results/results_utmost/utmost.fdr.gene.csv")






#____________________________________________________________________________________________________# FUSION
dir_tar <- "d:/GWAS/data_fusion/weight_file/Mancuso-GTEx/GTExv8.ALL2/"
dir_tar_file <- dir(dir_tar,full.names = T)

dir_untar <- "d:/GWAS/data_fusion/weight_file/Mancuso-GTEx/GTExv8.ALL2/"
dir_weight <- list.files(dir_untar,recursive = FALSE, full.names = TRUE)
file_weights <-  unlist(lapply(dir_weight, FUN = function(x){ dir(x,pattern = ".pos$",recursive = F,full.names = T)}))
file_weights <- file_weights[!grepl(".nofilter.pos",file_weights)]
df_weights <- data.frame(file_weights = file_weights,weights_dir = dirname(file_weights))
df_weights$tissues <- basename(df_weights$file_weights)
sumstat_data <- data.table::fread("./Sarcoidosis.sumstats",data.table=F)

## 使用循环，挨个组织执行FUSION分析
for (i in 1:nrow(df_weights)) {
  fusion_assoc(sumstat_data = sumstat_data,
               weights = df_weights$file_weights[i],
               weights_dir =df_weights$weights_dir[i],
               ref_ld_chr_prefix = "d:/GWAS/data_fusion/LDREF/EUR/1000G.EUR.",
               ref_ld_chr_num = 1:22,
               out_path ="./results/results_fusion/",
               out_prefix = df_weights$tissues[i])
}

res_fusion_files <- dir("./results/results_fusion/",pattern = "summary.csv",recursive = F,full.names = T)
res_fusion <- data.frame()
for (i in res_fusion_files) {
  
  if(!file.exists(i)){
    next
  }
  
  temp_res <- read.csv(i)
  temp_res <- temp_res[!is.na(temp_res$TWAS.P),]
  temp_res$TWAS.P.fdr <- p.adjust(temp_res$TWAS.P,method = "fdr")
  
  if(nrow(temp_res)>0){
    res_fusion <- rbind(res_fusion,temp_res)
  }
}
res_fusion$gene = sub("\\..*", "", res_fusion$ID)
res_fusion = res_fusion %>% left_join(Ensembl_GRCh38[,c("gene_id","gene_biotype","gene_name")],by = c("gene" = "gene_id"))
Fusion_gene <- res_fusion[res_fusion$TWAS.P.fdr < 0.05,] %>% pull(gene_name)
write.csv(res_fusion,file = "./results/results_fusion/Fusion.p_fdr.csv")
write.csv(Fusion_gene,file = "./results/results_fusion/Fusion_gene.csv")





#____________________________________________________________________________________________________# 基因取交集

utmost_gene %in% SCCA_gene
SCCA_gene %in% utmost_gene
table(utmost_gene %in% SCCA_gene)
utmost_gene[utmost_gene%in%SCCA_gene][utmost_gene%in%Fusion_gene]

# "RNF215"      "SF3B1"       "PLCL1"       "FAM117B"     "RFTN2"       "MIR4435-2HG"

res_fusion = read.csv("./results/results_fusion/Fusion.p_fdr.csv")
res_fusion = res_fusion %>% filter(TWAS.P.fdr < 0.05) %>% filter(gene_name %in% c("RFTN2","PLCL1","CCDC157","WDR12","MIR4435-2HG"))
write.csv(res_fusion,file = "./results/results_fusion/res_fusion.filter.csv")


#____________________________________________________________________________________________________# COJO分析
fusion = read.csv("results/results_fusion/Fusion.p_fdr.csv")
file_genesymbol = fusion %>% filter(gene_name %in% c("RFTN2","PLCL1","CCDC157","WDR12","MIR4435-2HG"))
file_genesymbol$ID = file_genesymbol$gene_name
write.csv(file_genesymbol,file = "results/cojo.input.csv")
fusion_cojo(input = "results/cojo.input.csv",
            out_path = "./results/results_fusion/cojo",
            sumstats = "./Sarcoidosis.sumstats",
            ref_ld_chr_prefix = "d:/GWAS/data_fusion/LDREF/EUR/1000G.EUR.",
            ref_ld_chr_num =  c(2,22),
            locus_win = 5e+05,
            p_adj_method = "none",
            plot = T,
            plot_legend = "joint",
            glist = "d:/GWAS/data_fusion/glist-hg19")







#____________________________________________________________________________________________________# MR分析

### 基于smr格式（.besd；.epi; .esi）的文件，查询用于mr分析的暴露输入数据
### 第一步：获取可用组织的文件名前缀
qtls_path <- "d:/GWAS/GTEx_V8_cis_eqtl_summary_lite/"
qtls_path_prefix <- tools::file_path_sans_ext(dir(qtls_path,pattern = ".epi",full.names = T))

### 第二步：准备用于查询的探针：对应.epi文件的v2列。
probe <- Ensembl_GRCh38[Ensembl_GRCh38$gene_name %in% c("RNF215","SF3B1","PLCL1","FAM117B","RFTN2","MIR4435-2HG"),"gene_id"]


### 第三步：批量查询
?smr_query
for (temp_qtls_path in qtls_path_prefix) {
  smr_query(smr_exe_path = "d:/GWAS/SMR/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe",
            qtls_path = temp_qtls_path,
            query_probe = probe,
            query_p = 1,
            cis_wind = 2000,
            #combine = T, # 更新R包之后，增加combine参数，可将查询到结果进行合并输出。作用，查询多个探针效率更高
            out_path = "./results/result_mr/data_prepare",
            out_prefix = gsub(".lite","",basename(temp_qtls_path)))
}


### 第四步：鉴于mr_local进行合并得到时候，会将相同Phenotype对应的行进行合并分析，故执行Phenotype列去除
query_files <- dir("results/result_mr/data_prepare/",recursive = T,full.names = T)
print(query_files)

238+98+20
14+19+8
138+8


mr_clump(
  local_exp_path =c(query_files),
  exp_p = 5e-5,
  maf = 0.01,
  Fvalue = 10,
  cis =F,
 # cis_trans_col = "cis_trans_hg38",
  clump =T,
  clump_kb = 10000,
  clump_r2 = 0.001,
  bfile = "d:/GWAS/data_ref/1kg.v3/EUR",
  plink_exe = get_plink_exe(),
  out_path = "./results/result_mr/")

mr_local_clumped(
  local_exp_path = "./results/result_mr/mr_clump_method1-5e-05-clumped-10000kb-0.001r2.csv",
  out_p = 0,
  local_out_path=c('./finngen_R12_D3_SARCOIDOSIS.txt'),
  bfile = "d:/GWAS/data_ref/1kg.v3/EUR",
  out_path = "./results/result_mr/")


prepare_finngen(file_path ="finngen_R12_D3_SARCOIDOSIS.gz",generate_mr = F,generate_smr = T)
smr_prepare_ma(file_path ="finngen_R12_D3_SARCOIDOSIS.ma",out_path = "./")                

qtl_besd <- dir("d:/GWAS/GTEx_V8_cis_eqtl_summary_lite/",pattern = ".besd$",full.names = T)
qtl_besd <- qtl_besd[grepl(paste(c("Thyroid", "Spleen", "Ovary", "Heart_Left", "Breast"), collapse = "|"), qtl_besd)]


qtl_besd <- gsub(".besd","",qtl_besd)
for (i in qtl_besd) {
  print(i)
  smr_qtl2gwas(smr_exe_path = "d:/GWAS/SMR/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe",
               bfile = "d:/GWAS/data_ref/1kg.v3/EUR",
               gwas_path = "finngen_R12_D3_SARCOIDOSIS_treated.ma",
               qtls_path = i,
               smr_multi_snp = T , # 基于多个SNP计算p值。
               smr_peqtl = 5e-08,
               smr_cis_wind = 1000,  # cis_wind值得是探针上游或者下游的碱基数量，即单侧碱基数量，单位为kb
               
               out_path = "./results/results_smr/",
               out_prefix = basename(i))
}


#____________________________________________________________________________________________________# coloc分析

##########################  如果是跑全部组织的共定位分析，可以使用如下代码：
# 拿到GTEx_v8，含有49个组织数据之后，可以使用如下代码进行批量解
# 1. 获取GTEx文件夹下的所有的权重文件，为.tar.gz的形式
dir_gz <- "d:/GWAS/GTEXcis/"
dir_gz_file <- dir(dir_gz,full.names = T)

# 2. 创建一个新的文件夹，用于存放加压后的文件
dir_unigz <- "d:/GWAS/GTEx_V8_cis_eqtl_summary_hg19_48G"
dir.create(dir_unigz)

# 3. 使用lapply，批量进行tar包解压
lapply(dir_gz_file, FUN = function(x){unzip(zipfile  =x,exdir = dir_unigz)})

library(openxlsx)
gtex_ss <- read.xlsx("d:/GWAS/GTEx_v8_science_2020 附图和附表-附表2.xlsx")
gtex_folder <- dir("d:/GWAS/GTEx_V8_cis_eqtl_summary_hg19_48G",full.names = T)
gtex_list <- data.frame(folder = gtex_folder,
                        Tissue = basename(gtex_folder))
gtex_list <- merge(gtex_list,gtex_ss[,c("Tissue","Samples")],by = "Tissue")
probe <- Ensembl_GRCh38[Ensembl_GRCh38$gene_name %in%c( "RNF215","SF3B1","PLCL1","FAM117B","RFTN2","MIR4435-2HG"),"gene_id"]


### 使用for循环，完成批量共定位
for (i in 1:nrow(gtex_list)) {
  print(i)
  # 结局GWAS样本量信息，20908 case and 312803 controls
  coloc_batch_smr(probes = probe,
                  smr_exe_path = "d:/GWAS/SMR/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe",
                  qtl_path = gtex_list$folder[i], # 这里传对应文件夹的路径即可
                  cis_wind = 1000,
                  
                  # prepare_method=2,按照SNP的rsid选择SNP，以GWAS1的染色体坐标为参考，不使用GWAS2的chr和pos信息。
                  prepare_method = 2,
                  type1 = "quant",
                  SS1 = gtex_list$Samples[i],
                  gwas_path2 = "./finngen_R12_D3_SARCOIDOSIS.txt",
                  type2 = "cc",
                  SS2 = 20908+312803, # 对应sample.size
                  NC2 = 20908,        # 对应ncase
                  bfile = "d:/GWAS/data_ref/1kg.v3/EUR",
                  coloc_plot = T,
                  coloc_plot_pph4 = 0.75,
                  coloc_plot_genome = "hg38",
                  out_path = paste0("./results/results_coloc/",gtex_list$Tissue[i]))
}


#____________________________________________________________________________________________# 绘制森林图
mr = read.csv("results/result_mr/table_s2_mr_results.csv")
mr$id.exposure <- sub("_ENSG.*", "", mr$id.exposure)
mr = mr %>% filter(pval_mr<0.05)
RFTN2 <- mr %>% filter(str_detect(exposure, "ENSG00000162944")) %>% 
  filter(id.exposure %in% c("Breast_Mammary_Tissue","Ovary","Spleen","Thyroid")) %>%
  mutate(Gene_symbol = "RFTN2")
FAM117B <- mr %>% filter(str_detect(exposure, "ENSG00000138439")) %>% 
  filter(id.exposure %in% c("Whole_Blood","Pancreas","Muscle_Skeletal","Spleen","Adipose_Visceral_Omentum","Adipose_Subcutaneous")) %>%
  mutate(Gene_symbol = "FAM117B")
RNF215 <- mr %>% filter(str_detect(exposure, "ENSG00000138439")) %>% 
  filter(id.exposure %in% c("Testis","Skin_Sun_Exposed_Lower_leg","Skin_Not_Sun_Exposed_Suprapubic","Brain_Caudate_basal_ganglia","Artery_Tibial")) %>%
  mutate(Gene_symbol = "RNF215")
mr = rbind(RFTN2,FAM117B,RNF215)

library(forestplot)
mr$pval_mr <- format(mr$pval_mr, scientific = TRUE, digits = 4)
p1 <- visualization_forest(data_input =mr,
                           col_labeltext = c("id.exposure","Gene_symbol","pval_mr"),
                           col_labeltext_rename =  c("Tissue","Gene_symbol","P value"),
                           xlog = F,
                           graph_pos = 3,
                           ci_pos = 3,
                           pdf_out = F,
                           height = 10,
                           ci_pos_sep = ",",
                           width = 8)
p2 <- p1 |> fp_add_lines(h_6 = gpar(lty = 2),
                         h_11 = gpar(lty = 2))
p2












#____________________________________________________________________________________________# Venn

utmost_gene = read.csv("d:/sarcoidosis_lung/results/results_utmost/utmost.fdr.gene.csv") %>%
  pull(x) %>%
  .[ . != "" ] %>%
  na.omit()
length(utmost_gene)
SCCA_gene = read.csv("d:/sarcoidosis_lung/results/results_fusion_scca/SCCA.fdr.gene.csv")%>% pull(x) %>%
  .[ . !=""]%>%
  na.omit()
length(SCCA_gene)

utmost_gene[utmost_gene %in% SCCA_gene] #

res_fusion = read.csv("d:/sarcoidosis_lung/results/results_fusion/Fusion.p_fdr.csv") %>% filter(TWAS.P.fdr < 0.05)

final_genes <- intersect(SCCA_gene[SCCA_gene %in% utmost_gene], res_fusion$gene_name) |> print()

data2 = res_fusion[res_fusion$gene_name %in% final_genes, c("gene_biotype", "CHR", "P0", "P1", "gene_name")] %>%
  dplyr::distinct(gene_name, .keep_all = TRUE)

data2 = data2 %>%
  mutate(Shape = case_when(gene_biotype == "protein_coding" ~ "box",
                           TRUE ~ "circle")) %>%
  mutate(color = case_when(gene_biotype == "protein_coding" ~ "6a3d9a",
                           TRUE ~ "33a02c")) %>%
  rename(Type = gene_name, Chr = CHR, Start = P0, End = P1) %>% 
  select(-gene_biotype)

write.table(data2,file = "data2.txt",sep = "\t",row.names = F)

#"RNF215","SF3B1","PLCL1","FAM117B","RFTN2","MIR4435-2HG"








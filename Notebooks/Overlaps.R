subgroup_regions_genes %>%
    mutate(MAF = ATAC_DE_MAF_vs_ND == RNA_DE_MAF_vs_ND & RNA_DE_MAF_vs_ND != 0,
                      MMSET = ATAC_DE_MMSET_vs_ND == RNA_DE_MMSET_vs_ND  & RNA_DE_MMSET_vs_ND != 0,
                       CCND1 = ATAC_DE_CCND1_vs_ND == RNA_DE_CCND1_vs_ND & RNA_DE_CCND1_vs_ND != 0,
                       HD = ATAC_DE_HD_vs_ND == RNA_DE_HD_vs_ND   & RNA_DE_HD_vs_ND != 0) %>%
     group_by(gene_id) %>%
     summarise_at(c("MMSET", "MAF", "HD", "CCND1"), any) %>% 
     mutate(count = factor(MMSET+MAF+HD+CCND1)) %>%
     summary()



 subgroup_regions_genes %>%
     mutate(MAF = ATAC_DE_MAF_vs_ND ==1 & RNA_DE_MAF_vs_ND != 0,
                       MMSET = ATAC_DE_MMSET_vs_ND == 1 & RNA_DE_MMSET_vs_ND != 0,
                       CCND1 = ATAC_DE_CCND1_vs_ND == 1 & RNA_DE_CCND1_vs_ND != 0,
                      HD = ATAC_DE_HD_vs_ND ==1 & RNA_DE_HD_vs_ND != 0) %>%
   # select(CCND1, HD, MMSET, MAF) %>%
     group_by(gene_id) %>%
     summarise_at(c("MMSET", "MAF", "HD", "CCND1"), any) %>% 
     mutate(count = factor(MMSET+MAF+HD+CCND1)) %>%
     summary()
 
 subgroup_regions_genes %>%
   mutate(MAF = ATAC_DE_MAF_vs_ND ==1 ,
          MMSET = ATAC_DE_MMSET_vs_ND == 1 ,
          CCND1 = ATAC_DE_CCND1_vs_ND == 1 ,
          HD = ATAC_DE_HD_vs_ND ==1 ) %>%
   # select(CCND1, HD, MMSET, MAF) %>%
   group_by(gene_id) %>%
   summarise_at(c("MMSET", "MAF", "HD", "CCND1"), any) %>% 
   mutate(count = factor(MMSET+MAF+HD+CCND1)) %>%
   summary()
 
 subgroup_regions_genes %>%
   mutate(MAF = ATAC_DE_MAF_vs_ND != 0 & RNA_DE_MAF_vs_ND != 0,
          MMSET = ATAC_DE_MMSET_vs_ND != 0 & RNA_DE_MMSET_vs_ND != 0,
          CCND1 = ATAC_DE_CCND1_vs_ND != 0 & RNA_DE_CCND1_vs_ND != 0,
          HD = ATAC_DE_HD_vs_ND !=0 & RNA_DE_HD_vs_ND != 0) %>%
   # select(CCND1, HD, MMSET, MAF) %>%
   group_by(gene_id) %>%
   summarise_at(c("MMSET", "MAF", "HD", "CCND1"), any) %>% 
   mutate(count = factor(MMSET+MAF+HD+CCND1)) %>%
   summary()
 
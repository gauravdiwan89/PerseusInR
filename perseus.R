library(tidyverse)
source("functions.R")

data <- read_tsv("report.pg_matrix.tsv")

comparisons_of_interest <- c("Arl13b_TurboID_vs_mTurbo", "Arl13b_TID_Ciliobrevin_vs_mTurbo")

test <- data %>% 
  select(`Protein IDs` = Protein.Group, `Gene names` = Genes, contains(".raw")) %>% 
  pivot_longer(cols = c(-`Protein IDs`, -`Gene names`), names_to = "mut") %>% 
  mutate(mut = sub(".*GD\\d+_", "", str_extract(mut, "([^\\\\]+)(?=\\.raw$)"))) # %>% 
# drop_na(value)

stat_test <- test %>% 
  group_by(`Protein IDs`) %>% 
  do({
    sub_tbl <- .
    strains <- unique(sub_tbl$mut) #figure out which mutations are in the proteingroups file
    expand.grid(g1 = strains, g2 = strains, stringsAsFactors = FALSE) %>% #make a combination of WT vs MUT
      filter(g1 != g2) %>% #filter out mutant vs mutant
      mutate(test = as.character(glue::glue("{g2}_vs_{g1}"))) %>% #create a column that tells us which two strains are being compared, like WT_vs_K0
      filter(test %in% comparisons_of_interest) %>% 
      rowwise() %>% #do the next steps for each combination
      do({
        gs <- .
        
        test_x <- sub_tbl %>% filter(mut == gs$g1) %>% drop_na(value) %>% pull(value) #take the intensity values for WT
        log_x <- log2(test_x) #do a log2
        test_y <- sub_tbl %>% filter(mut == gs$g2) %>% drop_na(value) %>% pull(value) #take the intensity values for the mutants from the same timepoint
        log_y <- log2(test_y) #do a log2
        #the zero intensity values will be converted to infinity in the earlier step
        if((sum(!is.infinite(log_x)) > 3 | sum(!is.infinite(log_y)) > 3) & (length(test_x) > 2 & length(test_y) > 2)) { 
          new_x <- sapply(test_x, function(y) if(y == 0) impute_norm(x = test_x, down_shift = 1.8, width = 0.5) else y) #here we impute missing values with values from a normal distribution. You could change the values for down shift and width here as you do in Perseus
          new_y <- sapply(test_y, function(y) if(y == 0) impute_norm(x = test_y, down_shift = 1.8, width = 0.5) else y) #here we impute missing values with values from a normal distribution. You could change the values for down shift and width here as you do in Perseus
          t_res <- t.test(log2(new_x), log2(new_y)) #do a t-test
          ratio <- log2(median(new_y)/median(new_x)) #find the ratio of the medians of others vs K0
          data.frame(test=gs$test, estimate=t_res$statistic, p_value=t_res$p.value, ratio = ratio, stringsAsFactors = FALSE) #export the t-test output with the ratio
        } else {
          data.frame(test=gs$test, estimate=NA, p_value=NA, ratio = NA,  stringsAsFactors = FALSE) #empty output if there were more than 2 zero values in both the sets
        }
      })
  }) %>%
  ungroup() %>%
  #clean up the statistics output table
  pivot_longer(cols = estimate:ratio, names_to ="statistic", values_to="val") %>% 
  unite(test, statistic, col = "test_statistic", sep = "_", remove = T) %>%
  pivot_wider(names_from = "test_statistic", values_from = "val")

##calculate the corrected p values and significance A values
max_ratio <- stat_test %>% select(ends_with("_ratio")) %>% drop_na() %>% pivot_longer(cols = everything()) %>% filter(!is.infinite(value)) %>% filter(value==max(value)) %>% pull(value)+1 ##calculate the max ratio to replace +Inf values with

min_ratio <- stat_test %>% select(ends_with("_ratio")) %>% drop_na() %>% pivot_longer(cols = everything()) %>% filter(!is.infinite(value)) %>% filter(value==min(value)) %>% pull(value)-1 ##calculate the min ratio to replace -Inf values with

post_test <- stat_test %>% 
  # group_by(time) %>% 
  # filter_if(is.double, any_vars(!is.na(.))) %>% #remove genes with no tests in all mutants
  mutate(across(contains("_p_value"), .fns = p.adjust, method = "BH", .names = "{gsub('_p_value', '', col)}_q_value")) %>% #do a Benjamini-Hochberg correction on the pvalue columns
  # mutate(across(contains("_ratio"), .fns = (-1*.))) %>%
  mutate(across(ends_with("_ratio"), .fns = significanceA, .names = "{gsub('_ratio', '', col)}_sigA")) %>% 
  mutate(across(ends_with("_ratio"), .fns = function(x) ifelse(is.infinite(x), ifelse(x<0, min_ratio, max_ratio), x)))

comparison_table <- post_test %>% 
  select(`Protein IDs`, contains("ratio"), contains("_q_value"), contains("_sigA")) %>% 
  # filter(across(.cols = -`Protein IDs`, .fns = ~ any(.x < 0.05))) %>% 
  pivot_longer(cols = -c(`Protein IDs`)) %>% 
  extract(col = name, into = c("comparison", "statistic"), regex = "(^.*_vs_.*)_(ratio|q_value|sigA$)") %>% 
  pivot_wider(names_from = statistic, values_from = value) %>% 
  mutate(Gene = data$Genes[match(`Protein IDs`, data$Protein.Group)], .after = `Protein IDs`) %>% 
  mutate(sign = case_when(
    q_value < 0.05 & sigA < 0.05 ~ "Tier 1",
    sigA < 0.05 ~ "Tier 2",
    q_value < 0.05 ~ "Tier 3",
    .default = "NS"
  ))

table_oi <- comparison_table

#volcano plots
mean_ratio <- mean(table_oi$ratio[!(table_oi$ratio == max_ratio | table_oi$ratio == min_ratio)], na.rm = T)
two_sd <- 2*sd(table_oi$ratio[!(table_oi$ratio == max_ratio | table_oi$ratio == min_ratio)], na.rm = T)

genes_to_label <- c("Kit
Ehd3
Ehd1
Snap29
Kif14
Snap25
Nphp4
Cep131
Bub3
Snap25
Kif11
Arl13b
Sstr3
Kif21a
Kif1a
IFT81
IFT80
IFT172
BBS2
IFT46
IFT140
INPP5E")
genes_to_label <- str_to_sentence(unique(unlist(strsplit(genes_to_label, "\n")))) 

volcano_plot <- table_oi %>% 
  filter(comparison == "Arl13b_TurboID_vs_mTurbo") %>% 
  mutate(across(.cols = q_value:sigA, .fns = ~ replace_na(.x, 1))) %>% 
  mutate(significance = if_else(q_value < 0.05 & sigA < 0.05, "Both positive", if_else(q_value < 0.05, "T-test significant", if_else(sigA < 0.05, "Significance A positive", "NS")))) %>% 
  mutate(label = if_else(Gene %in% genes_to_label, Gene, "")) %>% 
  ggplot(aes(x = ratio, y = -log10(q_value))) +
  geom_hline(yintercept = -log10(0.05), lty = 2) +
  geom_vline(xintercept = c(0), lty = 2) +
  geom_point(aes(fill = sign), colour = "gray70", size = 3, pch = 21) +  
  ggrepel::geom_label_repel(aes(label = label, fill = sign), size = 5, alpha = 0.7, show.legend = F, max.overlaps = 100000, min.segment.length = 0, force = 5) +
  scale_fill_manual(values = c("white", "#f05039", "#e57a77", "#eebab4")) +
  # scale_colour_manual(values = c("gray80", "#f05039", "#e57a77", "#eebab4")) +
  
  # geom_text(aes(label = label), size = 5, alpha = 0.7, show.legend = F) +
  facet_wrap(~comparison) +
  xlim(c(min_ratio, max_ratio)) +
  cowplot::theme_minimal_grid(font_size = 16) +
  # scale_color_brewer(palette = "Set1") +
  labs(x = "log2 (ratio)",
       y = "-log10 (p-values)") +
  theme(panel.spacing = unit(2, "lines"))

volcano_plot

cowplot::save_plot(filename = "volcano_plot_full.pdf", volcano_plot, base_width = 18, base_height =6)
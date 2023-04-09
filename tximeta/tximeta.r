

# Load packages
library(tidyverse)
library(BiocManager)
library(tximeta)
library(biomaRt)
library(plotly)

## Tximeta
# Load sample table
sampleData <- as.data.frame(read_tsv("/Users/keldaperumal/Desktop/MSc_stuff/data/sample_metadata.tsv"))

# Check if run_accession is supplied instead of sample_id
if ("run_accession" %in% colnames(sampleData)) {
  sampleData = sampleData %>% dplyr::rename(sample_id = run_accession)
}

# Get salmon quant file paths 
my.files <-  list.files(list.dirs(path = "/Users/keldaperumal/Desktop/MSc_stuff/work/decoy_quants/",
                                  full.names = TRUE, recursive = FALSE),
                        pattern = "quant.sf", full.names = TRUE)

# Append file paths to sampleData
sampleData$files = my.files

# Check if sample_title is supplied
if ("sample_title" %in% colnames(sampleData))
{
  sampleData = sampleData %>% dplyr::rename(names = sample_title)
} else{
  sampleData$names = sampleData$sample_id
}

# Run tximeta
se <- tximeta(sampleData)
saveRDS(se, "work/tximeta/se_decoy.rds")

# Summarize to gene level counts
gse <- summarizeToGene(se)
saveRDS(gse, "work/tximeta/gse_decoy.rds")

abundance <- as.data.frame(gse@assays@data@listData[["abundance"]])
saveRDS(abundance, "work/tximeta/tpm_decoy.rds")

counts <- as.data.frame(gse@assays@data@listData[["counts"]])
saveRDS(counts, "work/tximeta/counts_decoy.rds")

length <- as.data.frame(gse@assays@data@listData[["length"]])
saveRDS(length, "work/tximeta/length_decoy.rds")


##### PLOT TPM #####

# change data frame so each tpm value is associated with a sample name (i.e 2 columns) for plotting
abundance_longer <- abundance %>% pivot_longer(cols = c(1:9), names_to = 'sample', values_to = 'tpm')

# add column with sample descriptions
abundance_longer_with_condition <- abundance_longer %>%
  mutate(condition = case_when(
    abundance_longer$sample=="SRR6997095" ~ "UNTREATED",
    abundance_longer$sample=="SRR6997096" ~ "UNTREATED",
    abundance_longer$sample=="SRR6997123" ~ "UNTREATED",
    abundance_longer$sample=="SRR6997118" ~ "PMA-TREATED",
    abundance_longer$sample=="SRR6997119" ~ "PMA-TREATED",
    abundance_longer$sample=="SRR6997121" ~ "PMA-TREATED",
    abundance_longer$sample=="SRR6997122" ~ "VD3-TREATED",
    abundance_longer$sample=="SRR6997124" ~ "VD3-TREATED",
    abundance_longer$sample=="SRR6997125" ~ "VD3-TREATED",
  ))

# add column with replicate numbers
abundance_longer_with_condition_replicate <- abundance_longer_with_condition %>%
  mutate(replicate = case_when(
    abundance_longer$sample=="SRR6997095" ~ "Replicate 1",
    abundance_longer$sample=="SRR6997096" ~ "Replicate 2",
    abundance_longer$sample=="SRR6997123" ~ "Replicate 3",
    abundance_longer$sample=="SRR6997118" ~ "Replicate 1",
    abundance_longer$sample=="SRR6997119" ~ "Replicate 2",
    abundance_longer$sample=="SRR6997121" ~ "Replicate 3",
    abundance_longer$sample=="SRR6997122" ~ "Replicate 1",
    abundance_longer$sample=="SRR6997124" ~ "Replicate 2",
    abundance_longer$sample=="SRR6997125" ~ "Replicate 3",
  ))


ggplot(abundance_longer_with_condition_replicate, aes(x=condition, y=tpm, colour=replicate)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1,
                                             dodge.width = 0.7), aes(fill = replicate), pch = 16)+ 
  scale_colour_manual(values = c("purple","aquamarine","yellow")) +
  theme_classic() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

  
  
  





















library(tidyverse)
tpm <- read_rds("work/tximeta/tpm_decoy.rds")
tpm_rowname <- rownames_to_column(tpm, var = "gene_id")
ggplot(tpm_rowname, aes(x=gene_id, y=thp1_1))+
  geom_point()



#' ============================================================================
#' 26-Week Fluxomics with M+0 Sum Normalization
#' ============================================================================
#'
#' Manuscript: Multi-omics Reveals How Diet and Exercise Improve Peripheral Neuropathy
#' Authors:    Eid SA, Guo K, Noureldein MH, Jang D-G, Allouch AM, Miller CM, Kiriluk CP,
#   Lentz W, Boutary S, Mule JJ, Pennathur S, Bhatt DK, Feldman EL, Hur J
#'
#' Contact:
#   Dr. Kai Guo        <guokai8@gmail.com>
#   Dr. Junguk Hur     <junguk.hur@med.UND.edu>
#'
#' Description:
#'   26-week isotopologue analysis with M+0 normalization after summing isotopologue
#   values. First sums isotopologue intensities per animal/group, then normalizes
#   by M+0 (unlabeled) fraction. Generates heatmaps excluding M+0 values.
#'
#' Input:
#'   EX01317.xlsx (26-week isotopologue data)
#'
#' Output:
#'   M+0-sum-normalized results: t-test CSVs, heatmaps (excluding M+0),
#   visualization PDFs, fluxomics_26wk_final_m0_avg.csv
#'
#' Related figures:
#'   Fig. 6, Fig. S12
#'
#' ============================================================================
## -----------------------------------------------------------------------------
## Setup
## -----------------------------------------------------------------------------
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  pth <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) NULL)
  if (!is.null(pth) && nzchar(pth)) setwd(dirname(pth))
}

library(ggplot2)
library(ggsignif)
library(ggpubr)
library(ggprism)
library(ggtext)
library(here)
# knitr::knit_hooks$set(webgl = hook_webgl)
library(magick)
library(readxl)
library(lubridate)
library(here)
require(ggpubr)
require(rstatix)
library(tidyverse)
Glucose_dat<-
  read_excel(file.path(getwd(), "EX01317.xlsx"),
             range="EX01317!B8:AH203") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>% filter(!grepl("Pool",Group))%>%
  mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>%
  relocate(ID, Diet, Glucose, Group)
####
Serine_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                       range="EX01317!B224:AH419") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>%   filter(!grepl("Pool",Group))%>% mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>%relocate(ID, Diet, Glucose, Group)
###
Fructose_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                         range="EX01317!B440:AH635") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>% filter(!grepl("Pool",Group))%>% mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>%  relocate(ID, Diet, Glucose, Group)
####
Glutamic_acid_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                              range="EX01317!B656:AH851") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>% filter(!grepl("Pool",Group))%>%  mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>% relocate(ID, Diet, Glucose, Group)
###
Aspartic_acid_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                              range="EX01317!B872:AH1067") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>% filter(!grepl("Pool",Group))%>%   mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>%relocate(ID, Diet, Glucose, Group)
####
Lactic_acid_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                            range="EX01317!B1086:AH1281") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>% filter(!grepl("Pool",Group))%>%   mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>%relocate(ID, Diet, Glucose, Group)
####
R5P_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                    range="EX01317!B1302:AH1497") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>% filter(!grepl("Pool",Group))%>%  mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>% relocate(ID, Diet, Glucose, Group)
###Erythrose 4-phosphate
E4P_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                    range="EX01317!B1511:AH1706") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>% filter(!grepl("Pool",Group))%>%  mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>% relocate(ID, Diet, Glucose, Group)
#####
F6P_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                    range="EX01317!B1724:AH1919") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>% filter(!grepl("Pool",Group))%>%  mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>% relocate(ID, Diet, Glucose, Group)
#####Sedoheptulose 7-phosphate
S7P_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                    range="EX01317!B1939:AH2134") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>% filter(!grepl("Pool",Group))%>%  mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>% relocate(ID, Diet, Glucose, Group)
####Glutathione
Glutathione_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                            range="EX01317!B2151:AH2346") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>% filter(!grepl("Pool",Group))%>%  mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>% relocate(ID, Diet, Glucose, Group)
###sn-Glycero-3-phosphate
sn_Glycero_3_phosphate_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                                       range="EX01317!B2366:AH2561") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>% filter(!grepl("Pool",Group))%>%   mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>%relocate(ID, Diet, Glucose, Group)
####
NAD_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                        range="EX01317!B2574:AH2768") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>%filter(!grepl("Pool",Group))%>%   mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>% relocate(ID, Diet, Glucose, Group)
###
Pyruvic_acid_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                    range="EX01317!B2786:AH2981") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>% filter(!grepl("Pool",Group))%>% mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>%  relocate(ID, Diet, Glucose, Group)
###Oxidized glutathione
Phosphoserine_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                                     range="EX01317!B2999:AH3194") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>% filter(!grepl("Pool",Group))%>% mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>%  relocate(ID, Diet, Glucose, Group)
####
AMP_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                         range="EX01317!B3209:AH3404") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>% filter(!grepl("Pool",Group))%>%  mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>% relocate(ID, Diet, Glucose, Group)
###
Oxidized_glutathione_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                           range="EX01317!B3425:AH3620") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>% filter(!grepl("Pool",Group))%>%  mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>% relocate(ID, Diet, Glucose, Group)
####
Succinic_acid_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                                 range="EX01317!B3638:AH3833") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>% filter(!grepl("Pool",Group))%>%  mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>% relocate(ID, Diet, Glucose, Group)
####
Malic_acid_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                                   range="EX01317!B3854:AH4049") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>% filter(!grepl("Pool",Group))%>%   mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>%relocate(ID, Diet, Glucose, Group)
###
Oxoglutaric_acid_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                                     range="EX01317!B4070:AH4265") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>% filter(!grepl("Pool",Group))%>%  mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>% relocate(ID, Diet, Glucose, Group)
###
X6_Phosphogluconic_acid_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                     range="EX01317!B4279:AH4474") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>%  filter(!grepl("Pool",Group))%>% mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>% relocate(ID, Diet, Glucose, Group)
###
Phosphoglyceric_acid_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                    range="EX01317!B4495:AH4690") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>% filter(!grepl("Pool",Group))%>% mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>%  relocate(ID, Diet, Glucose, Group)
###
#########
NADP_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                                          range="EX01317!B4711:AH4906") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>% filter(!grepl("Pool",Group))%>%  mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>% relocate(ID, Diet, Glucose, Group)
#####################
Fumaric_acid_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                                      range="EX01317!B4926:AH5121") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>% filter(!grepl("Pool",Group))%>%  mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>% relocate(ID, Diet, Glucose, Group)
####
NADH_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                    range="EX01317!B5142:AH5336") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>% filter(!grepl("Pool",Group))%>%  mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>% relocate(ID, Diet, Glucose, Group)
####
ADP_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                    range="EX01317!B5356:AH5551") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>% filter(!grepl("Pool",Group))%>% mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>%  relocate(ID, Diet, Glucose, Group)

Fructose_1_6_bisphosphate_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                    range="EX01317!B5570:AH5765") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>% filter(!grepl("Pool",Group))%>%  mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>% relocate(ID, Diet, Glucose, Group)

Citric_Isocitric_acid_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                    range="EX01317!B5785:AH5980") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>% filter(!grepl("Pool",Group))%>%  mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>% relocate(ID, Diet, Glucose, Group)

############
FAD_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                    range="EX01317!B6001:AI6196") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>%  filter(!grepl("Pool",Group))%>% mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>% relocate(ID, Diet, Glucose, Group)

NADPH_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                    range="EX01317!B6217:AH6408") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>% filter(!grepl("Pool",Group))%>% mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>%  relocate(ID, Diet, Glucose, Group)


ATP_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                    range="EX01317!B6429:AH6624") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>% filter(!grepl("Pool",Group))%>%  mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>% relocate(ID, Diet, Glucose, Group)


Acetyl_CoA_dat<-read_excel(file.path(getwd(), "EX01317.xlsx"),
                    range="EX01317!B6643:AH6835") %>% 
  select(!all_of(starts_with("..."))) %>% 
  mutate(ID=sub(".*-S([0-9]*)-([0-9]*).*", "\\2", `File Name`),
         Diet=sub('26wk ','',sub(".*wk ([A-Z]*) .*", "\\1",sub('\\ -.*','',`SampleGroups`))),
         Glucose=if_else(grepl("glc", SampleGroups, ignore.case=T), "Glucose", "Palm"),
         Group=sub(' ','',sub(".*-(.*)", "\\1", `SampleGroups`))) %>% filter(!grepl("Pool",Group))%>%  mutate(s = substr(sub(".*S", "S", `File Name`), 1, 9)) %>% 
  group_by(s) %>% 
  filter(row_number() == ifelse(n() == 1, 1, 2)) %>% 
  ungroup()%>% relocate(ID, Diet, Glucose, Group)


###########
dat_list<-ls()
dat_list<-dat_list[endsWith(dat_list, "_dat")]
dat<-lapply(dat_list, function(x)get(x))
names(dat)<-dat_list
# ###only Glucose
# datg<-lapply(dat, function(x)subset(as.data.frame(x),Glucose=="Glucose"))
# datgs<-lapply(datg, function(x)x[,c("Diet","Group","m+0")])
# datgss<-lapply(datgs, function(x)x%>%mutate(m=ifelse(`m+0`==0,"Yes","No"))%>%group_by(Diet,Group,m)%>%summarise(m0=n()))
# datgss<-lapply(names(datgss), function(x)datgss[[x]]%>%mutate(isotype=x))
# datgss<-do.call(rbind,datgss)
# datgss<-as.data.frame(datgss)
# datgss<-datgss%>%select(Group,m,m0)%>%spread(Group,m0)%>%column_to_rownames(var="m")%>%t()
# rownames(datgss)<-sub('_dat','',rownames(datgss))
# datgss%>%na.omit()%>%pheatmap::pheatmap(display_numbers = T,cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("cyan4","white","darkorange"))(128),border_color = "white",number_format = "%.0f")


###########################################################
####first calculate the sum and then calculate the ratio
###########################################################
datu<-lapply(dat, function(x)x%>%mutate(group=ifelse(Group%in%c("nd","np"),"Nerve",Group))%>%select(ID:Glucose,contains("m+"),group)%>%
               gather(m,val,-ID,-Diet,-Glucose,-group)%>%group_by(ID,Diet,group,Glucose,m)%>%
                summarise(su=sum(val))%>%spread(m,su))
####
cal_per<-function(dat){
  d<-dat%>%mutate(across(starts_with("m"),~.x/`m+0`))
  return(d)
}
#dats<-lapply(dat, function(x)filter(x,Abund>0))
dats<-lapply(datu, function(x)filter(x,`m+0`>0))
datpu<-lapply(dats, function(x)cal_per(x))
#### merge np and nd with average value
#datpu<-lapply(datp, function(x)x%>%mutate(group=ifelse(Group%in%c("nd","np"),"Nerve",Group))%>%select(ID:Glucose,contains("m+"),group)%>%
#                gather(m,val,-ID,-Diet,-Glucose,-group)%>%group_by(ID,Diet,group,Glucose,m)%>%
#                summarise(mu=mean(val))%>%spread(m,mu))

#######
cal_ttest<-function(x,g){
  x<-ungroup(x)
  mid<-ungroup(x)%>%filter(Glucose=="Glucose",Diet%in%g)%>%select(Diet,starts_with("m"))%>%select_if(~ is.character(.) || any(. > 0))%>%
    gather(m,val,-Diet)%>%group_by(Diet,m)%>%
    ###### filter 1% at least in 1/3 samples
    summarise(count=n(),cc=sum(val>0.01))%>%filter(cc/count>0.3)%>%pull(m)%>%unique()
  mid<-setdiff(mid,"m+0")
  if(length(mid)==0){
    return(NULL)
  }
  res<-x%>%filter(Glucose=="Glucose",Diet%in%g)%>%select(Diet,starts_with("m"))%>%select_if(~ is.character(.) || any(. > 0))%>%
    gather(m,val,-Diet)%>%filter(m%in%mid)%>%group_by(m)%>%t_test(val~Diet)%>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj")
  resf<-x%>%filter(Glucose=="Glucose",Diet%in%g)%>%select(Diet,starts_with("m"))%>%
    select_if(~ is.character(.) || any(. > 0))%>%
    gather(m,val,-Diet)%>%group_by(Diet,m)%>%filter(m%in%mid)%>%summarise(mu=mean(val))%>%
    spread(Diet,mu)
  return(cbind(res,resf))
}
####
daths<-lapply(datpu, function(x)cal_ttest(x,g=c("HFD","SD")))
dathr<-lapply(datpu, function(x)cal_ttest(x,g=c("DR","HFD")))
dathex<-lapply(datpu, function(x)cal_ttest(x,g=c("EX","HFD")))
dathrex<-lapply(datpu, function(x)cal_ttest(x,g=c("DR-EX","HFD")))
datrex<-lapply(datpu, function(x)cal_ttest(x,g=c("DR","DR-EX")))
datex<-lapply(datpu, function(x)cal_ttest(x,g=c("EX","DR-EX")))
datrrex<-lapply(datpu, function(x)cal_ttest(x,g=c("DR","EX")))
#######
sapply(names(daths), function(x)write.csv(daths[[x]],file=paste0(x,"_glucose_HFD_SD_ttest_m0_sum.csv")))
sapply(names(dathr), function(x)write.csv(dathr[[x]],file=paste0(x,"_glucose_HFD_DR_ttest_m0_sum.csv")))
sapply(names(dathex), function(x)write.csv(dathex[[x]],file=paste0(x,"_glucose_HFD_EX_ttest_m0_sum.csv")))
sapply(names(dathrex), function(x)write.csv(dathrex[[x]],file=paste0(x,"_glucose_HFD_DR_EX_ttest_m0_sum.csv")))
sapply(names(datex), function(x)write.csv(datex[[x]],file=paste0(x,"_glucose_EX_DR_EX_ttest_m0_sum.csv")))
sapply(names(datrrex), function(x)write.csv(datrrex[[x]],file=paste0(x,"_glucose_DR_EX_ttest_m0_sum.csv")))
######
##########
hs<-do.call(rbind,daths)
hs$Group<-sub('_dat.*','',rownames(hs))
hr<-do.call(rbind,dathr)
hr$Group<-sub('_dat.*','',rownames(hr))
hex<-do.call(rbind,dathex)
hex$Group<-sub('_dat.*','',rownames(hex))
hrex<-do.call(rbind,dathrex)
hrex$Group<-sub('_dat.*','',rownames(hrex))
ex<-do.call(rbind,datex)
ex$Group<-sub('_dat.*','',rownames(ex))
rrex<-do.call(rbind,datrrex)
rrex$Group<-sub('_dat.*','',rownames(rrex))
#################
hs$log2FC<-log2(hs$HFD/hs$SD)
hr$log2FC<-log2(hr$DR/hr$HFD)
hex$log2FC<-log2(hex$EX/hex$HFD)
hrex$log2FC<-log2(hrex$`DR-EX`/hrex$HFD)
ex$log2FC<-log2(ex$`DR-EX`/ex$EX)
rrex$log2FC<-log2(rrex$DR/rrex$EX)
#########
res<-list(HFDvsSD=hs,DRvsHFD=hr,EXvsHFD=hex,DREXvsHFD=hrex,DREXvsEX=ex,DRvsEX=rrex)
###
sapply(names(res), function(x)write.csv(res[[x]],file=paste0(x,"_ttest_logFC_m0_sum.csv")))
ress<-lapply(res, function(x)subset(x,p<0.05))
###
resd<-lapply(ress, function(x)x%>%select(m,p,log2FC))
resdd<-do.call(rbind,resd)
resdd$Group<-sub('\\..*','',rownames(resdd))
resdd$Glucose<-sub('_dat','',sub('.*\\.','',sub('_dat\\.\\d+','',rownames(resdd))))
ggplot(resdd,aes(m,Group,color=log2FC))+geom_point(size=5,pch=20)+geom_point(size=7)+
  scale_color_gradient2(low="cyan4",mid="grey",high = "red",midpoint = 0)+
  geom_text(aes(label=round(log2FC,2)),size=2.5,color="white")+
  theme_light(base_size = 14)+xlab("")+ylab("")+
  theme(axis.text.x=element_text(angle=90,hjust = 0.5,vjust = 1))+facet_wrap(~Glucose,ncol = 5)
dev.print(pdf,file="all_sig_dot_m0_sum.pdf")
###############
library(rcolors)
mycol<-c("#B3D6AD","#A17DB4")
names(mycol)<-c("SD","HFD")

### manually added by JH - not sure if these are accurate colors for
#distcolor <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00")

distcolor <- c(
  "SD"   = "#CD338B",  "HFD"  = "#389182",  "DR"   = "#AA7322",
  "EX"   = "#29A1CA",  "DR-EX"= "#E17100"
)

###
genef<-function(x){
  d<-ungroup(datpu[[x]])
  p<-d%>%filter(Glucose=="Glucose")%>%select(Diet,starts_with("m"))%>%
    select_if(~ is.character(.) || any(. > 0))%>%gather(m,val,-Diet)%>%
    # Set Diet order for x axis AND legend
    mutate(Diet = factor(Diet, levels = c("SD","HFD","DR","EX","DR-EX"))) %>%
    mutate(or=sub('.*\\+','',m))%>%mutate(m=fct_reorder(m,as.numeric(or)))%>%
    ggplot(aes(Diet,val,fill=Diet))+geom_boxplot()+facet_wrap(~m,scale="free")+
    theme_classic()+scale_fill_manual(values=distcolor[1:5])+ylab("Percentage")+xlab("")
  p+xlim(c("SD","HFD","DR","EX","DR-EX"))
}
####
for(i in names(datpu)){
  genef(i)
  ggsave(paste0(i,"_glucose_boxplot_sum.pdf"),width = 10.8,height = 8)
}
#######
heat<-function(x){
  x<-ungroup(x)
  d<-x%>%filter(Diet!="Pool",Glucose=="Glucose")%>%select(Diet,starts_with("m"),-`m+0`)%>%
    gather(m,val,-Diet)%>%group_by(Diet,m)%>%summarise(mu=mean(val))%>%spread(Diet,mu)%>%separate(m,c("sam","val"),sep="\\+")%>%
    arrange(as.numeric(val))%>%mutate(m=paste0(sam,"+",val))%>%select(-sam,-val)%>%column_to_rownames(var="m")
  d<-d[,intersect(c("SD","HFD","DR","EX","DR-EX"),colnames(d))]
  pheatmap::pheatmap(d,scale="row",cluster_cols = F,color=colorRampPalette(c("deepskyblue","white","darkorange"))(128),
                     border_color = "white",cluster_rows = F,cellwidth = 18,cellheight = 18)
}
###
for(i in names(datpu)){
  heat(datpu[[i]])
  dev.print(pdf,paste0(i,"_glucose_heat_wo_sum.pdf"))
}
# #####
# heatp<-function(x,text=TRUE){
#   d<-datt[[x]]
#   d<-d%>%setNames(make.names(names(.), unique = TRUE))
#   d<-d%>%select(group,m,logFC,p)
#   d$label<-scientific(d$logFC, digits = 2)
#   d<-d%>%mutate(sig=ifelse(p>0.05|is.na(p),"","*"))
#   if(isTRUE(text)){
#     d$label<-paste0(d$label,d$sig)
#   }else{
#     d$label<-d$sig
#   }
#   d[is.na(d$logFC),"logFC"]<-0
#   d[is.infinite(d$logFC)&d$logFC<0,"logFC"]<-min(d$logFC[!is.infinite(d$logFC)])-1
#   d[is.infinite(d$logFC)&d$logFC>0,"logFC"]<-max(d$logFC[!is.infinite(d$logFC)])+1
#   p<-d%>%mutate(or=sub('.*\\+','',m))%>%
#     mutate(m=fct_reorder(m,as.numeric(or)))%>%ggplot(aes(group,m,fill=logFC))+
#     geom_tile(color="white",linewidth=0.2)+geom_text(aes(label=label),size=4)+
#     scale_fill_gradient2(high="red",low="cyan4",mid="white",midpoint = 0)+theme_minimal()
#   p+xlab("")+ylab("")
# }
# ####
# for(i in names(datt)){
#   heatp(i)+xlab("")+ylab("")+labs(fill="log2FC")
#   ggsave(paste0(i,"_glucose_heat.pdf"))
# }

###########mfuzz

dat_avg<-function(name){
  x<-ungroup(datpu[[name]])
  mid<-x%>%filter(Glucose=="Glucose")%>%select(Diet,starts_with("m"))%>%select_if(~ is.character(.) || any(. > 0))%>%
    gather(m,val,-Diet)%>%group_by(Diet,m)%>%
    ###### filter 1% at least in 1/3 samples
    summarise(count=n(),cc=sum(val>0.01))%>%filter(cc/count>0.3)%>%pull(m)%>%unique()
  res<-x%>%filter(Glucose=="Glucose")%>%select(Diet,starts_with("m"))%>%select_if(~ is.character(.) || any(. > 0))%>%
    gather(m,val,-Diet)%>%filter(m%in%mid)%>%group_by(Diet,m)%>%summarise(mu=mean(val))%>%spread(m,mu)%>%column_to_rownames(var="Diet")
  colnames(res)<-paste0(name,"_",sub('\\+','_',colnames(res)))
  return(res)
}
du<-do.call(cbind,lapply(names(datpu), function(x)dat_avg(x)))
####NP only
dat_avgp<-function(name){
  x<-ungroup(datp[[name]])
  mid<-x%>%filter(Glucose=="Glucose",Group=="np")%>%select(Diet,starts_with("m"))%>%select_if(~ is.character(.) || any(. > 0))%>%
    gather(m,val,-Diet)%>%group_by(Diet,m)%>%
    ###### filter 1% at least in 1/3 samples
    summarise(count=n(),cc=sum(val>0.01))%>%filter(cc/count>0.3)%>%pull(m)%>%unique()
  res<-x%>%filter(Glucose=="Glucose",Group=="np")%>%select(Diet,starts_with("m"))%>%select_if(~ is.character(.) || any(. > 0))%>%
    gather(m,val,-Diet)%>%filter(m%in%mid)%>%group_by(Diet,m)%>%summarise(mu=mean(val))%>%spread(m,mu)%>%column_to_rownames(var="Diet")
  colnames(res)<-paste0(name,"_",sub('\\+','_',colnames(res)))
  return(res)
}
dup<-do.call(cbind,lapply(names(datp), function(x)dat_avgp(x)))
write.csv(dup,file="fluxomics_avg_np.csv")
dat_avgd<-function(name){
  x<-ungroup(datp[[name]])
  mid<-x%>%filter(Glucose=="Glucose",Group=="nd")%>%select(Diet,starts_with("m"))%>%select_if(~ is.character(.) || any(. > 0))%>%
    gather(m,val,-Diet)%>%group_by(Diet,m)%>%
    ###### filter 1% at least in 1/3 samples
    summarise(count=n(),cc=sum(val>0.01))%>%filter(cc/count>0.3)%>%pull(m)%>%unique()
  res<-x%>%filter(Glucose=="Glucose",Group=="nd")%>%select(Diet,starts_with("m"))%>%select_if(~ is.character(.) || any(. > 0))%>%
    gather(m,val,-Diet)%>%filter(m%in%mid)%>%group_by(Diet,m)%>%summarise(mu=mean(val))%>%spread(m,mu)%>%column_to_rownames(var="Diet")
  colnames(res)<-paste0(name,"_",sub('\\+','_',colnames(res)))
  return(res)
}
dud<-do.call(cbind,lapply(names(datp), function(x)dat_avgd(x)))
write.csv(dup,file="fluxomics_avg_nd.csv")


###
library(Mfuzz)
#############Mfuzz
du<-t(du)
mname<-c("SD","HFD","DR","EX","DR-EX")
mat<-du[,mname]
library(Mfuzz)
eset <- new("ExpressionSet", exprs=as.matrix(mat))
eset <- standardise(eset)
m <- mestimate(eset)
#######
###### Function to calculate minimum distance between centroids
min_dist_between_centroids <- function(data, c, m) {
  cl <- mfuzz(data, c = c, m = m)
  centroids <- cl$centers
  dist_matrix <- as.matrix(dist(centroids))
  diag(dist_matrix) <- NA
  return(min(dist_matrix, na.rm = TRUE))
}
# Calculate minimum distances
min_distances <- sapply(2:10, function(c) min_dist_between_centroids(eset, c, m))
# Plot the minimum distances
plot(2:10, min_distances, type = "b", xlab = "Number of clusters", ylab = "Minimum Distance Between Centroids")
dev.print(pdf,file="DREX28_minimum_distance_mfuzz.pdf")

######
library(clusterSim)
library(fpc)
# Function to calculate Dunn's index
dunn_index <- function(data, c, m) {
  cl <- mfuzz(data, c = c, m = m)
  # Extract clustering membership
  cluster_membership <- apply(cl$membership, 1, which.max)
  return(cluster.stats(d = dist(exprs(data)), clustering = cluster_membership)$dunn)
}
set.seed(1234)
# Calculate Dunn's index for different cluster numbers
dunn_indices <- sapply(2:10, function(c) dunn_index(eset, c, m))

# Plot the Dunn's indices
plot(2:10, dunn_indices, type = "b", xlab = "Number of clusters", ylab = "Dunn's Index")
dev.print(pdf,file="DREX_Dunn_distance_mfuzz.pdf")
#####
c <- 5  # Number of clusters
cl <- mfuzz(eset, c=c, m=m)
mfuzz.plot(eset, cl=cl, mfrow=c(2,2), time.labels=mname, new.window=FALSE)
dev.print(pdf,file="DREX_mfuzz_cluster.pdf",width=6.2,height=4.2)
mat<-as.data.frame(mat)
mat$cluster<-cl$cluster[rownames(mat)]
matg<-split(rownames(mat),mat$cluster)
write.csv(mat,file="mfuzz_cluster.csv")
######

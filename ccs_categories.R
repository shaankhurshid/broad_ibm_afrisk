# Script to generate multi-level CCS categories from single entries

# Dependencies
library(data.table)

# Load CCS master
ccs_master <- fread(file='/Volumes/medpop_afib/skhurshid/broad_afrisk/CCS_2015_unifed_august_1_2017.csv')
setnames(ccs_master,c('code','ccs','code_type'))

# Create hierarchy
## Names
ccs_master[,':='(ccs_name_level1 = ifelse((ccs %in% 1:10),'infectious_parasitic_diseases',
                 ifelse((ccs %in% 11:47),'neoplasms',
                 ifelse((ccs %in% 48:58),'endocrine_nutritional_metabolic_immunity',
                 ifelse((ccs %in% 59:64),'blood_diseases',
                 ifelse((ccs %in% c(650:663,670)),'mental_illess',
                 ifelse((ccs %in% 76:95),'nervous_system',
                 ifelse((ccs %in% 96:121),'circulatory',
                 ifelse((ccs %in% 122:134),'respiratory',
                 ifelse((ccs %in% 135:155),'digestive',
                 ifelse((ccs %in% 156:175),'genitourinary',
                 ifelse((ccs %in% 176:196),'pregnancy_related',
                 ifelse((ccs %in% 197:200),'dermatology',
                 ifelse((ccs %in% 201:212),'msk',
                 ifelse((ccs %in% 213:217),'congenital',
                 ifelse((ccs %in% 218:224),'perinatal',
                 ifelse((ccs %in% 225:244),'injury_poisoning',
                 ifelse((ccs %in% 245:258),'symptoms_signs_illdefined','residual'))))))))))))))))))]

## Codes
ccs_master[,':='(ccs_level1 = ifelse((ccs %in% 1:10),1,
  ifelse((ccs %in% 11:47),2,
  ifelse((ccs %in% 48:58),3,
  ifelse((ccs %in% 59:64),4,
  ifelse((ccs %in% c(650:663,670)),5,
  ifelse((ccs %in% 76:95),6,
  ifelse((ccs %in% 96:121),7,
  ifelse((ccs %in% 122:134),8,
  ifelse((ccs %in% 135:155),9,
  ifelse((ccs %in% 156:175),10,
  ifelse((ccs %in% 176:196),11,
  ifelse((ccs %in% 197:200),12,
  ifelse((ccs %in% 201:212),13,
  ifelse((ccs %in% 213:217),14,
  ifelse((ccs %in% 218:224),15,
  ifelse((ccs %in% 225:244),16,
  ifelse((ccs %in% 245:258),17,18))))))))))))))))))]

# Create sub-hierarchy (THIS NEEDS TO BE BATCHED BECAUSE TOO MUCH DATA)
## Names
ccs_master[,':='(ccs_name_level2 = 
  ifelse((ccs %in% c(1:3,9)),'bacterial_infection',
  ifelse((ccs %in% 4),'mycoses',
  ifelse((ccs %in% 5:7),'viral_infection',
  ifelse((ccs %in% 8),'other_infection',
  ifelse((ccs %in% 10),'immunization_id_screening',
  ifelse((ccs %in% 14:15),'colon_cancer',
  ifelse((ccs %in% c(12,13,16:18)),'other_gi_cancer',
  ifelse((ccs %in% 19),'lung_cancer',
  ifelse((ccs %in% 22:23),'skin_cancer',
  ifelse((ccs %in% 24),'breast_cancer',
  ifelse((ccs %in% 25:26),'uterine_cancer',
  ifelse((ccs %in% 27:28),'other_gyn_cancer',
  ifelse((ccs %in% 29:31),'male_gu_cancer',
  ifelse((ccs %in% 32:34),'urinary_cancer',
  ifelse((ccs %in% 37:40),'lymph_heme_cancer',
  ifelse((ccs %in% c(11,20,21,35,36,41)),'other_cancer',
  ifelse((ccs %in% 42),'secondary_cancer',
  ifelse((ccs %in% 43),'cancer_nosite',
  ifelse((ccs %in% 44),'cancer_unspecified',
  ifelse((ccs %in% 45),'cancer_chemo_xrt',
  ifelse((ccs %in% 46:47),'benign_neoplasm',
  ifelse((ccs %in% 48),'thyroid',
  ifelse((ccs %in% 49),'dm_no_complication',
  ifelse((ccs %in% 50),'dm_complication', 
  ifelse((ccs %in% 51),'other_endocrine',
  ifelse((ccs %in% 52),'nutritional',
  ifelse((ccs %in% 53),'lipid',
  ifelse((ccs %in% 54),'gout',
  ifelse((ccs %in% 55),'fluid_electrolytes',
  ifelse((ccs %in% 56),'cystic_fibrosis',
  ifelse((ccs %in% 57),'immunity_disorder',              
  ifelse((ccs %in% 58),'other_nutrtional_endocrine',
  ifelse((ccs %in% 59:61),'anemia',
  ifelse((ccs %in% 62),'coagulation_hemorrhagic',
  ifelse((ccs %in% 63),'wbc_disease',                             
  ifelse((ccs %in% 64),'rbc_disease',
  ifelse((ccs %in% 650),'adjustment_disorders',
  ifelse((ccs %in% 651),'anxiety_disorders',
  ifelse((ccs %in% 652),'attention_conduct_disorders',              
  ifelse((ccs %in% 653),'deliirum_dementia',
  ifelse((ccs %in% 654),'developmental_disorders',
  ifelse((ccs %in% 655),'childhood_disorders',
  ifelse((ccs %in% 656),'unspecified_impulse_disorders','placeholder'))))))))))))))))))))))))))))))))))))))))))))]
  
ccs_master[,':='(ccs_name_level2 = 
  ifelse((ccs %in% 657),'mood_disorders',   
  ifelse((ccs %in% 658),'personality_disorders',   
  ifelse((ccs %in% 659),'schizophrenia',   
  ifelse((ccs %in% 660),'alcohol_use_disorders',   
  ifelse((ccs %in% 661),'substance_use_disorders',   
  ifelse((ccs %in% 662),'suicide',   
  ifelse((ccs %in% 663),'screening_history_mental_illness',   
  ifelse((ccs %in% 670),'misc_mental_disorders',   
  ifelse((ccs %in% 76:78),'cns_infection',
  ifelse((ccs %in% 79:81),'cns_hereditary_degenerative',   
  ifelse((ccs %in% 82),'paralysis',   
  ifelse((ccs %in% 83),'epilepsy_convulsions',   
  ifelse((ccs %in% 84),'headache',   
  ifelse((ccs %in% 85),'coma',   
  ifelse((ccs %in% 86:91),'eye_disorders',   
  ifelse((ccs %in% 92:94),'ear_disorders',   
  ifelse((ccs %in% 95),'cns_other',   
  ifelse((ccs %in% 98:99),'hypertension',   
  ifelse((ccs %in% 96:108),'heart_diseases',   
  ifelse((ccs %in% 109:113),'cerebrovascular_diseases',
  ifelse((ccs %in% 114:117),'arterial_diseases',   
  ifelse((ccs %in% 118:121),'vein_lymphatic_diseases',   
  ifelse((ccs %in% 122:126),'respiratory_infection',   
  ifelse((ccs %in% 127),'copd',   
  ifelse((ccs %in% 128),'asthma',   
  ifelse((ccs %in% 129),'aspiration',   
  ifelse((ccs %in% 130),'pleurisy_ptx',   
  ifelse((ccs %in% 131),'respiratory_failure_adult',   
  ifelse((ccs %in% 132),'lung_disease_external_agents',
  ifelse((ccs %in% 133),'other_lower_respiratory_disease',   
  ifelse((ccs %in% 134),'other_upper_respiratory_disease',   
  ifelse((ccs %in% 135),'intestinal_infection',   
  ifelse((ccs %in% 136),'teeth_jaw_disorders',   
  ifelse((ccs %in% 137),'mouth_disease',   
  ifelse((ccs %in% 138:141),'upper_gi_disorders',   
  ifelse((ccs %in% 143),'hernia',(ccs_name_level2))))))))))))))))))))))))))))))))))))))]   
  
ccs_master[,':='(ccs_name_level2 = 
  ifelse((ccs %in% c(142,144:148)),'lower_gi_disorders',   
  ifelse((ccs %in% 149),'biliary_disease',
  ifelse((ccs %in% 150:151),'liver_disease',   
  ifelse((ccs %in% 152),'pancreatic_disease',   
  ifelse((ccs %in% 153),'gi_hemorrhage',   
  ifelse((ccs %in% 154),'non_infectious_gastroenteritis',   
  ifelse((ccs %in% 155),'other_gi_disorders',   
  ifelse((ccs %in% 156:163),'urinary_disease',   
  ifelse((ccs %in% 164:166),'male_genital_disease',   
  ifelse((ccs %in% 167:175),'female_genital_disease',   
  ifelse((ccs %in% 176),'contraceptive_management',   
  ifelse((ccs %in% 177:179),'abortion_related_disorders',       
  ifelse((ccs %in% 180:186),'pregnancy_disorders',  
  ifelse((ccs %in% 187:191),'labor_delivery',  
  ifelse((ccs %in% 192:194),'labor_complications',  
  ifelse((ccs %in% 195),'other_birth_complications',  
  ifelse((ccs %in% 196),'normal_pregnancy',
  ifelse((ccs %in% 197),'skin_infection',  
  ifelse((ccs %in% 198),'skin_inflammation',
  ifelse((ccs %in% 199),'skin_ulcer',
  ifelse((ccs %in% 200),'other_skin_disorder',
  ifelse((ccs %in% 201),'infective_arthritis',  
  ifelse((ccs %in% 202:204),'non_traumatic_arthritis',  
  ifelse((ccs %in% 205),'spondylosis',  
  ifelse((ccs %in% 206),'osteoporosis',  
  ifelse((ccs %in% 207),'pathologic_fracture',  
  ifelse((ccs %in% 208:209),'acquired_deformity',  
  ifelse((ccs %in% 210),'lupus',  
  ifelse((ccs %in% 211),'other_connective_tissue_disease',  
  ifelse((ccs %in% 212),'other_bone_disease',(ccs_name_level2))))))))))))))))))))))))))))))))] 
  
ccs_master[,':='(ccs_name_level2 = 
  ifelse((ccs %in% 213),'cardiac_congenital_anomaly',  
  ifelse((ccs %in% 214),'gi_congenital_anomaly',  
  ifelse((ccs %in% 215),'gu_congenital_anomaly',  
  ifelse((ccs %in% 216),'cns_congenital_anomaly',  
  ifelse((ccs %in% 217),'other_congenital_anomaly',  
  ifelse((ccs %in% 218),'liveborn',  
  ifelse((ccs %in% 219),'short_gestation',  
  ifelse((ccs %in% 220),'intrauterine_hypoxia',  
  ifelse((ccs %in% 221),'respiratory_distress_syndrome',  
  ifelse((ccs %in% 222),'perinatal_jaundice',
  ifelse((ccs %in% 223),'birth_trauma',  
  ifelse((ccs %in% 224),'other_perinatal',  
  ifelse((ccs %in% 225),'joint_dislocation',  
  ifelse((ccs %in% c(226,228:231)),'fracture',
  ifelse((ccs %in% 227),'spinal_cord_injury',  
  ifelse((ccs %in% 233),'intracranial_injury',
  ifelse((ccs %in% 234),'crush_injury', 
  ifelse((ccs %in% 235:236),'open_wound',  
  ifelse((ccs %in% 232),'sprain_strain',  
  ifelse((ccs %in% 239),'contusion',  
  ifelse((ccs %in% 240),'burn',  
  ifelse((ccs %in% 237:238),'complication',  
  ifelse((ccs %in% 241:243),'poisoning',  
  ifelse((ccs %in% 244),'other_injury',  
  ifelse((ccs %in% 245:253),'symptoms_signs_ill_defined',  
  ifelse((ccs %in% 254:258),'healthcare_factors',
  ifelse((ccs %in% 259:260),'residual',
  ifelse((ccs_level1 == 18),'residual',(ccs_name_level2))))))))))))))))))))))))))))))]

## Codes
ccs_master[,':='(ccs_level2 = as.character(
ifelse((ccs %in% c(1:3,9)),1.1,
ifelse((ccs %in% 4),1.2,
ifelse((ccs %in% 5:7),1.3,
ifelse((ccs %in% 8),1.4,
ifelse((ccs %in% 10),1.5,
ifelse((ccs %in% 14:15),2.1,
ifelse((ccs %in% c(12,13,16:18)),2.2,
ifelse((ccs %in% 19),2.3,
ifelse((ccs %in% 22:23),2.4,
ifelse((ccs %in% 24),2.5,
ifelse((ccs %in% 25:26),2.6,
ifelse((ccs %in% 27:28),2.7,
ifelse((ccs %in% 29:31),2.8,
ifelse((ccs %in% 32:34),2.9,
ifelse((ccs %in% 37:40),2.10,
ifelse((ccs %in% c(11,20,21,35,36,41)),2.11,
ifelse((ccs %in% 42),2.12,
ifelse((ccs %in% 43),2.13,
ifelse((ccs %in% 44),2.14,
ifelse((ccs %in% 45),2.15,
ifelse((ccs %in% 46:47),2.16,
ifelse((ccs %in% 48),3.1,
ifelse((ccs %in% 49),3.2,
ifelse((ccs %in% 50),3.3, 
ifelse((ccs %in% 51),3.4,
ifelse((ccs %in% 52),3.5,
ifelse((ccs %in% 53),3.6,
ifelse((ccs %in% 54),3.7,
ifelse((ccs %in% 55),3.8,
ifelse((ccs %in% 56),3.9,
ifelse((ccs %in% 57),3.10,              
ifelse((ccs %in% 58),3.11,
ifelse((ccs %in% 59:61),4.1,
ifelse((ccs %in% 62),4.2,
ifelse((ccs %in% 63),4.3,                             
ifelse((ccs %in% 64),4.4,
ifelse((ccs %in% 650),5.1,
ifelse((ccs %in% 651),5.2,
ifelse((ccs %in% 652),5.3,              
ifelse((ccs %in% 653),5.4,
ifelse((ccs %in% 654),5.5,
ifelse((ccs %in% 655),5.6,
ifelse((ccs %in% 656),5.7,NA)))))))))))))))))))))))))))))))))))))))))))))]

ccs_master[,':='(ccs_level2 = as.character(
ifelse((ccs %in% 657),5.8,   
ifelse((ccs %in% 658),5.9,   
ifelse((ccs %in% 659),5.10,   
ifelse((ccs %in% 660),5.11,   
ifelse((ccs %in% 661),5.12,   
ifelse((ccs %in% 662),5.13,   
ifelse((ccs %in% 663),5.14,   
ifelse((ccs %in% 670),5.15,   
ifelse((ccs %in% 76:78),6.1,
ifelse((ccs %in% 79:81),6.2,   
ifelse((ccs %in% 82),6.3,   
ifelse((ccs %in% 83),6.4,   
ifelse((ccs %in% 84),6.5,   
ifelse((ccs %in% 85),6.6,   
ifelse((ccs %in% 86:91),6.7,   
ifelse((ccs %in% 92:94),6.8,   
ifelse((ccs %in% 95),6.9,   
ifelse((ccs %in% 98:99),7.1,   
ifelse((ccs %in% 96:108),7.2,   
ifelse((ccs %in% 109:113),7.3,
ifelse((ccs %in% 114:117),7.4,   
ifelse((ccs %in% 118:121),7.5,   
ifelse((ccs %in% 122:126),8.1,   
ifelse((ccs %in% 127),8.2,   
ifelse((ccs %in% 128),8.3,   
ifelse((ccs %in% 129),8.4,   
ifelse((ccs %in% 130),8.5,   
ifelse((ccs %in% 131),8.6,   
ifelse((ccs %in% 132),8.7,
ifelse((ccs %in% 133),8.8,   
ifelse((ccs %in% 134),8.9,   
ifelse((ccs %in% 135),9.1,   
ifelse((ccs %in% 136),9.2,   
ifelse((ccs %in% 137),9.3,   
ifelse((ccs %in% 138:141),9.4,   
ifelse((ccs %in% 143),9.5,(ccs_level2)))))))))))))))))))))))))))))))))))))))]   

ccs_master[,':='(ccs_level2 = as.character(
ifelse((ccs %in% c(142,144:148)),9.6,   
ifelse((ccs %in% 149),9.7,
ifelse((ccs %in% 150:151),9.8,   
ifelse((ccs %in% 152),9.9,   
ifelse((ccs %in% 153),9.10,   
ifelse((ccs %in% 154),9.11,   
ifelse((ccs %in% 155),9.12,   
ifelse((ccs %in% 156:163),10.1,   
ifelse((ccs %in% 164:166),10.2,   
ifelse((ccs %in% 167:175),10.3,   
ifelse((ccs %in% 176),11.1,   
ifelse((ccs %in% 177:179),11.2,       
ifelse((ccs %in% 180:186),11.3,  
ifelse((ccs %in% 187:191),11.4,  
ifelse((ccs %in% 192:194),11.5,  
ifelse((ccs %in% 195),11.6,  
ifelse((ccs %in% 196),11.7,
ifelse((ccs %in% 197),12.1,  
ifelse((ccs %in% 198),12.2,
ifelse((ccs %in% 199),12.3,
ifelse((ccs %in% 200),12.4,
ifelse((ccs %in% 201),13.1,  
ifelse((ccs %in% 202:204),13.2,  
ifelse((ccs %in% 205),13.3,  
ifelse((ccs %in% 206),13.4,  
ifelse((ccs %in% 207),13.5,  
ifelse((ccs %in% 208:209),13.6,  
ifelse((ccs %in% 210),13.7,  
ifelse((ccs %in% 211),13.8,  
ifelse((ccs %in% 212),13.9,(ccs_level2)))))))))))))))))))))))))))))))))] 

ccs_master[,':='(ccs_level2 = as.character(
ifelse((ccs %in% 213),14.1,  
ifelse((ccs %in% 214),14.2,  
ifelse((ccs %in% 215),14.3,  
ifelse((ccs %in% 216),14.4,  
ifelse((ccs %in% 217),14.5,  
ifelse((ccs %in% 218),15.1,  
ifelse((ccs %in% 219),15.2,  
ifelse((ccs %in% 220),15.3,  
ifelse((ccs %in% 221),15.4,  
ifelse((ccs %in% 222),15.5,
ifelse((ccs %in% 223),15.6,  
ifelse((ccs %in% 224),15.7,  
ifelse((ccs %in% 225),16.1,  
ifelse((ccs %in% c(226,228:231)),16.2,
ifelse((ccs %in% 227),16.3,  
ifelse((ccs %in% 233),16.4,
ifelse((ccs %in% 234),16.5,      
ifelse((ccs %in% 235:236),16.6,  
ifelse((ccs %in% 232),16.7,  
ifelse((ccs %in% 239),16.8,  
ifelse((ccs %in% 240),16.9,  
ifelse((ccs %in% 237:238),16.10,  
ifelse((ccs %in% 241:243),16.11,  
ifelse((ccs %in% 244),16.12,  
ifelse((ccs %in% 245:253),17.1,  
ifelse((ccs %in% 254:258),17.2,
ifelse((ccs %in% 259:260),18,
ifelse((ccs_level1 == 18),18,(ccs_level2)))))))))))))))))))))))))))))))]

# Add names for raw ccs classes
## Load lookup table
ccs_lookup <- fread(file='/Volumes/medpop_afib/skhurshid/broad_afrisk/CCS_2015_categories_august_1_2017.csv')
setDF(ccs_lookup); setDF(ccs_master)

## Loop over and select names
ccs_master$ccs_name <- rep(NA,nrow(ccs_master))
for (i in 1:nrow(ccs_master)){
  exit_j <- 0
  for (j in 1:length(ccs_lookup$V1)){
    if (exit_j == 1) {
      break
    } else if (ccs_master$ccs[i] == ccs_lookup$V1[j]) {
      ccs_master$ccs_name[i] <- ccs_lookup$V2[j]
      exit_j <- 1
    }}
}

# Format ICD codes for PHS 
## Put in the period
ccs_master$code <- with(ccs_master,ifelse(nchar(code)==3,code,
                            ifelse(nchar(code)==4,paste0(substr(code,1,3),'.',substr(ccs_master$code,4,4)),
                                   ifelse(nchar(code)==5,paste0(substr(code,1,3),'.',substr(ccs_master$code,4,5)),
                                          ifelse(nchar(code)==6,paste0(substr(code,1,3),'.',substr(ccs_master$code,4,6)),
                                                 ifelse(nchar(code)==7,paste0(substr(code,1,3),'.',substr(ccs_master$code,4,7)),NA))))))

## Reformat code types
ccs_master[,code_type := paste0('ICD',as.character(code_type))]

## Add indicator for level of CCS for disambiguation
ccs_master[,':='(ccs_name = paste0('ccs3_',ccs_name),
                 ccs_name_level1 = paste0('ccs1_',ccs_name_level1),
                 ccs_name_level2 = paste0('ccs2_',ccs_name_level2))]

# Parse into subcomponents for covariate make files
ccs_individual <- ccs_master[,c('ccs_name','code','code_type')]
setnames(ccs_individual,c('ccs_name','code','code_type'),c('cov','cov.code','cov.code.type'))

ccs_subgroup <- ccs_master[,c('ccs_name_level2','code','code_type')]
setnames(ccs_subgroup,c('ccs_name_level2','code','code_type'),c('cov','cov.code','cov.code.type'))

ccs_group <- ccs_master[,c('ccs_name_level1','code','code_type')]
setnames(ccs_group,c('ccs_name_level1','code','code_type'),c('cov','cov.code','cov.code.type'))

ccs_all <- rbind(ccs_group,ccs_subgroup,ccs_individual)

# Save out
save(ccs_master,file='/Volumes/medpop_afib/skhurshid/broad_afrisk/ccs_master_102519.RData')
write.csv(ccs_master,file='/Volumes/medpop_afib/skhurshid/broad_afrisk/ccs_master_102519.csv')

# Save out subcomponents for covariate make files
write.csv(ccs_individual,file='/Volumes/homedir$/MGH Research/IBM_Broad_afrisk_validation/ccs_individual_102519.csv')
write.csv(ccs_subgroup,file='/Volumes/homedir$/MGH Research/IBM_Broad_afrisk_validation/ccs_subgroup_102519.csv')
write.csv(ccs_group,file='/Volumes/homedir$/MGH Research/IBM_Broad_afrisk_validation/ccs_group_102519.csv')
write.csv(ccs_all,file='/Volumes/homedir$/MGH Research/IBM_Broad_afrisk_validation/ccs_all_102519.csv')
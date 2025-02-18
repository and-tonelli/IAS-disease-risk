
require(tidyverse)
require(magrittr)
'%!in%' <- function(x,y)!('%in%'(x,y))

# A bit of housekeeping first
DAMA_tab <- read.csv("data/DataS1/DAMA_table.csv")
DAMA_tab$CommonName[DAMA_tab$CommonName == "Cotton-Top Tamarins"] <- "Cotton-Top Tamarin" 

Genbank_IAS_common_I <- read_csv("data/Genbank_IAS_common_viruses.csv") %>% select(-1) %>% filter(!is.na(species))
Genbank_IAS_common_I$species[DAMA_tab$species == "Cotton-Top Tamarins"] <- "Cotton-Top Tamarin"
Genbank_IAS_common_II <- read_csv("data/Genbank_IAS_common_virusesII.csv") %>% select(-1) %>% filter(!is.na(species))

Genbank_IAS_I <- read.csv("data/Genbank_IAS_spp_viruses.csv") %>% select(-1) %>% filter(!is.na(species))
Genbank_IAS_II <- read.csv("data/Genbank_IAS_spp_virusesII.csv") %>% select(-1) %>% filter(!is.na(species))
Genbank_IAS_III <- read.csv("data/Genbank_IAS_spp_virusesIII.csv") %>% select(-1) %>% filter(!is.na(species))

Genbank_IAS_species <- rbind(Genbank_IAS_I, Genbank_IAS_II, Genbank_IAS_III) %>% merge.data.frame(DAMA_tab[c(4, 5)], by = 1)
Genbank_IAS_species <- Genbank_IAS_species %>% rename("Binomial" = species)
Genbank_IAS_common <- rbind(Genbank_IAS_common_I, Genbank_IAS_common_II)%>% merge.data.frame(DAMA_tab[c(4, 5)], by.x = 1, by.y = 2)
Genbank_IAS_common <- Genbank_IAS_common %>% rename("CommonName" = species)

Genbank_IAS_final <- rbind(Genbank_IAS_common, Genbank_IAS_species) %>% relocate(Binomial)
write.csv(Genbank_IAS_final, "data/Genbank_IAS_raw.csv", row.names = F)
Genbank_IAS_final <- read.csv("data/Genbank_IAS_raw.csv")

length(unique(Genbank_IAS_final$Binomial))

# Refine countries of Genbank records
# NAs
length(Genbank_IAS_final$Binomial[is.na(Genbank_IAS_final$Location)]) # 64,023 records
Genbank_IAS_final %>% filter(is.na(Location)) %>% View() 

# Dataframe of countries where mammals have been found
Countries_mammals <- read.csv("data/Mammals_GADMcountries.csv")

cross_tab <- read.csv("data/country_crosswalk.csv")
Countries_mammals <- Countries_mammals %>% separate_rows(countries, sep = ", ")
Countries_mammals %<>% merge.data.frame(cross_tab, ., by.y = 2, by.x = 1)

# Refine countries of Genbank records
# Unmatched countries
length(Genbank_IAS_final$Binomial[Genbank_IAS_final$Location %!in% Countries_mammals$NAME_0]) #562,022 records
unique(Genbank_IAS_final$Location[Genbank_IAS_final$Location %!in% Countries_mammals$NAME_0]) #1,344 locations
# try <- "Portugal: Madeira Island"
# gsub(":.*", "", try)
Genbank_IAS_final$Location <- gsub(":.*", "", Genbank_IAS_final$Location)
unique(Genbank_IAS_final$Location[Genbank_IAS_final$Location %!in% Countries_mammals$NAME_0])

Countries_mammals$NAME_0[Countries_mammals$NAME_0 %in% "MÃ©xico"] <-  "Mexico"
Countries_mammals$NAME_0[Countries_mammals$NAME_0 %in% "CÃ´te d'Ivoire"] <-  "Cote d'Ivoire"
Countries_mammals$NAME_0[Countries_mammals$NAME_0 %in% "RÃ©union"] <- "Reunion"

# Let's fix individual names
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "USA"] <- "United States"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "Svalbard"] <- "Svalbard and Jan Mayen"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "Czech Republic"] <- "Czechia"
# Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "Saint Kitts and Nevis"] <- "St. Kitts and Nevis"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "Yugoslavia"] <- NA
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "Viet Nam"] <- "Vietnam"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "Macedonia"] <- "North Macedonia"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "Antarctica"] <- "Norway" # actually in Norway (https://link.springer.com/article/10.1186/1743-422X-2-79#Sec9) 
# Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "French Guiana"] # Technically France. Records of endemic species, shouldn't be a problem.
# Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "Dominican Republic"] <- "Dominican Rep."
# Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "Republic of the Congo"] <- "Congo"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "19.259 S 146.8169 E"] <- "Australia"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "19.259000 S 146.816900 E"] <- "Vietnam"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "53.596582 N 9.938249 E"] <- NA # zoo in Hamburg, Germany
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "45.750000 N 126.630000 E"] <- "China"
# Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "Cote d'Ivoire"] <- "CÃ´te d'Ivoire"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "22.500000 S 42.320000 W"] <- "Brazil"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "58.15304556 N 68.55260255 E"] <- "Russia"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "60.89483413 N 43.21704367 E"] <- "Russia"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "52.46206496 N 44.19864799 E"] <- "Russia"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "55.66796300 N 36.52816906 E"] <- "Russia"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "56.85553255 N 41.91563866 E"] <- "Russia"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "57.02081794 N 53.85595761 E"] <- "Russia"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "51.80429930 N 45.78379134 E"] <- "Russia"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "51.53968392 N 45.95777524 E"] <- "Russia"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "54.19112853 N 49.48268493 E"] <- "Russia"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "57.258725 N 65.129592 E"] <- "Russia"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "54.19112853 N 49.48268493 E"] <- "Russia"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "56.973581 N 65.074813 E"] <- "Russia"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "57.625918 N 66.002593 E"] <- "Russia"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "57.571226 N 67.176107 E"] <- "Russia"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "51.773315 N 45.621112 E"] <- "Russia"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "Cape Verde"] <- "Cabo Verde"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "Reunion"] # Technically France. These are records from invasive species!
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "1.3548 N 103.7763 E"] <- "Singapore"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "1.3552 N 103.7972 E"] <- "Singapore"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "1.3020 N 103.8971 E"] <- "Singapore"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "1.2819 N 103.8239 E"] <- "Singapore"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "1.2914 N 103.7667 E"] <- "Singapore"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "1.2728 N 103.8366 E"] <- "Singapore"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "1.3051 N 103.8509 E"] <- "Singapore"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "1.3312 N 103.8691 E"] <- "Singapore"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "1.30823 N 103.79868 E"] <- "Singapore"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "1.30808 N 103.79884 E"] <- "Singapore"
Genbank_IAS_final$Location[Genbank_IAS_final$Location %in% "1.27712 N 103.83399 E"] <- "Singapore"

# Filter for records outside of the native range of the species
# first checking that name match
Countries_mammals$species <- str_replace(Countries_mammals$species, "_", " ")
setdiff(Genbank_IAS_final$Binomial, Countries_mammals$species) # These domesticated species remain out: "Camelus dromedarius", "Lama glama", "Ovis orientalis", "Bubalus bubalis"

# actually filtering Genbank records
Genbank_IAS_distinct <- Genbank_IAS_final %>% distinct(Binomial, CommonName, Virus, VirusGenus, VirusFamily, Location) # 4,053
Genbank_IAS_distinct <- Genbank_IAS_distinct %>% filter(Location %in% c(Countries_mammals$NAME_0)) #3,388
Genbank_IAS_distinct <- Genbank_IAS_distinct %>% rename("NAME_0" = Location)
Countries_mammals <- Countries_mammals %>% rename("Binomial" = species)
Genbank_filtered_records <- anti_join(Genbank_IAS_distinct, Countries_mammals, by = c("Binomial", "NAME_0"))

# Ok, now these contain records in lab and zoos around the world. Let's double check against the DAMA database.
View(Genbank_filtered_records %>% filter(Binomial %!in% c("Camelus dromedarius", "Lama glama", "Ovis orientalis", "Bubalus bubalis")))

dama_countries <- read.csv("data/DAMA_GADM.csv")

setdiff(dama_countries$GID_0, Countries_mammals$GID_0)

# HongKong_mammals <- Countries_mammalsold %>% filter(countries == "Hong Kong") %>% 
#   mutate(GID_0 = "HKG",
#          COUNTRY = "Hong Kong",
#          CONTINENT = "Asia") %>% 
#   rename("NAME_0" = countries,
#          "Binomial" = species) %>% 
#   mutate(Binomial = str_replace(Binomial, "_", " "))
# 
# Countries_mammals <- rbind(Countries_mammals %>% rename("Binomial" = species), HongKong_mammals)

dama_countries2 <- dama_countries %>% 
  distinct(Binomial, GID_0) %>% 
  mutate(in_dama = "Y")

setdiff(Genbank_filtered_records$NAME_0, cross_tab$NAME_0)

cross_tab$NAME_0[cross_tab$NAME_0 %in% "MÃ©xico"] <-  "Mexico"
cross_tab$NAME_0[cross_tab$NAME_0 %in% "CÃ´te d'Ivoire"] <-  "Cote d'Ivoire"
cross_tab$NAME_0[cross_tab$NAME_0 %in% "RÃ©union"] <- "Reunion"

Genbank_filtered_records %<>% merge.data.frame(., cross_tab[c(1, 2)], by.x = 6, by.y = 2)

Genbank_filtered_records_dama <- left_join(Genbank_filtered_records %>% filter(Binomial %!in% c("Camelus dromedarius", "Lama glama", "Ovis orientalis", "Bubalus bubalis")), 
                                           dama_countries2, by = c("Binomial", "GID_0")) %>% 
  filter(in_dama == "Y")

# 215 viruses and 35 species

Genbank_filtered_records_dama %>% 
  group_by(Binomial) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n))

Genbank_filtered_records_dama %>% 
  group_by(Virus) %>% 
  summarise(n = n())

Genbank_filtered_records_dama %>% 
  group_by(GID_0) %>% 
  summarise(n = n())

write.csv(Genbank_filtered_records_dama, "data/Genbank_records_alienranges.csv", row.names = F)
Genbank_filtered_records_dama <- read.csv("data/Genbank_records_alienranges.csv")

virion_m_clean <- read_csv("data/virion_m_clean.csv")

setdiff(unique(str_to_lower(Genbank_filtered_records_dama$Virus)), str_to_lower(virion_m_clean$Virus)) 
not_found_Virus <- setdiff(unique(str_to_lower(Genbank_filtered_records_dama$Virus)), str_to_lower(virion_m_clean$VirusOriginal)) 

Genbank_filtered_records_dama$Virus <- str_to_lower(Genbank_filtered_records_dama$Virus)
virion_m_clean$VirusOriginal <- str_to_lower(virion_m_clean$VirusOriginal)

unique(Genbank_filtered_records_dama$Virus)

# Matching VIRION names of viruses
Genbank_filtered_records_dama$Virus[Genbank_filtered_records_dama$Virus == "aleutian mink disease parvovirus"] <- "aleutian mink disease virus"
Genbank_filtered_records_dama$Virus[Genbank_filtered_records_dama$Virus == "bovine-like parainfluenza virus 3"] <-  "Bovine parainfluenza virus 3"
Genbank_filtered_records_dama$Virus[Genbank_filtered_records_dama$Virus == "paslahepevirus balayani"] <- "Hepatitis E virus"
Genbank_filtered_records_dama$Virus[Genbank_filtered_records_dama$Virus %in% c("lyssavirus rabies", "rabies virus ontario fox")] <- "rabies lyssavirus"
Genbank_filtered_records_dama$Virus[Genbank_filtered_records_dama$Virus == "morbillivirus canis"] <- "canine distemper virus"
Genbank_filtered_records_dama$Virus[Genbank_filtered_records_dama$Virus %in% c("rabbit hemorrhagic disease virus 2", "rabbit calicivirus")] <- "rabbit hemorrhagic disease virus"
Genbank_filtered_records_dama$Virus[Genbank_filtered_records_dama$Virus == "alphapolyomavirus callosciuri"] <- "Callosciurus erythraeus polyomavirus 1"
Genbank_filtered_records_dama$Virus[Genbank_filtered_records_dama$Virus == "canine parvovirus 2c"] <- "canine parvovirus"
Genbank_filtered_records_dama$Virus[Genbank_filtered_records_dama$Virus == "genet fecal theilovirus"] <- "Cardiovirus B"
Genbank_filtered_records_dama$Virus[Genbank_filtered_records_dama$Virus == "hepatitis b virus genotype d"] <- "Hepatitis B virus"
Genbank_filtered_records_dama$Virus[Genbank_filtered_records_dama$Virus == "norovirus gvi"] <- "norwalk virus"
Genbank_filtered_records_dama$Virus[Genbank_filtered_records_dama$Virus == "canine parvovirus"] <- "protoparvovirus carnivoran1"
Genbank_filtered_records_dama$Virus[Genbank_filtered_records_dama$Virus %in% c("rotavirus a rva/raccoon-tc/jpn/rac-311/2011/g34p[17]", "lapine rotavirus", "simian rotavirus a strain rrv")] <- "Rotavirus A"
Genbank_filtered_records_dama$Virus[Genbank_filtered_records_dama$Virus == "cervus papillomavirus 2"] <- "Epsilonpapillomavirus 2"
Genbank_filtered_records_dama$Virus[Genbank_filtered_records_dama$Virus %in% c("simian adenovirus 58", "simian adenovirus 52", "simian adenovirus 51", "simian adenovirus 52")] <- "human mastadenovirus g"
Genbank_filtered_records_dama$Virus[Genbank_filtered_records_dama$Virus == "simian sapelovirus 3"] <- "sapelovirus b"
Genbank_filtered_records_dama$Virus[Genbank_filtered_records_dama$Virus == "betapolyomavirus macacae"] <- "betapolyomavirus macaca mulatta polyomavirus 1"

Genbank_filtered_records_dama$Virus <- str_to_lower(Genbank_filtered_records_dama$Virus)

Records_of_alien_species <- Genbank_filtered_records_dama %>% rename("VirusOriginal" = Virus) %>% merge.data.frame(., distinct(virion_m_clean[c(2, 17)]), by.x = 4, by.y = 2)

# Identifying alien viruses (found in an alien mammal in a country outside of its native range, not found in other mammals in the country)
vir_tab_tot <- NULL

for (i_sp in unique(Records_of_alien_species$Binomial)){
  
  print(paste0("Species ", which(unique(Records_of_alien_species$Binomial) == i_sp), "/", length(unique(Records_of_alien_species$Binomial))))
  
  unique(Records_of_alien_species %>% 
           filter(Binomial == i_sp) %>% 
           pull(Virus)) -> sp_viruses

  for (i_vir in c(sp_viruses)){
    
    unique(Records_of_alien_species %>% 
             filter(Binomial == i_sp, Virus == i_vir) %>% 
             pull(GID_0)) -> countries_arrived
    
    unique(virion_m_clean %>% 
             filter(Virus == i_vir,
                    Host != i_sp) %>% 
             pull(Host)) -> native_spp_w_virus
    
    unique(Countries_mammals %>% 
             filter(Binomial %in% native_spp_w_virus) %>% 
             pull(GID_0)) -> virus_native_countries
    
    if (length(virus_native_countries) > 0){
    setdiff(countries_arrived, virus_native_countries) -> new_country
    }
    
    if (length(new_country) > 0){
      
      vir_tab <- tibble(spp = i_sp,
                        virus = i_vir,
                        countries = paste(new_country, collapse = ', ')) # mind the space
      
      vir_tab_tot <- rbind(vir_tab_tot, vir_tab)
      
    }
    
  }
  
}

View(vir_tab_tot)

# It's possible that an alien virus has already spilled over into some endemic mammals (thus filtered out in the previous step). Double check and include such cases.
Tab_potential_spillover_tot <- NULL

for (i_sp in unique(Records_of_alien_species$Binomial)){
  
  print(paste0("Species ", which(unique(Records_of_alien_species$Binomial) == i_sp), "/", length(unique(Records_of_alien_species$Binomial))))
  
  unique(Records_of_alien_species %>% 
           filter(Binomial == i_sp) %>% 
           pull(Virus)) -> sp_viruses
  
  for (i_vir in c(sp_viruses)){
    
    unique(Records_of_alien_species %>% 
             filter(Binomial == i_sp, Virus == i_vir) %>% 
             pull(GID_0)) -> countries_arrived
    
    unique(virion_m_clean %>% 
             filter(Virus == i_vir,
                    Host != i_sp) %>% 
             pull(Host)) -> native_spp_w_virus
    
    Countries_mammals %>% 
      distinct(GID_0, NAME_0, Binomial) %>% 
      filter(Binomial %in% native_spp_w_virus,
             GID_0 %in% countries_arrived) -> Tab_potential_spillover
    
    
    if (nrow(Tab_potential_spillover) > 0){
      
      Tab_potential_spillover$Virus <- i_vir
      Tab_potential_spillover$Alien <- i_sp
      Tab_potential_spillover_tot <- rbind(Tab_potential_spillover_tot, Tab_potential_spillover)
      
    }
    
  }
  
}

Tab_potential_spillover_tot %>% 
  group_by(Virus, Alien) %>% 
  summarise(n = n()) %>% 
  View(.)
 
# Alien viruses that have already spilled over:
# carnivore amdoparvovirus 1 - alien everywhere but USA
# rabbit hemorrhagic disease virus - purposefully introduced in Australia
# myxoma virus - alien everywhere but Mexico and USA
# european brown hare syndrome virus - alien in Sweden (from Lepus europaeus introduced in late 19th century)
# monkeypox virus - alien in USA (by Cricetomys gambianus, it infected prairie dogs)

rescued_associations <- Tab_potential_spillover_tot %>% 
  filter((Virus %in% c("carnivore amdoparvovirus 1", "rabbit hemorrhagic disease virus", "myxoma virus")) | 
         (Virus == "european brown hare syndrome virus" & Alien == "Lepus eropaeus" & NAME_0 == "Sweden") |
         (Virus == "monkeypox virus" & Alien == "Cricetomys gambianus" & NAME_0 == "United States")) %>% 
  filter(Alien %!in% c("Vulpes vulpes"))


write.csv(rescued_associations, "data/already_spilled_associations.csv", row.names = F)

rescued_associations %<>% 
  distinct(Virus, Alien, GID_0) %>% 
  rename("countries" = GID_0,
         "spp" = Alien,
         "virus" = Virus)

vir_tab_long_final <- rbind(rescued_associations, vir_tab_tot_long)
merge.data.frame(vir_tab_long_final, cross_tab %>% distinct(GID_0, NAME_0), by.x = 3, by.y = 1)
vir_tab_long_final %<>% filter(spp %!in% c("Macaca mulatta", "Chlorocebus sabaeus")) #they come from exp inoculations in the US
write.csv(vir_tab_long_final, "data/vir_tab_long_final.csv", row.names = F)

alien_viruses <- unique(vir_tab_long_final$virus)
  
virion_m_clean %>% 
  filter(Virus %in% alien_viruses) %>% 
  distinct(Host, Virus) %>% 
  group_by(Virus) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(x = reorder(Virus, -n), y = n))+
  geom_bar(stat = "identity")+
  geom_bar(fill = "#E74F4F", stat = "identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank())+
  labs(y = "Number of hosts in VIRION")

ggsave("m/Alien_viruses.jpeg", width = 6, height = 5, dpi = 600)

dama_countries2 %>%
  group_by(countries) %>%
  summarise("n alien mammals" = n()) -> tab_alien_per_country

Records_of_alien_species %>%
  group_by(countries) %>%
  summarise("n viral records associated to alien mammals" = n()) -> tab_viralrecords_per_country

vir_tab_tot_long %>%
  group_by(countries) %>%
  summarise("n alien viruses" = n()) -> tab_alienviruses_per_country


# Higher taxonomic ranks:
unique(vir_tab_long_final %>% 
  pull(virus)) -> alien_viruses

unique(virion_m_clean %>% 
  filter(Virus %in% alien_viruses) %>% 
  distinct(Host, Virus, VirusFamily) %>% pull(VirusFamily)) -> target_virfams

unique(virion_m_clean %>% 
         filter(Virus %in% alien_viruses) %>% 
         distinct(Host, Virus, VirusGenus) %>% pull(VirusGenus)) -> target_virgenera

virion_m_clean %>% 
  filter(VirusGenus %in% target_virgenera) %>% 
  distinct(Host, Virus, VirusGenus) %>% 
  group_by(VirusGenus) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(x = reorder(VirusGenus, -n), y = n))+
  geom_bar(stat = "identity")+
  geom_bar(fill = "#E74F4F", stat = "identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank())+
  labs(y = "Number of unique associations in VIRION")

ggsave("m/Alien_viruses_genuslevel_associations.jpeg", width = 6, height = 5, dpi = 600)

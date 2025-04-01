#####Required Packages#################
library(dplyr)
library(readr)
library(lubridate)
library(wbstats)
library(stringr)
library(plm)
library(tidyr)
library(ggplot2)
library(broom)
########################################
rm(list = ls())
########################################
# First Step, Policy_Index Estimation ##
########################################

# Load the dataset (Update the file path if necessary)

df <- read.csv(file = "./Document_Data_Download-2025-03-17.csv")

# Extract year from the date column
df <- df %>%
  mutate(First.event.in.timeline = ymd_hms(First.event.in.timeline, tz = "UTC")) %>%  # Convert to Date-Time
  mutate(year = year(First.event.in.timeline))

# A silly check to remove outliers 
df <- subset(df, year > "10")
colnames(df)
# Clean column names (replace spaces with underscores and convert to lowercase)
colnames(df) <- gsub(" ", "_", tolower(colnames(df)))

# Compute quantitative indicators by country
climate_policy_summary <- df %>%
  group_by(geographies, year) %>%
  summarise(
    num_policies = n(),  # Total number of policies
    num_policy_types = n_distinct(document.type),  # Unique policy types
    num_sectors = n_distinct(sector),  # Unique sectors covered
    num_mitigation_policies = sum(grepl("Mitigation", topic.response, ignore.case = TRUE), na.rm = TRUE),
    num_adaptation_policies = sum(grepl("Adaptation", topic.response, ignore.case = TRUE), na.rm = TRUE),
    num_policy_instruments = n_distinct(instrument),  # Unique policy instruments
    num_recent_policies = sum(date.added.to.system > "2020-01-01", na.rm = TRUE)  # Policies added after 2020
  ) %>%
  arrange(desc(num_policies))  # Sort by number of policies

# Min-Max Scaling (Range: 0 to 1)
climate_policy_summary_scaled <- climate_policy_summary %>%
  mutate(across(c(num_policies, num_recent_policies, num_policy_types, 
                  num_policy_instruments, num_sectors, 
                  num_mitigation_policies, num_adaptation_policies), 
                ~ (. - min(.)) / (max(.) - min(.)), 
                .names = "scaled_{.col}"))

# Z-Score Normalization
climate_policy_summary_scaled <- climate_policy_summary %>%
  mutate(across(c(num_policies, num_recent_policies, num_policy_types, 
                  num_policy_instruments, num_sectors, 
                  num_mitigation_policies, num_adaptation_policies), 
                ~ scale(.), .names = "scaled_{.col}"))

# Constructing the Global Policy Index
climate_policy_summary <- climate_policy_summary_scaled %>%
  rowwise() %>%
  mutate(policy_index = 0.175 * scaled_num_policies +
           0.175 * scaled_num_recent_policies +
           0.1375 * scaled_num_policy_types +
           0.15 * scaled_num_policy_instruments +
           0.115 * scaled_num_sectors +
           0.125 * scaled_num_mitigation_policies +
           0.1225 * scaled_num_adaptation_policies) %>%
  ungroup()
#
climate_policy_summary <- subset(climate_policy_summary,
                                 geographies != "No Geography")
#
########################################
# Second Step, Health Data##############
########################################

deaths_data <- wb_data(indicator = c(
  # Death rate attributed to household and ambient air pollution (per 100,000 population)
  "CC.SH.AIRP.AMB",
  
  # Air pollution, mean annual exposure (micrograms per cubic meter)
  "CC.SH.AIRP.AIR",
  
  # PM2.5 air pollution, mean annual exposure (micrograms per cubic meter)
  "EN.ATM.PM25.MC.M3",
  
  # PM2.5 air pollution, population exposed to levels exceeding WHO guideline value (% of total)
  "EN.ATM.PM25.MC.ZS",
  
  # Mortality rate attributed to household and ambient air pollution, age-standardized (per 100,000 population)
  "SH.STA.AIRP.P5",
  
  # Mortality rate attributed to household and ambient air pollution, age-standardized, female (per 100,000 female population)
  "SH.STA.AIRP.FE.P5",
  
  # Mortality rate attributed to household and ambient air pollution, age-standardized, male (per 100,000 male population)
  "SH.STA.AIRP.MA.P5"
), 
start_date = 1980, 
end_date = 2025)
######################################

# Normalize each indicator column using Z-Score
normalized_data <- deaths_data %>%
  mutate(across(
    c(CC.SH.AIRP.AIR, CC.SH.AIRP.AMB, EN.ATM.PM25.MC.M3, 
      EN.ATM.PM25.MC.ZS, SH.STA.AIRP.FE.P5, SH.STA.AIRP.MA.P5, SH.STA.AIRP.P5),
    ~ scale(.),  # Built-in z-score normalization
    .names = "norm_{.col}"
  ))

# Aggregate by country (mean across years)
aggregated_country <- normalized_data %>%
  group_by(country, iso3c) %>%
  summarise(
    across(
      starts_with("norm_"),
      ~ mean(., na.rm = TRUE),
      .names = "agg_{.col}"
    ),
    .groups = "drop"
  )

# Aggregate by year (mean across countries)
aggregated_year <- normalized_data %>%
  group_by(date) %>%
  summarise(
    across(
      starts_with("norm_"),
      ~ mean(., na.rm = TRUE),
      .names = "agg_{.col}"
    ),
    .groups = "drop"
  )

# Create composite index per country-year (mean of normalized indicators)
composite_data <- normalized_data %>%
  rowwise() %>%
  mutate(
    composite_index = mean(c_across(starts_with("norm_")), na.rm = TRUE)
  ) %>%
  ungroup()

######################################
deaths_merged <- composite_data %>%
  # Rename columns to match the merged_data dataset
  rename(
    geographies = country,
    year = date,
    ISO = iso3c
  )
#
########################################
# Third Step, Merge the Datasets########
########################################
#
#
merged_data <- climate_policy_summary %>%
  left_join(deaths_merged, by = c("geographies", "year"))
#
#
########################################
# Forth Step, Model Development#########
########################################
panel_data <- merged_data %>%
  # Keep rows where at least one of the key variables is not NA
  filter(
    !is.na(policy_index) |  # Either policy_index is not NA
      !is.na(composite_index)  # Or ambient air pollution is not NA
  ) %>%
  mutate(
    geographies = as.factor(geographies),
    year = as.numeric(as.character(year))
  )

# Fit the fixed effects panel data model
fixed_effects_model <- plm(
  policy_index ~ composite_index,  # Dependent and independent variables
  data = panel_data,
  model = "within",  # Fixed effects model
  effect = "twoway"  # Allows for both individual and time fixed effects
)

# Summary of the model
summary(fixed_effects_model)


########################################
# Fifth Step, Region-Based Model########
########################################
omid <- read.csv(file = "./final.csv") # Data include the regions;
#this dataset also includes the ghg emission per capita for each region
#
#
merged_data <- subset(merged_data, geographies != "European Union")
merged_data <- merged_data %>%
  mutate(geographies = recode(geographies,
                              "United States of America" = "USA"))
#
merged_data <- merged_data %>%
  left_join(omid %>% select(ISO, region) %>% distinct(), by = "ISO")
###
#merged_data <- merged_data %>%
# filter(!is.na(region))
merged_data <- merged_data %>%
  mutate(region = ifelse(is.na(region), "Far-Lands", region))

####
region_fe_model <- plm(policy_index ~ composite_index + factor(region), 
                       data = merged_data, 
                       index = c("ISO", "year"), 
                       model = "within")
summary(region_fe_model)
#####
interaction_model <- plm(policy_index ~ composite_index * factor(region), 
                         data = merged_data, 
                         index = c("ISO", "year"), 
                         model = "within")
summary(interaction_model)
#####

region_models <- merged_data %>%
  group_by(region) %>%
  do(model = lm(policy_index ~ composite_index, data = .)) %>%
  mutate(tidied = list(tidy(model)))
region_models %>% unnest(tidied)
#####

ggplot(merged_data, aes(x = composite_index, y = policy_index, color = region)) +
  geom_point() + theme_bw() + xlab("Health Index") +
  ylab("Policy Index") +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ region) + 
  theme(
    legend.position = "none",
    legend.background = element_rect(fill = "white", color = "black"),
    legend.text = element_text(size = 12),  # Size of legend text
    legend.title = element_text(size = 14),  # Size of legend title
    axis.title = element_text(size = 14),   # Size of axis titles
    axis.text = element_text(size = 12),    # Size of axis text
    strip.text = element_text(size = 12),   # Size of facet strip text
    plot.title = element_text(size = 16, hjust = 0.5)  # Size of plot title
  ) +  
  labs(title = "Relationship between Composite Index and Policy Index by Region",
       color = "Region") + 
  guides(color = guide_legend(nrow = 3))

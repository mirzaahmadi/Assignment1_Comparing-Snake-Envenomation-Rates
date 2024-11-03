library(BiocManager)
library(remotes)
library(bold)
library(tidyverse)
library(readxl)
library(maps)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(viridis)
library(ggfortify)

# Loading and filtering data ----
raw_external_envenomations_df <- read_excel("Estimation of the total number of snake bite envenomings by country.xlsx")
# Data inspection is done frequently to ensure the structure, type, organization of data set is understood before working with the data
dim(raw_external_envenomations_df)
names(raw_external_envenomations_df)

# Only keep 'Country' and 'Incidence' columns as they are relevant for analysis
external_envenomations_df <- raw_external_envenomations_df %>%
  select(
    Country,
    `Incidence per 100 000 population (Lower)`,
    `Incidence per 100 000 population (Higher)`
  )
names(external_envenomations_df)

# Create average incidence rate from the higher and lower bounds in order to have a single, representative measure of envenoming incidence for each measure
filtered_external_envenomations_df <- external_envenomations_df %>%
  mutate(Average_Incidence_rate_per_100k = (`Incidence per 100 000 population (Lower)` + `Incidence per 100 000 population (Higher)`) / 2) %>%
  select(Country, Average_Incidence_rate_per_100k)
names(filtered_external_envenomations_df)

# squamata_df <- as_tibble(bold_specimens(taxon='squamata'))
squamata_df <- read_tsv("squamata_data.tsv")
dim(squamata_df)
names(squamata_df)

# Subset data to only include relevant venomous snake families for analysis
vipers_elapids_df <- subset(squamata_df, family_name == c("Viperidae", "Elapidae"))
dim(vipers_elapids_df)

# Filter data to include only country and family_name for analyzing distinct family sampling by country
filtered_countries_from_BOLD <- vipers_elapids_df %>%
  select(country, family_name) %>%
  distinct(country) %>%
  mutate(across(where(is.character), ~ na_if(.x, ""))) %>%
  filter(!is.na(country)) # Filter NA values as they are not useful for geographic mapping and visualizations related to geographic distance
names(filtered_countries_from_BOLD)


# Load in venom protein proportion data in order to distinguish differing venom types of different snake families in analysis
unique_species_df <- vipers_elapids_df %>%
  distinct(species_name, .keep_all = TRUE)
dim(unique_species_df)

venom_proportions_df <- read_excel("Venom Protein Proportion Data.xlsx") %>%
  inner_join(unique_species_df, by = c("SPECIES" = "species_name")) %>%
  select(SPECIES, PLA2, SVSP, SVMP, `3FT`, LAAO, CRiSP, `CTL/SNACLEC`, DIS, NP, KUN, VEGF, CYS, DEF, MPi, VT, family_name)
dim(venom_proportions_df)
names(venom_proportions_df)



# Figure 1: Create envenoming incidence map by country ----
merged_country_incidence_df <- filtered_external_envenomations_df %>%
  inner_join(filtered_countries_from_BOLD, by = c("Country" = "country"))

# Merge country lat and long data to create a map of countries in order to visually represent envenoming incidence per country
world_map <- map_data("world")
merged_map_data <- world_map %>%
  left_join(merged_country_incidence_df, by = c("region" = "Country")) %>%
  select(-subregion)

# The envenoming incidence rate is log-transformed because it makes it easier to interpret changes in incidence rates
map <- ggplot(data = merged_map_data, aes(x = long, y = lat, group = group, fill = log(Average_Incidence_rate_per_100k + 1))) +
  geom_polygon(color = "grey") +
  scale_fill_viridis(option = "B", na.value = "grey80") +
  labs(
    title = "Global Snake Bite Envenoming Incidence Rate per 100,000 Population",
    fill = "Log-Transformed Envenoming Incidence Rate"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
    legend.position = c(0.17, 0.4),
    legend.title = element_text(size = 12, face = "plain", hjust = 0, margin = margin(r = 10)),
    legend.text = element_text(size = 12, hjust = 0.5),
    legend.key.width = unit(2, "cm"),
    legend.box.margin = margin(10, 10, 10, 10)
  )

print(map)



# Figure 2: Stacked bar chart of the most prevalent snake families in the top ten countries with the highest envenoming incidence ----
BOLD_sampling_per_country_df <- vipers_elapids_df %>%
  group_by(country) %>%
  summarise(Count = n()) %>%
  na.omit(country)

# Select only the 10 most sampled countries to compare viper and elapid proportions in order to keep the graph clear and uncluttered
most_sampled_countries <- BOLD_sampling_per_country_df %>%
  arrange(desc(Count)) %>%
  head(15) # The top 15 countries are specified, but only 10 appear in the figure because 5 of them are missing from the BOLD viper/elapid data set.


family_proportions_per_country_df <- vipers_elapids_df %>%
  select(country, family_name) %>%
  mutate(across(where(is.character), ~ na_if(.x, ""))) %>%
  filter(!is.na(country)) %>%
  group_by(country) %>%
  summarise(
    Proportion_of_vipers_sampled = sum(family_name == "Viperidae") / (sum(family_name == "Viperidae") + sum(family_name == "Elapidae")),
    Proportion_of_elapids_sampled = sum(family_name == "Elapidae") / (sum(family_name == "Elapidae") + sum(family_name == "Viperidae"))
  ) %>%
  ungroup() %>%
  inner_join(merged_country_incidence_df, by = c("country" = "Country")) %>%
  arrange(desc(Average_Incidence_rate_per_100k)) %>%
  select(-Average_Incidence_rate_per_100k) %>%
  filter(country %in% most_sampled_countries$country) %>%
  mutate(country = factor(country, levels = most_sampled_countries$country))


# Reshape data from wide to long format is done because it makes the data more manageable and ensuring compatibility with analytical tools and visualization libraries, leading to clearer insight
data_long <- family_proportions_per_country_df %>%
  pivot_longer(
    cols = c(Proportion_of_vipers_sampled, Proportion_of_elapids_sampled),
    names_to = "Snake_Family",
    values_to = "Proportion"
  ) %>%
  mutate(country = factor(country, levels = most_sampled_countries$country)) # # Set country as a factor to maintain the same order when displayed on x axis


# Create the stacked bar chart
ggplot(data_long, aes(x = country, y = Proportion, fill = Snake_Family)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = c("#D32F2F", "#1976D2"),
    labels = c("Elapidae", "Viperidae")
  ) + # Set colors for viper and elapid
  labs(
    title = "Proportion of Sampled Vipers and Elapids by Country",
    x = "Top 10 Most Sampled Countries",
    y = "Proportion",
    fill = "Snake Family"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(family = "bold", size = 16),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
    legend.title = element_text(size = 20, face = "plain", hjust = 0, margin = margin(r = 10)),
    legend.text = element_text(size = 16, hjust = 0.5),
  )


# Figure 3: Create PCA plot with venom protein proportion data ----
# Centering is important for PCA as it prevents bias towards variables with larger means, ensuring accurate representation of variance.
venom_protein_proportions_df.pca <- prcomp(venom_proportions_df[, c(2:15)], center = TRUE) # Filter to retain only numeric columns, as PCA requires numeric input.

summary(venom_protein_proportions_df.pca)

pca_scores <- as.data.frame(venom_protein_proportions_df.pca$x)

# ensure the pca_scores data frame has an according species column as each point in the plot will represent a viper or elapid species
pca_scores$SPECIES <- venom_proportions_df$SPECIES

# Different colours are chosen to represent various snake venom types
autoplot(venom_protein_proportions_df.pca,
  data = venom_proportions_df,
  colour = "VT",
  shape = "family_name"
) +
  geom_point(aes(color = venom_proportions_df$VT, shape = venom_proportions_df$family_name), size = 5) +
  scale_color_manual(values = c("green", "cyan", "brown", "red", "blue", "magenta", "black"), name = "Venom Type") +
  scale_shape(name = "Family") +
  ggtitle("Principal Component Analysis of 14 Venom Protein Proportions") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 25),
    axis.title.x = element_text(face = "bold", size = 20),
    axis.title.y = element_text(face = "bold", size = 20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 20, face = "plain"),
    legend.text = element_text(size = 16, hjust = 0),
  )

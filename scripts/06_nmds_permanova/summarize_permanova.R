# permanova_tables --------------------------------------------------------

library(qs)
library(dplyr)
library(tidyr)
library(stringr)
library(glue)

# Set the directory containing the PERMANOVA results
path <- "outputs/permanova_2026"
files <- list.files(path, pattern = "_10percent_allyrs\\.qs$", full.names = TRUE)
# Filter only CanESM2 files
fls <- files[grepl("INM-CM4", files)]


# Function to extract statistics from each file
extract_permanova <- function(file) {
  result <- qread(file)
  result_df <- as.data.frame(result)
  
  # Remove extension and suffix
  filename <- basename(file) %>%
    str_remove("\\.qs$") %>%
    str_remove("_10percent_allyrs$")
  
  # Split parts by "_"
  parts <- str_split(filename, "_")[[1]]
  
  # GCM is always the last element
  gcm <- parts[length(parts)]
  
  # Region is everything before the GCM
  region_name <- paste(parts[1:(length(parts) - 1)], collapse = "_")
  
  cat("Processing region:", region_name, "GCM:", gcm, "\n")
  
  result_clean <- result_df %>%
    mutate(
      Term = rownames(result_df),
      GCM = gcm,
      Region = region_name,
      df = Df,
      SS = round(`Sum Sq`, 2),
      MS = round(`Mean Sq`, 4),
      F = round(`F value`, 2),
      `p-value` = case_when(is.na(`Pr(>F)`) ~ NA_character_,
                            `Pr(>F)` < 0.01 ~ "<0.01",
                            TRUE ~ sprintf("%.3f", `Pr(>F)`)
      )) %>%
    select(GCM, Region, Term, df, SS, MS, F, `p-value`)
  
  return(result_clean)
}

# Process all files
permanova_table <- bind_rows(lapply(fls, extract_permanova))

permanova_table <- permanova_table %>%
  mutate(
    Region = case_when(
      Region == "Low_Arctic_north" ~ "Low Arctic north",
      Region == "Mid_Boreal" ~ "Mid-Boreal",
      TRUE ~ gsub("_", " ", Region)
    )
  )

# View and save
print(permanova_table)
write.csv(permanova_table, glue("./tables/Table4_PERMANOVA_summary _{gcm}.csv"), row.names = FALSE)


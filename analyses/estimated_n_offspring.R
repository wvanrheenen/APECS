#' loading libraries and functions for src folder seems a bit complicated. However, makes sure the script can be run from any location as long as the file/directory structure is maintained as in the github repository.
#' Could be improved...

#' Get the script's path using commandArgs
script_path = normalizePath(sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]))
script_dir  = dirname(script_path)

#' load libraries and functions
source(file.path(script_dir, "/../src/libraries_simPed.R"))
source(file.path(script_dir, "/../src/functions_simPed.R"))

#' locate and read data
fertility_path = file.path(script_dir, "/../data/fertility_rate/gapminder_children_per_woman_total_fertility.csv")
fertility = read.table(fertility_path, header=F, sep=",",  quote="\"", comment.char="", row.names = 1) %>% 
  tibble::rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value) %>% 
  rename("Year" = "country")

mean_fertility_NL = fertility %>%
                      mutate(interval = (Year - min(Year)) %/% 25 + 1) %>%  # Create 25-year interval groups
                      group_by(interval) %>%
                      summarize(
                        start_year = min(Year),
                        end_year = max(Year),
                        mean_fertility_rate = mean(Netherlands, na.rm = TRUE)
                      )

ggplot(fertility, aes(x=Year, y=Netherlands)) +
  geom_line(colour="black") +
  geom_segment(data = mean_fertility_NL, aes(x = start_year, y = mean_fertility_rate, xend = end_year, yend = mean_fertility_rate), colour="red") +
  xlab("Year") +
  ylab("Children per woman (in NL)") +
  theme_classic()


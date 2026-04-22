# Average age of mothers at childbirth - Data package

This data package contains the data that powers the chart ["Average age of mothers at childbirth"](https://ourworldindata.org/grapher/period-average-age-of-mothers?tab=chart&country=SWE~NLD~Europe+%28UN%29&v=1&csvType=filtered&useColumnShortNames=false) on the Our World in Data website. It was downloaded on February 04, 2025.

### Active Filters

A filtered subset of the full data was downloaded. The following filters were applied:
country: SWE, NLD, Europe (UN)
tab: chart

## CSV Structure

The high level structure of the CSV file is that each row is an observation for an entity (usually a country or region) and a timepoint (usually a year).

The first two columns in the CSV file are "Entity" and "Code". "Entity" is the name of the entity (e.g. "United States"). "Code" is the OWID internal entity code that we use if the entity is a country or region. For normal countries, this is the same as the [iso alpha-3](https://en.wikipedia.org/wiki/ISO_3166-1_alpha-3) code of the entity (e.g. "USA") - for non-standard countries like historical countries these are custom codes.

The third column is either "Year" or "Day". If the data is annual, this is "Year" and contains only the year as an integer. If the column is "Day", the column contains a date string in the form "YYYY-MM-DD".

The final column is the data column, which is the time series that powers the chart. If the CSV data is downloaded using the "full data" option, then the column corresponds to the time series below. If the CSV data is downloaded using the "only selected data visible in the chart" option then the data column is transformed depending on the chart type and thus the association with the time series might not be as straightforward.

## Metadata.json structure

The .metadata.json file contains metadata about the data package. The "charts" key contains information to recreate the chart, like the title, subtitle etc.. The "columns" key contains information about each of the columns in the csv, like the unit, timespan covered, citation for the data etc..

## About the data

Our World in Data is almost never the original producer of the data - almost all of the data we use has been compiled by others. If you want to re-use data, it is your responsibility to ensure that you adhere to the sources' license and to credit them correctly. Please note that a single time series may have more than one source - e.g. when we stich together data from different time periods by different producers or when we calculate per capita metrics using population data from a second source.

## Detailed information about the data


## Mean age at childbearing
Mean age of mothers at the birth of their children if women were subject throughout their lives to the age-specific fertility rates observed in a given year.
Last updated: December 18, 2024  
Next update: December 2025  
Date range: 1891–2023  
Unit: years  


### How to cite this data

#### In-line citation
If you have limited space (e.g. in data visualizations), you can use this abbreviated in-line citation:  
Human Fertility Database (2024); UN, World Population Prospects (2024) – processed by Our World in Data

#### Full citation
Human Fertility Database (2024); UN, World Population Prospects (2024) – processed by Our World in Data. “Mean age at childbearing” [dataset]. Human Fertility Database, “Human Fertility Database”; United Nations, “World Population Prospects” [original data].
Source: Human Fertility Database (2024), UN, World Population Prospects (2024) – processed by Our World In Data

### Sources

#### Human Fertility Database
Retrieved on: 2024-11-19  
Retrieved from: https://www.humanfertility.org/Home/Index  

#### United Nations – World Population Prospects
Retrieved on: 2024-07-11  
Retrieved from: https://population.un.org/wpp/Download/  

#### Notes on our processing step for this indicator
This indicator is constructed by combining data from multiple sources:

- Before 1950: Historical estimates by Human Fertility Database (2024).

- 1950-2023: Population records by the UN World Population Prospects (2024 revision).


    
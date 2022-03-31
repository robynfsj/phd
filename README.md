# PhD data analysis and visualisations

This repository contains R scripts for the data analyses and visualisations of 
my PhD results. The data are of diatom species assemblages from Lake Ioannina 
and Lake Ohrid. The data were generated by manually counting diatom valves, 
which were preserved in lake sediment cores, under a microscope. They span a 
time period roughly equivalent to Marine Isotope Stages 7–9 (approximately 
320,000 to 190,000 years ago).


## Scripts

[ioannina_read_data.R](https://github.com/robynfsj/phd/blob/master/scripts/ioannina_read_data.R) 
– reads in the data from an Excel spreadsheet  
  
[ioannina_zonation.R](https://github.com/robynfsj/phd/blob/master/scripts/ioannina_zonation.R) 
– constrained incremental sum of squares (CONISS) cluster analysis to group 
similar adjacent samples together to form diatom assemblage zones  
  
[ioannina_ordination.R](https://github.com/robynfsj/phd/blob/master/scripts/ioannina_ordination.R)
– dimensionality reduction of the multivariate data to pick out the main trends  
  
[ioannina_stratigraphic_diagram.R](https://github.com/robynfsj/phd/blob/master/scripts/ioannina_stratigraphic_diagram.R)
– plots the stratigraphic diagram


## Details of analyses

[Determining the number of diatom valves to count per sample of the Lake Ioannina MIS 7–9 diatom record](https://robynfsj.github.io/phd/count-justification.html) - TODO: review markdown file and publish

[Zonation of the Lake Ioannina diatom data](https://robynfsj.github.io/phd/zonation.html)  

[Ordination of the Lake Ioannina MIS 7–9 diatom data](https://robynfsj.github.io/phd/ordination.html)

[Numerical analyses of the Lake Ohrid MIS 7 diatom data]() - TODO: add Ohrid data and scripts to repo, then review markdown file and publish

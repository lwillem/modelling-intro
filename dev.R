require('rmarkdown')

project_dir <- path.expand(".")
Rmd_dir <- file.path(project_dir, "Rmd")
html_dir <- file.path(project_dir, "docs")

Rmd_files <- c("index","intro_to_r", "intro_to_SIR",
               "intro_to_RStudio",'intro_to_IBM_coding',
               'intro_to_IBM_walk','intro_to_IBM_location',
               'HTA_flu_basic','HTA_flu_CEAresults',
               'HTA_uncertainty') 

for (Rmd_file in Rmd_files) {
  render(file.path(Rmd_dir, sprintf("%s.Rmd", Rmd_file)),
         output_dir = html_dir)
}


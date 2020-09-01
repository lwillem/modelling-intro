require('rmarkdown')

project_dir <- path.expand(".")
Rmd_dir <- file.path(project_dir, "Rmd")
html_dir <- file.path(project_dir, "docs")

Rmd_files <- c("index","intro_to_r_lw", "intro_to_SIR_lw") 

for (Rmd_file in Rmd_files) {
  render(file.path(Rmd_dir, sprintf("%s.Rmd", Rmd_file)),
         output_dir = html_dir)
}



# clear workspace
rm(list = ls())

# load package
require('rmarkdown')

# set directories
project_dir = path.expand(".")
Rmd_dir     = file.path(project_dir, "Rmd")
html_dir    = file.path(project_dir, "docs")

# Rmd files that use show_results to toggle solution visibility
rmd_files = dir(Rmd_dir, full.names = FALSE, pattern = '.Rmd')

# Rmd files that use show_results to toggle solution visibility
rmd_files_results = c("CEA_flu_basic",
                      "CEA_flu_uncertainty",
                      "CEA_RSV_dynamic",
                      "PSA_parameters")

# define help function
knit_with_show_results = function(rmd_name, show_results) {
  e = new.env(parent = globalenv())
  e$show_results = show_results
  output_file = if (show_results) paste0(rmd_name, "_results.html") else paste0(rmd_name, ".html")
  output_file = gsub(".Rmd","",output_file)
  rmd_name = if (grepl(".Rmd",rmd_name)) rmd_name else paste0(rmd_name, ".Rmd")
  render(file.path(Rmd_dir, rmd_name),
         output_file = output_file,
         output_dir  = html_dir,
         envir       = e)
}

# generate all html files
for (rmd_name in rmd_files) {
  knit_with_show_results(rmd_name, show_results = FALSE)
}

# generate specific html files with visible results
for (rmd_name in rmd_files_results) {
  knit_with_show_results(rmd_name, show_results = TRUE)
}

require('rmarkdown')

project_dir = path.expand(".")
Rmd_dir     = file.path(project_dir, "Rmd")
html_dir    = file.path(project_dir, "docs")

# Rmd files that use show_results to toggle solution visibility
rmd_files = c("CEA_flu_basic",
              "CEA_flu_uncertainty",
              "CEA_RSV_dynamic",
              "PSA_parameters")

knit_with_show_results = function(rmd_name, show_results) {
  e = new.env(parent = globalenv())
  e$show_results = show_results
  output_file = if (show_results) paste0(rmd_name, "_results.html") else paste0(rmd_name, ".html")
  render(file.path(Rmd_dir, paste0(rmd_name, ".Rmd")),
         output_file = output_file,
         output_dir  = html_dir,
         envir       = e)
}

for (rmd_name in rmd_files) {
  knit_with_show_results(rmd_name, show_results = TRUE)
  knit_with_show_results(rmd_name, show_results = FALSE)
}

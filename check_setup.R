##############################################################################
# Setup check for the Health Economics course
#
# Run this script before the course to verify that your R version and all
# required packages are correctly installed.
##############################################################################

# ---- R version check --------------------------------------------------------

r_version_required = "4.1.0"

if (getRversion() < r_version_required) {
  stop(paste("R version", r_version_required, "or higher is required.",
             "You are running", getRversion(),
             "- please update R at https://cran.r-project.org"))
} else {
  message("R version OK: ", getRversion())
}

# ---- Package installation and check -----------------------------------------

# these packages are really required
required_packages = c(
  "mgcv",         # gam-function for EVPPI 
  "EpiDynamics"   # epidemic dynamics utilities 
)

# if this package cannot be installed, we provide a workaround
optional_packages = c(
  'BCEA' # for cost-effectiveness calculations with uncertainty
)

install_if_missing = function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing: ", pkg)
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

invisible(lapply(required_packages, install_if_missing))

invisible(lapply(optional_packages, install_if_missing))

message("Check passed.")

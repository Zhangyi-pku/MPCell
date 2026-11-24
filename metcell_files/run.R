# ============================================================
# run.R — MetCell pipeline, 所有参数从环境变量读取
# ============================================================
library(MetCell)

.env_str  <- function(name, default) {
  v <- Sys.getenv(name, unset = "")
  if (nchar(v) == 0) default else v
}
.env_num  <- function(name, default) as.numeric(.env_str(name, as.character(default)))
.env_int  <- function(name, default) as.integer(.env_str(name, as.character(default)))
.env_bool <- function(name, default) {
  v <- .env_str(name, as.character(default))
  tolower(v) %in% c("true", "t", "1", "yes")
}

# ---------- 路径 ----------
d_path           <- .env_str("METCELL_D_PATH", "")
if (d_path == "") stop("METCELL_D_PATH is required")
output_dir       <- .env_str("METCELL_OUTPUT_DIR", "/data/output")
script_dir       <- .env_str("METCELL_SCRIPT_DIR", "/data/scripts")
time_limit_table <- .env_str("METCELL_TIME_LIMIT_TABLE",
                             file.path(script_dir, "time_limit_table.csv"))
bruker_lib       <- .env_str("METCELL_BRUKER_LIB",
                             file.path(script_dir, "libtimsdata.so"))

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
wd <- output_dir
setwd(wd)

cat("===== MetCell run.R =====\n")
cat("d_path           :", d_path, "\n")
cat("output_dir       :", output_dir, "\n")
cat("time_limit_table :", time_limit_table, "\n")
cat("bruker_lib       :", bruker_lib, "\n")

# ---------- 参数 ----------
n_slaves      <- .env_int("METCELL_N_SLAVES", 2)
ion_mode      <- .env_str("METCELL_ION_MODE", "positive")
query_pool    <- .env_int("METCELL_QUERY_POOL_SIZE", 2000)
marker_mz     <- .env_num("METCELL_MARKER_MZ", 760.5845)
marker_mob    <- .env_num("METCELL_MARKER_MOBILITY", 1.402)
mz_tol_ppm    <- .env_num("METCELL_MZ_TOL_PPM", 20)
mob_range     <- .env_num("METCELL_MOBILITY_RANGE", 0.05)
int_upper     <- .env_num("METCELL_INT_THRES_UPPER", 600000)
int_lower     <- .env_num("METCELL_INT_THRES_LOWER", 60000)
marker_span   <- .env_int("METCELL_MARKER_EIC_PEAK_SPAN", 5)
super_mz_tol  <- .env_num("METCELL_SUPERPOSITION_MZ_TOL", 20)
cell_number   <- .env_int("METCELL_CELL_NUMBER", 100)
super_pool    <- .env_int("METCELL_SUPERPOSITION_POOL", 40)
search_mz_tol <- .env_num("METCELL_SEARCH_PEAK_MZ_TOL", 20)
search_n_skip <- .env_int("METCELL_SEARCH_N_SKIP", 1)
peak_len      <- .env_int("METCELL_PEAK_TARGET_LENGTH", 15)
smooth_window <- .env_int("METCELL_SMOOTH_WINDOW", 10)
span_det      <- .env_int("METCELL_PEAK_SPAN_DETECTION", 21)
span_intg     <- .env_int("METCELL_PEAK_SPAN_INTEGRATION", 27)
signal_sd     <- .env_num("METCELL_SIGNAL_SD_THRESHOLD", 0.35)
chg_slope     <- .env_num("METCELL_SINGLE_CHARGE_SLOPE", 0.0009)
chg_intercept <- .env_num("METCELL_SINGLE_CHARGE_INTERCEPT", 0.35)
extract_mz    <- .env_num("METCELL_EXTRACT_MZ_TOL", 20)
extract_mob   <- .env_num("METCELL_MOBILITY_TOL", 0.04)
min_cell_num  <- .env_int("METCELL_MIN_CELL_NUMBER", 10)
extract_pool  <- .env_int("METCELL_EXTRACT_POOL", 150)
iso_mz_tol    <- .env_num("METCELL_ISO_MZ_TOL", 20)
iso_mob_tol   <- .env_num("METCELL_ISO_MOBILITY_TOL", 0.01)
iso_delta     <- .env_num("METCELL_ISO_DELTA", 1.003355)
iso_max_num   <- .env_int("METCELL_ISO_MAX_NUM", 4)
iso_ratio_chk <- .env_bool("METCELL_ISO_INT_RATIO_CHECK", TRUE)
iso_ratio_cut <- .env_num("METCELL_ISO_INT_CUTOFF", 500)
sc_type_ccs   <- .env_str("METCELL_SC_TYPE_CCS", "percentage")
sc_tol_ccs    <- .env_num("METCELL_SC_TOL_CCS", 4)
sc_tol_mz     <- .env_num("METCELL_SC_TOL_MZ", 20)
sc_adduct     <- .env_str("METCELL_SC_ADDUCT", "[M+H]+")
eic_mz_tol    <- .env_num("METCELL_EIC_MZ_TOL", 20)
eic_chunk     <- .env_int("METCELL_EIC_CHUNK_SIZE", 500)
eic_n_threads <- .env_int("METCELL_EIC_N_THREADS", 4)

#### Experiment setup ##########################################################
exp <- Experiment(wd = wd, nSlaves = n_slaves, ion_mode = ion_mode,
                  experiment_type = 'SingleCell')
tims_data <- TimsData(exp)

#### Single-cell frame detection ###############################################
param <- QueryTimsDataSegmentParam(
  time_limit_table = time_limit_table,
  pool_size = query_pool, rerun = FALSE)
sc_param <- DiscoverSCEventsParam(
  marker_mz = marker_mz, marker_mobility = marker_mob,
  mz_tolerance_ppm = mz_tol_ppm, mobility_range = mob_range,
  intensity_abs_threshold_upper = int_upper,
  intensity_abs_threshold_lower = int_lower,
  marker_eic_peak_span = marker_span)
tims_data <- QueryTimsDataSegment(tims_data, param, sc_param)

#### Cell superposition ########################################################
param <- UniteCellSuperpositionFrameParam(
  mz_tolerance_combine = super_mz_tol, cell_number = cell_number,
  pool_size = super_pool, rerun = FALSE)
tims_data <- UniteCellSuperpositionFrame(tims_data, param)

#### Peak detection ############################################################
param <- SearchPeakTargetParam(mz_tol_ppm = search_mz_tol, n_skip = search_n_skip,
                               peak_target_length = peak_len, rerun = FALSE)
tims_data <- SearchPeakTargets(tims_data, param)

param <- DetectEIMPeaksParam(smooth_window = smooth_window,
                             peak_span_eim_detection = span_det,
                             peak_span_eim_integration = span_intg,
                             signal_sd_threshold = signal_sd,
                             skip_invalid_peaks = TRUE,
                             single_charge_line_slope = chg_slope,
                             single_charge_line_intercept = chg_intercept,
                             rerun = FALSE)
tims_data <- DetectEIMPeaks(tims_data, param)

#### Peak quantification #######################################################
param <- ExtractIndividualCellParam(mz_tol_ppm = extract_mz,
                                    mobility_tol = extract_mob,
                                    minimal_cell_number = min_cell_num,
                                    pool_size = extract_pool, rerun = FALSE)
tims_data <- ExtractIndividualCell(tims_data, param)

#### Metabolite annotation #####################################################
iso_param <- IsotopeParam(mz_tol_ppm = iso_mz_tol, mobility_tol = iso_mob_tol,
                          isotope_delta = iso_delta, isotope_max_num = iso_max_num,
                          isotope_int_ratio_check = iso_ratio_chk,
                          isotope_int_ratio_cutoff = iso_ratio_cut, rerun = FALSE)
param <- SCMatchParam(typeCCS = sc_type_ccs, toleranceCCS = sc_tol_ccs,
                      tolerancemz = sc_tol_mz, adduct_for_id = sc_adduct,
                      rerun = FALSE)
tims_data <- IdentifySCPeaks(tims_data, iso_param, param)

#### Full-frame EIC extraction #################################################
source(file.path(script_dir, "extract_eic_full.R"))

ft <- read.csv(file.path(wd, "results/01_feature_table.csv"), stringsAsFactors = FALSE)
cat(sprintf("\n===== Extracting EIC for ALL %d features =====\n", nrow(ft)))

res <- tryCatch({
  extract_feature_eic_full(
    tims_data,
    d_path             = d_path,
    feature_table_path = file.path(wd, "results/01_feature_table.csv"),
    output_prefix      = file.path(wd, "results/all_features_EIC_full"),
    mz_tol_ppm         = eic_mz_tol,
    chunk_size         = eic_chunk,
    n_threads          = eic_n_threads,
    feature_names      = NULL,
    output_format      = "npy",
    bruker_lib         = bruker_lib,
    verbose            = TRUE
  )
}, error = function(e) { cat("FULL EIC FAILED:", conditionMessage(e), "\n"); NULL })

if (!is.null(res)) {
  cat("\n===== FULL EIC DONE =====\n")
  cat("intensity dim :", paste(dim(res$intensity), collapse = " x "), "\n")
  cat("RT range (s)  :",
      sprintf("%.2f - %.2f", min(res$retention_time, na.rm=TRUE),
              max(res$retention_time, na.rm=TRUE)), "\n")
}

quit(save = "no", status = 0)
# =============================================================
# extract_eic_full.R
# Full-frame EIC extraction: reads ALL MS1 frames directly from
# the raw .d file, chunk by chunk, and aggregates intensity per
# candidate peak. Output is a dense (n_features x n_frames) npy
# matrix, plus retention-time vector and feature metadata.
#
# 依赖：opentimsr + 本地 Bruker 转换库 libtimsdata.so
#   （默认容器内路径 /data/bruker_lib/libtimsdata.so，可用 bruker_lib 参数覆盖）
#
# Usage in run.R (at the very end, AFTER the full pipeline):
#   source("extract_eic_full.R")
#   extract_feature_eic_full(tims_data, mz_tol_ppm = 20)
# =============================================================

# ---- Pure-R NPY writer (no external packages) ----
.write_npy <- function(x, file, dtype = NULL) {
  if (is.null(dtype)) {
    dtype <- if (is.double(x)) "<f8" else if (is.integer(x)) "<i4" else
      stop("Unsupported type for NPY writer: ", class(x))
  }
  if (is.null(dim(x))) {
    shape <- paste0("(", length(x), ",)")
    arr <- x
  } else {
    shape <- paste0("(", paste(dim(x), collapse = ", "), ")")
    arr <- as.vector(t(x))
  }
  header_base <- sprintf(
    "{'descr': '%s', 'fortran_order': False, 'shape': %s, }",
    dtype, shape
  )
  total_wo_pad <- 10L + nchar(header_base, type = "bytes") + 1L
  pad <- (64L - (total_wo_pad %% 64L)) %% 64L
  header <- paste0(header_base, strrep(" ", pad), "\n")
  hdr_len <- nchar(header, type = "bytes")
  
  magic <- as.raw(c(0x93, 0x4E, 0x55, 0x4D, 0x50, 0x59))
  ver   <- as.raw(c(1, 0))
  con <- file(file, "wb")
  on.exit(close(con), add = TRUE)
  writeBin(magic, con)
  writeBin(ver,   con)
  writeBin(as.integer(hdr_len), con, size = 2, endian = "little")
  writeBin(charToRaw(header), con)
  if (dtype == "<f8") {
    writeBin(as.double(arr), con, size = 8, endian = "little")
  } else if (dtype == "<i4") {
    writeBin(as.integer(arr), con, size = 4, endian = "little")
  } else {
    stop("dtype not implemented: ", dtype)
  }
  invisible(file)
}

# ---- Robust column picker for opentimsr query results ----
.pick_col <- function(df, candidates) {
  for (cc in candidates) if (cc %in% colnames(df)) return(df[[cc]])
  NULL
}

# ---- 激活 tof -> mz 转换器（必须在 OpenTIMS 之前调用） ----
.activate_converter <- function(bruker_lib = "/data/bruker_lib/libtimsdata.so",
                                verbose = TRUE) {
  ok <- FALSE
  
  # 策略 1: 本地 Bruker 专有库（已验证可用）
  if (!is.null(bruker_lib) && nzchar(bruker_lib) && file.exists(bruker_lib)) {
    tryCatch({
      opentimsr::setup_bruker_so(bruker_lib)
      ok <- TRUE
      if (verbose) cat("[eic_full] Converter: setup_bruker_so(", bruker_lib, ") OK\n")
    }, error = function(e) {
      if (verbose) cat("[eic_full] setup_bruker_so() failed:", conditionMessage(e), "\n")
    })
  } else if (verbose) {
    cat("[eic_full] Bruker lib not found at:", bruker_lib, "\n")
  }
  
  # 策略 2: opentimsr 命名空间内的 setup_opensource（旧版本可能没有）
  if (!ok) {
    ns <- asNamespace("opentimsr")
    if (exists("setup_opensource", envir = ns, inherits = FALSE)) {
      tryCatch({
        get("setup_opensource", envir = ns)()
        ok <- TRUE
        if (verbose) cat("[eic_full] Converter: opentimsr namespace setup_opensource() OK\n")
      }, error = function(e) {
        if (verbose) cat("[eic_full] setup_opensource() via namespace failed:",
                         conditionMessage(e), "\n")
      })
    }
  }
  
  # 策略 3: MsBackendTimsTof 包装器（若安装）
  if (!ok && requireNamespace("MsBackendTimsTof", quietly = TRUE)) {
    tryCatch({
      MsBackendTimsTof::setup_converter_library(opensource = TRUE)
      ok <- TRUE
      if (verbose) cat("[eic_full] Converter: MsBackendTimsTof::setup_converter_library() OK\n")
    }, error = function(e) {
      if (verbose) cat("[eic_full] MsBackendTimsTof failed:", conditionMessage(e), "\n")
    })
  }
  
  if (!ok && verbose) {
    cat("[eic_full] WARNING: no tof->mz converter activated; query(mz) will return zeros.\n")
  }
  invisible(ok)
}

# ---- 把 query 结果统一整理成 (frames_list, used_ids) ----
# 返回 list(frames_list, used_ids)，frames_list 中每个元素是一帧的 data.frame
.normalize_query_result <- function(res, chunk_ids, verbose_diag = FALSE) {
  frame_candidates <- c("frame", "Frame", "FRAME", "frame_id", "FrameId", "Id")
  
  if (is.data.frame(res)) {
    # 合并的 data.frame：必须有 frame 列
    fc <- .pick_col(res, frame_candidates)
    if (is.null(fc)) {
      stop("query() returned data.frame without a 'frame' column; ",
           "cannot assign points to frames. Columns returned: ",
           paste(colnames(res), collapse = ", "))
    }
    sp <- split(res, fc)
    used_ids <- as.character(names(sp))
    if (verbose_diag) {
      cat("[DIAG5] query() returned merged data.frame; split into",
          length(sp), "frames by 'frame' column.\n")
    }
    return(list(frames_list = sp, used_ids = used_ids))
  }
  
  if (is.list(res)) {
    # list 有两种情况：
    #  (a) 长度为 1，且元素里有 frame 列 -> 是合并表
    #  (b) 每个元素是一帧的 data.frame
    if (length(res) == 1 && is.data.frame(res[[1]])) {
      fc <- .pick_col(res[[1]], frame_candidates)
      if (!is.null(fc)) {
        sp <- split(res[[1]], fc)
        used_ids <- as.character(names(sp))
        if (verbose_diag) {
          cat("[DIAG5] query() returned list of length 1 with 'frame' column; split into",
              length(sp), "frames.\n")
        }
        return(list(frames_list = sp, used_ids = used_ids))
      }
    }
    
    frames_list <- res
    used_ids <- vapply(seq_along(frames_list), function(i) {
      f <- frames_list[[i]]
      if (is.data.frame(f)) {
        fc <- .pick_col(f, frame_candidates)
        if (!is.null(fc) && length(fc) > 0) return(as.character(fc[1]))
      }
      chunk_ids[i]
    }, character(1))
    if (verbose_diag) {
      cat("[DIAG5] query() returned list of", length(frames_list),
          "frames; used per-element frame column when available.\n")
    }
    return(list(frames_list = frames_list, used_ids = used_ids))
  }
  
  stop("Unexpected query result type: ", class(res)[1])
}

extract_feature_eic_full <- function(
    tims_data,
    d_path             = NULL,
    feature_table_path = "results/01_feature_table.csv",
    output_prefix      = "results/all_features_EIC_full",
    mz_tol_ppm         = 20,
    feature_names      = NULL,
    chunk_size         = 500L,
    n_threads          = 4L,
    bruker_lib         = "/data/bruker_lib/libtimsdata.so",
    output_format      = c("npy", "both"),
    verbose            = TRUE
) {
  output_format <- match.arg(output_format)
  
  # ---- 1. Candidate peak table ----
  if (!file.exists(feature_table_path)) {
    stop("Feature table not found: ", feature_table_path)
  }
  feature_table <- read.csv(feature_table_path, stringsAsFactors = FALSE)
  if (!"mz" %in% colnames(feature_table)) {
    stop("Column 'mz' is missing from feature_table.")
  }
  if (!is.null(feature_names)) {
    if (!"name" %in% colnames(feature_table))
      stop("Column 'name' is missing from feature_table.")
    feature_table <- feature_table[feature_table$name %in% feature_names, , drop = FALSE]
    if (nrow(feature_table) == 0) stop("No peaks remain after feature_names filter.")
  }
  
  feature_mz   <- as.numeric(feature_table$mz)
  feature_tol  <- feature_mz * mz_tol_ppm * 1e-6
  feature_name <- if ("name" %in% colnames(feature_table))
    as.character(feature_table$name)
  else paste0("feature_", seq_along(feature_mz))
  n_features <- length(feature_mz)
  if (verbose) cat(sprintf("[eic_full] Features: %d\n", n_features))
  
  # ---- 2. Resolve raw .d path ----
  if (is.null(d_path)) {
    d_path <- tryCatch(tims_data@files[1], error = function(e) NULL)
  }
  if (is.null(d_path) || !dir.exists(d_path)) {
    stop("Raw .d folder not found. Pass d_path explicitly.\n",
         "  tried: ", if (is.null(d_path)) "<NULL>" else d_path)
  }
  if (verbose) cat(sprintf("[eic_full] Raw .d: %s\n", d_path))
  
  # ---- 3. OpenTIMS（必须先激活转换器） ----
  if (verbose) cat("[eic_full] Opening .d:", d_path, "\n")
  tryCatch(opentimsr::opentims_set_threads(as.integer(n_threads)), error = function(e) NULL)
  
  # 关键：转换器必须在 OpenTIMS 之前设定
  .activate_converter(bruker_lib = bruker_lib, verbose = verbose)
  
  D <- opentimsr::OpenTIMS(d_path)
  on.exit({ tryCatch(opentimsr::CloseTIMS(D), error = function(e) NULL) }, add = TRUE)
  if (verbose) cat("[eic_full] .d opened\n")
  
  # ---- DIAG1: raw .d 帧数 ----
  cat("[DIAG1] D@frames nrow:", nrow(D@frames), "\n")
  cat("[DIAG1] MS1(D) length:", length(opentimsr::MS1(D)), "\n")
  cat("[DIAG1] retention_times(D) length:", length(opentimsr::retention_times(D)), "\n")
  
  # ---- 4. Frame list ----
  ns <- asNamespace("opentimsr")
  get_fun <- function(nm) if (exists(nm, envir = ns, inherits = FALSE)) get(nm, envir = ns) else NULL
  
  all_frame_ids <- NULL; all_frame_times <- NULL
  
  # 优先：直接用 D@frames
  if (!is.null(D@frames) && is.data.frame(D@frames) && nrow(D@frames) > 0) {
    idcol <- intersect(c("Id","id","Frame","frame"), colnames(D@frames))
    tcol  <- intersect(c("Time","time","RetentionTime"), colnames(D@frames))
    idcol <- if (length(idcol)) idcol[1] else 1
    tcol  <- if (length(tcol))  tcol[1]  else 2
    all_frame_ids   <- as.integer(D@frames[[idcol]])
    all_frame_times <- as.numeric(D@frames[[tcol]])
    if (verbose) cat(sprintf("[eic_full] from D@frames: %d frames, cols %s / %s\n",
                             length(all_frame_ids), idcol, tcol))
  }
  
  # 回退 1: frames() / retention_times()
  if (is.null(all_frame_ids) || is.null(all_frame_times)) {
    frames_fun <- get_fun("frames")
    if (!is.null(frames_fun)) {
      val <- tryCatch(frames_fun(D), error = function(e) {
        if (verbose) cat("[eic_full] frames() FAILED |", class(e)[1], "|", conditionMessage(e), "\n"); NULL })
      if (!is.null(val)) {
        if (is.data.frame(val)) { all_frame_ids <- as.integer(val[[1]])
        } else if (is.atomic(val)) all_frame_ids <- as.integer(val)
        if (verbose && !is.null(all_frame_ids)) cat("[eic_full] frame IDs:", length(all_frame_ids), "\n")
      }
    }
    
    rt_fun <- get_fun("retention_times")
    if (!is.null(rt_fun)) {
      val <- tryCatch(rt_fun(D), error = function(e) {
        if (verbose) cat("[eic_full] retention_times() FAILED |", class(e)[1], "|", conditionMessage(e), "\n"); NULL })
      if (!is.null(val) && is.atomic(val)) {
        all_frame_times <- as.numeric(val)
        if (verbose) cat("[eic_full] frame times:", length(all_frame_times), "\n")
      }
    }
  }
  
  # 回退 2: table2df('Frames')
  if (is.null(all_frame_ids) || is.null(all_frame_times)) {
    if (verbose) cat("[eic_full] Fallback to table2df('Frames')\n")
    fd <- tryCatch(opentimsr::table2df(D, "Frames"), error = function(e) NULL)
    if (!is.null(fd)) {
      if (is.list(fd) && length(fd) >= 1 && !is.data.frame(fd)) fd <- fd[[1]]
      if (is.data.frame(fd) && nrow(fd) > 0 && ncol(fd) >= 2) {
        idcol <- intersect(c("Id","id","Frame","frame"), colnames(fd))
        tcol  <- intersect(c("Time","time","RetentionTime"), colnames(fd))
        idcol <- if (length(idcol)) idcol[1] else 1
        tcol  <- if (length(tcol))  tcol[1]  else 2
        if (is.null(all_frame_ids))   all_frame_ids   <- as.integer(fd[[idcol]])
        if (is.null(all_frame_times)) all_frame_times <- as.numeric(fd[[tcol]])
        if (verbose) cat("[eic_full] used cols:", idcol, "|", tcol, "\n")
      }
    }
  }
  
  if (is.null(all_frame_ids) || is.null(all_frame_times))
    stop("Cannot get frame list; see [eic_full] output above for details.")
  if (length(all_frame_ids) != length(all_frame_times))
    stop("frame length mismatch: ", length(all_frame_ids), " vs ", length(all_frame_times))
  
  # ---- DIAG2 ----
  cat("[DIAG2] all_frame_ids length:", length(all_frame_ids), "\n")
  cat("[DIAG2] all_frame_times length:", length(all_frame_times), "\n")
  cat("[DIAG2] head(all_frame_ids):", paste(head(all_frame_ids, 5), collapse = ", "), "\n")
  
  # ---- MS1 过滤 ----
  ms1_ids <- tryCatch(opentimsr::MS1(D), error = function(e) NULL)
  if (!is.null(ms1_ids) && length(ms1_ids) > 0) {
    keep <- all_frame_ids %in% ms1_ids
    ms1_mask <- keep
    if (verbose) cat(sprintf("[eic_full] MS1(D) keep: %d / %d\n", sum(keep), length(keep)))
  } else {
    ms1_mask <- rep(TRUE, length(all_frame_ids))
    fd2 <- tryCatch(opentimsr::table2df(D, "Frames"), error = function(e) NULL)
    if (!is.null(fd2)) {
      if (is.list(fd2) && length(fd2) >= 1 && !is.data.frame(fd2)) fd2 <- fd2[[1]]
      if (is.data.frame(fd2) && nrow(fd2) == length(all_frame_ids)) {
        msms_col <- intersect(c("MsMsType","msms_type","MsMsType"), colnames(fd2))
        if (length(msms_col)) {
          ms1_mask <- fd2[[msms_col[1]]] == 0
          if (verbose) cat(sprintf("[eic_full] %s==0: %d MS1 frames\n", msms_col[1], sum(ms1_mask)))
        }
      }
    }
  }
  
  frame_ids   <- as.character(all_frame_ids[ms1_mask])
  frame_times <- as.numeric(all_frame_times[ms1_mask])
  n_frames    <- length(frame_ids)
  if (n_frames == 0) stop("No MS1 frames after filtering.")
  
  # ---- DIAG3 ----
  cat("[DIAG3] ms1_mask length:", length(ms1_mask), "\n")
  cat("[DIAG3] sum(ms1_mask):", sum(ms1_mask), "\n")
  cat("[DIAG3] frame_ids length:", length(frame_ids), "\n")
  cat("[DIAG3] frame_times length:", length(frame_times), "\n")
  
  if (verbose) cat(sprintf("[eic_full] MS1 frames: %d\n", n_frames))
  if (verbose) cat(sprintf("[eic_full] Result matrix: %d x %d = %.2f GB (float64)\n",
                           n_features, n_frames, n_features * n_frames * 8 / 1024^3))
  
  # ---- 5. Preallocate ----
  intensity_mat <- matrix(0, nrow = n_features, ncol = n_frames)
  
  # ---- 6. Chunked read ----
  chunk_size <- max(1L, as.integer(chunk_size))
  starts <- seq(1L, n_frames, by = chunk_size)
  t0 <- Sys.time()
  diag_done <- FALSE
  
  for (start in starts) {
    end <- min(start + chunk_size - 1L, n_frames)
    chunk_ids <- frame_ids[start:end]
    chunk_ids_int <- as.integer(chunk_ids)
    
    # 关键修复：显式请求 frame 列，保证每点能归属到正确的帧
    res <- tryCatch(
      opentimsr::query(D, frames = chunk_ids_int,
                       columns = c("frame", "scan", "mz", "intensity")),
      error = function(e) {
        # 回退：某些版本可能不接受 "frame" 列名，去掉再试
        tryCatch(
          opentimsr::query(D, frames = chunk_ids_int,
                           columns = c("scan", "mz", "intensity")),
          error = function(e2) {
            warning("query failed: ", conditionMessage(e2))
            NULL
          }
        )
      }
    )
    if (is.null(res)) {
      if (!diag_done && verbose) cat("[DIAG5] query returned NULL (first chunk)\n")
      next
    }
    
    # 统一整理成 frames_list + used_ids
    norm <- tryCatch(
      .normalize_query_result(res, chunk_ids, verbose_diag = (!diag_done && verbose)),
      error = function(e) {
        warning("normalize_query_result failed: ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(norm)) next
    frames_list <- norm$frames_list
    used_ids    <- norm$used_ids
    
    # ---- DIAG5：第一个 chunk 的结构（详细） ----
    if (!diag_done && verbose) {
      cat("[DIAG5] res class:", class(res)[1],
          "| n frames_list:", length(frames_list), "\n")
      cat("[DIAG5] chunk_ids length:", length(chunk_ids), "\n")
      cat("[DIAG5] used_ids head:", paste(head(used_ids, 5), collapse = ", "), "\n")
      f1 <- frames_list[[1]]
      if (is.data.frame(f1)) {
        cat("[DIAG5] first frame cols:", paste(colnames(f1), collapse = ", "), "\n")
        cat("[DIAG5] first frame nrow:", nrow(f1), "\n")
        mzv  <- .pick_col(f1, c("mz","Mz","m/z"))
        intv <- .pick_col(f1, c("intensity","Intensity","Intensities"))
        if (!is.null(mzv))  cat("[DIAG5] mz range:", paste(range(mzv, na.rm = TRUE), collapse = " - "), "\n")
        if (!is.null(intv)) cat("[DIAG5] intensity range:", paste(range(intv, na.rm = TRUE), collapse = " - "), "\n")
      }
    }
    
    mz_list  <- lapply(frames_list, function(f) {
      if (is.null(f) || !is.data.frame(f)) return(numeric(0))
      v <- .pick_col(f, c("mz","Mz","m/z")); if (is.null(v)) numeric(0) else as.numeric(v)
    })
    int_list <- lapply(frames_list, function(f) {
      if (is.null(f) || !is.data.frame(f)) return(numeric(0))
      v <- .pick_col(f, c("intensity","Intensity","Intensities")); if (is.null(v)) numeric(0) else as.numeric(v)
    })
    n_per <- sapply(mz_list, length)
    
    all_mz    <- unlist(mz_list,  use.names = FALSE)
    all_int   <- unlist(int_list, use.names = FALSE)
    all_frame <- rep(used_ids, n_per)
    
    ok <- is.finite(all_mz) & is.finite(all_int) & all_int > 0
    all_mz <- all_mz[ok]; all_int <- all_int[ok]; all_frame <- all_frame[ok]
    
    if (length(all_mz) > 0) {
      ord <- order(all_mz)
      all_mz <- all_mz[ord]; all_int <- all_int[ord]; all_frame <- all_frame[ord]
      
      for (i in seq_len(n_features)) {
        mz_i  <- feature_mz[i]
        tol_i <- feature_tol[i]
        lo <- findInterval(mz_i - tol_i - 1e-12, all_mz)
        hi <- findInterval(mz_i + tol_i,        all_mz)
        
        # ---- DIAG4：第一个 chunk 每个特征的匹配点数 ----
        if (!diag_done && verbose) {
          cat(sprintf("[DIAG4] feature %s mz=%.5f tol=%.5f -> %d raw points\n",
                      feature_name[i], mz_i, tol_i, max(0, hi - lo)))
        }
        
        if (hi <= lo) next
        idx <- seq.int(lo + 1L, hi)
        agg <- tapply(all_int[idx], all_frame[idx], sum)
        pos <- match(names(agg), chunk_ids) + (start - 1L)
        okp <- !is.na(pos)
        if (any(okp)) intensity_mat[i, pos[okp]] <- as.numeric(agg[okp])
      }
    } else if (!diag_done && verbose) {
      cat("[DIAG4] first chunk: all_mz length = 0 (no mz extracted)\n")
    }
    
    if (start == starts[1]) diag_done <- TRUE
    
    if (verbose) {
      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      eta <- elapsed / (end / n_frames) - elapsed
      cat(sprintf("  frames %d-%d / %d  |  elapsed %.1fs  ETA %.1fs\n",
                  start, end, n_frames, elapsed, max(0, eta)))
    }
  }
  
  # ---- 7. Save ----
  out_dir <- dirname(output_prefix)
  if (out_dir != "." && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  feat_meta <- data.frame(
    feature_name = feature_name,
    mz           = feature_mz,
    mz_tol_da    = feature_tol,
    stringsAsFactors = FALSE
  )
  
  .write_npy(intensity_mat, paste0(output_prefix, "_intensity.npy"), dtype = "<f8")
  .write_npy(frame_times,   paste0(output_prefix, "_retention_time.npy"), dtype = "<f8")
  write.csv(feat_meta, paste0(output_prefix, "_features.csv"), row.names = FALSE)
  
  if (verbose) {
    cat(sprintf("[eic_full] Saved:\n  %s_intensity.npy       shape=(%d, %d)\n  %s_retention_time.npy  shape=(%d,)\n  %s_features.csv\n",
                output_prefix, n_features, n_frames,
                output_prefix, n_frames, output_prefix))
  }
  
  if (output_format == "both") {
    if (verbose) cat("[eic_full] Writing sparse CSV (non-zero only)...\n")
    nz_idx <- which(intensity_mat > 0, arr.ind = TRUE)
    if (nrow(nz_idx) > 0) {
      sparse_df <- data.frame(
        feature_name       = feature_name[nz_idx[, 1]],
        mz                 = feature_mz[nz_idx[, 1]],
        frame              = frame_ids[nz_idx[, 2]],
        retention_time_sec = frame_times[nz_idx[, 2]],
        intensity          = intensity_mat[nz_idx],
        stringsAsFactors   = FALSE
      )
      sparse_df <- sparse_df[order(sparse_df$feature_name, sparse_df$retention_time_sec), ]
      write.csv(sparse_df, paste0(output_prefix, "_sparse.csv"), row.names = FALSE)
      if (verbose) cat(sprintf("[eic_full] Sparse CSV rows: %d\n", nrow(sparse_df)))
    } else {
      if (verbose) cat("[eic_full] No non-zero entries; sparse CSV not written.\n")
    }
  }
  
  invisible(list(
    intensity      = intensity_mat,
    retention_time = frame_times,
    frame_ids      = frame_ids,
    features       = feat_meta,
    output_prefix  = output_prefix
  ))
}
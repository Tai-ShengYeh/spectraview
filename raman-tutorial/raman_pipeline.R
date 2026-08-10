#!/usr/bin/env Rscript
# =====================================================================
# raman_pipeline.R — 顏料拉曼光譜分析流程（教學版，R 語言）
# ---------------------------------------------------------------------
# 與 raman_pipeline.py 功能對等：
#   讀檔 → 扣基線(ALS) → 平滑(SG) → 找峰 → 比對譜庫 → NNLS 解混 → 出圖出表
#
# 只用 base R + Matrix（Matrix 是 R 官方隨附的 recommended 套件，
# 不必額外安裝），所以在任何一台裝了 R 的電腦上都能直接跑。
#
# 用法：
#   Rscript raman_pipeline.R 樣品.txt
#   Rscript raman_pipeline.R 資料夾/ -o 輸出資料夾
#   Rscript raman_pipeline.R 未知樣.txt --refs 孔雀藍_1.txt 鈦白_1.txt
#   Rscript raman_pipeline.R --list-library
#
# 作者：為「大一顏料科學」課程編寫，2026-08
# =====================================================================

suppressPackageStartupMessages(library(Matrix))

# ---------------------------------------------------------------------
# 第 0 部分：內建顏料譜庫
#   bands = 文獻帶位；key = 關鍵帶；main = 文獻上最強的那一帶
#   文獻依據：Burgio & Clark (2001) Spectrochim. Acta A 57, 1491–1521
# ---------------------------------------------------------------------
LIBRARY <- list(
  "辰砂 cinnabar (HgS)" = list(
    bands = c(253, 284, 343), key = c(253, 343), main = 253,
    note = "硃砂／銀硃。古代至今通用的紅色顏料。"),
  "金紅石 rutile (TiO2)" = list(
    bands = c(144, 232, 447, 609), key = c(447, 609), main = 447,
    note = "鈦白 PW6 的一種晶型。顏料級 1938 年起商業生產。"),
  "銳鈦礦 anatase (TiO2)" = list(
    bands = c(143, 396, 516, 639), key = c(143, 396), main = 143,
    note = "鈦白的另一晶型。顏料級 1918 年起商業生產。143 帶極強。"),
  "酞菁藍 PB15 (CuPc)" = list(
    bands = c(592, 681, 747, 952, 1007, 1106, 1339, 1450, 1527, 1591),
    key = c(747, 1527), main = 1527,
    note = "銅酞菁。1935 年 11 月以 Monastral Blue 之名上市。"),
  "酞菁綠 PG7 (Cl-CuPc)" = list(
    bands = c(685, 742, 776, 818, 977, 1080, 1215, 1281, 1339, 1538),
    key = c(776, 1538), main = 1538,
    note = "氯化銅酞菁。與 PB15 最大差別在 776 與 1538。"),
  "普魯士藍 PB27" = list(
    bands = c(276, 538, 950, 2091, 2154), key = c(2154), main = 2154,
    note = "亞鐵氰化鐵。1704 年發明。2154 落在拉曼靜默區，極好認。"),
  "群青 ultramarine PB29" = list(
    bands = c(258, 548, 808, 1096), key = c(548), main = 548,
    note = "天然為青金石，1828 年起有合成品。"),
  "石青 azurite" = list(
    bands = c(250, 403, 545, 839, 1098, 1580), key = c(403, 839), main = 403,
    note = "鹼式碳酸銅。傳統藍色礦物顏料。"),
  "鉛丹 red lead (Pb3O4)" = list(
    bands = c(121, 149, 313, 390, 548), key = c(313, 548), main = 548,
    note = "四氧化三鉛。橘紅色。"),
  "赤鐵礦 hematite (Fe2O3)" = list(
    bands = c(225, 292, 411, 613), key = c(292, 411), main = 292,
    note = "氧化鐵紅／代赭。"),
  "方解石 calcite (CaCO3)" = list(
    bands = c(282, 712, 1086), key = c(1086, 712), main = 1086,
    note = "常見填料／地仗層材料。"),
  "硫酸鋇 barite (BaSO4)" = list(
    bands = c(453, 617, 988, 1083, 1140), key = c(988), main = 988,
    note = "重晶石／立德粉成分。988 為最強帶，沒有它就不是。"),
  "石膏 gypsum (CaSO4·2H2O)" = list(
    bands = c(414, 493, 619, 1008, 1135), key = c(1008), main = 1008,
    note = "常見地仗層材料。"),
  "聚苯乙烯系樹脂" = list(
    bands = c(621, 795, 1001, 1031, 1450, 1583, 1602, 2905, 3055),
    key = c(1001, 1602), main = 1001,
    note = "現代合成黏合劑。1001 是單取代苯環呼吸振動，又尖又強。"),
  "油／樹脂類黏合劑 (C-H, C=O)" = list(
    bands = c(1440, 1740, 2850, 2930), key = c(2930), main = 2930,
    note = "泛指有機黏合劑，不具專一性，只能當輔助資訊。")
)

# ---------------------------------------------------------------------
# 第 1 部分：讀檔（自動判斷分隔符與表頭）
# ---------------------------------------------------------------------
load_spectrum <- function(path) {
  head5 <- readLines(path, n = 5, warn = FALSE)
  ok <- function(line, sep) {
    p <- if (is.null(sep)) strsplit(trimws(line), "[[:space:]]+")[[1]]
         else strsplit(trimws(line), sep, fixed = TRUE)[[1]]
    if (length(p) < 2) return(FALSE)
    v <- suppressWarnings(as.numeric(p[1:2]))
    all(!is.na(v))
  }
  seps <- list(",", "\t", ";", NULL)
  delim <- NA; skip <- 0; found <- FALSE
  for (s in 0:1) {
    for (d in seps) {
      if (length(head5) > s && ok(head5[s + 1], d)) {
        delim <- d; skip <- s; found <- TRUE; break
      }
    }
    if (found) break
  }
  if (!found) stop(sprintf("看不懂這個檔案的格式：%s", path))
  df <- if (is.null(delim))
    read.table(path, skip = skip, header = FALSE, colClasses = "numeric",
               comment.char = "")
  else
    read.table(path, sep = delim, skip = skip, header = FALSE,
               colClasses = "numeric", comment.char = "")
  x <- df[[1]]; y <- df[[2]]
  o <- order(x)                 # 有些儀器由高波數往低存，統一排序
  list(x = x[o], y = y[o])
}

# ---------------------------------------------------------------------
# 第 2 部分：ALS 基線 + Savitzky-Golay 平滑
# ---------------------------------------------------------------------
als_baseline <- function(y, lambda = 1e5, p = 0.01, n_iter = 12) {
  L <- length(y)
  # 二階差分矩陣 D（L-2 × L），懲罰項 lambda * D'D 控制基線的平滑度
  D <- bandSparse(L - 2, L,
                  k = c(0, 1, 2),
                  diagonals = list(rep(1, L - 2), rep(-2, L - 2), rep(1, L - 2)))
  DtD <- lambda * crossprod(D)
  w <- rep(1, L); z <- y
  for (i in seq_len(n_iter)) {
    W <- Diagonal(x = w)
    z <- as.numeric(solve(W + DtD, w * y))
    # 高於基線的點（也就是峰）權重壓到 p，讓基線沿著谷底走
    w <- p * (y > z) + (1 - p) * (y < z)
  }
  z
}

sg_filter <- function(y, window = 9, poly = 3) {
  if (window %% 2 == 0) window <- window + 1
  if (window < 5 || window >= length(y)) return(y)
  half <- (window - 1) / 2
  t <- -half:half
  A <- outer(t, 0:poly, "^")           # Vandermonde 矩陣
  # 最小平方的帽子矩陣中間那一列，就是平滑用的卷積係數
  coef <- solve(crossprod(A), t(A))[1, ]
  n <- length(y)
  ypad <- c(rep(y[1], half), y, rep(y[n], half))
  out <- stats::filter(ypad, rev(coef), sides = 2)
  as.numeric(out[(half + 1):(half + n)])
}

preprocess <- function(x, y, lambda = 1e5, p = 0.01) {
  base <- als_baseline(y, lambda, p)
  yc <- sg_filter(y - base, 9, 3)
  noise <- sd(diff(yc)) / sqrt(2)      # 用相鄰點差估雜訊，不需要空白區間
  list(yc = yc, base = base, noise = noise)
}

# ---------------------------------------------------------------------
# 第 3 部分：找峰（自行實作 prominence，對應 scipy.signal.find_peaks）
# ---------------------------------------------------------------------
peak_prominence <- function(y, i) {
  # 突出高度：從峰頂往左右走，各自走到「遇到更高的點」為止，
  # 取這段路上的最低點；兩邊最低點取較高者，峰高減掉它就是 prominence。
  n <- length(y); h <- y[i]
  j <- i; lmin <- h
  while (j > 1) { j <- j - 1; if (y[j] > h) break; lmin <- min(lmin, y[j]) }
  j <- i; rmin <- h
  while (j < n) { j <- j + 1; if (y[j] > h) break; rmin <- min(rmin, y[j]) }
  h - max(lmin, rmin)
}

detect_peaks <- function(x, yc, noise, snr = 6, xmin = 150, xmax = 1800,
                         min_sep = 4) {
  # xmin=150：低於此處是雷射濾光片的截止邊緣，會有一個假的「峰」，必須排除
  m <- which(x >= xmin & x <= xmax)
  xs <- x[m]; ys <- yc[m]; n <- length(ys)
  if (n < 5) return(data.frame())
  cand <- which(ys[2:(n - 1)] > ys[1:(n - 2)] & ys[2:(n - 1)] >= ys[3:n]) + 1
  if (!length(cand)) return(data.frame())
  prom <- vapply(cand, function(i) peak_prominence(ys, i), numeric(1))
  keep <- prom >= snr * noise
  cand <- cand[keep]; prom <- prom[keep]
  if (!length(cand)) return(data.frame())
  o <- order(prom, decreasing = TRUE)   # 由強到弱，強的優先佔位
  sel <- integer(0)
  for (i in o) if (all(abs(cand[i] - cand[sel]) >= min_sep)) sel <- c(sel, i)
  data.frame(position = xs[cand[sel]], height = ys[cand[sel]],
             prominence = prom[sel], snr = prom[sel] / noise)
}

# ---------------------------------------------------------------------
# 第 4 部分：比對內建譜庫
# ---------------------------------------------------------------------
band_height <- function(x, yc, center, tol = 6) {
  m <- which(x > center - tol & x < center + tol)
  if (!length(m)) NA_real_ else max(yc[m])
}

match_library <- function(x, yc, noise, peak_pos, tol = 6,
                          strong = 8, weak = 4) {
  # 一個文獻帶要算「出現」，必須 (1) 強度夠、(2) 該處真的偵測到一個峰。
  # 只看強度會被別的顏料強帶的肩部騙到（例：普魯士藍 532 的右肩誤中群青 548）。
  smax <- max(yc[x >= 150 & x <= min(max(x), 3200)])
  near <- function(b) length(peak_pos) > 0 && min(abs(peak_pos - b)) <= tol
  out <- list()
  for (nm in names(LIBRARY)) {
    info <- LIBRARY[[nm]]
    marks <- character(0); hs <- numeric(0); bs <- numeric(0); in_range <- 0
    for (b in info$bands) {
      if (b < min(x) + tol || b > max(x) - tol) {
        bs <- c(bs, b); hs <- c(hs, NA); marks <- c(marks, "界外"); next
      }
      in_range <- in_range + 1
      h <- band_height(x, yc, b, tol)
      mk <- if (h > strong * noise && near(b)) "O"
            else if (h > weak * noise) "~" else "X"
      bs <- c(bs, b); hs <- c(hs, h); marks <- c(marks, mk)
    }
    hits <- sum(marks == "O"); partial <- sum(marks == "~")
    score <- if (in_range) (hits + 0.5 * partial) / in_range else 0
    mb <- info$main
    main_ratio <- if (is.null(mb) || mb < min(x) + tol || mb > max(x) - tol ||
                      smax <= 0) NA_real_
                  else band_height(x, yc, mb, tol) / smax
    ks <- vapply(info$key, function(k) {
      if (k < min(x) + tol || k > max(x) - tol) NA
      else (band_height(x, yc, k, tol) > strong * noise) && near(k)
    }, logical(1))
    key_ok <- any(!is.na(ks)) && all(ks[!is.na(ks)])
    strong_main <- !is.na(main_ratio) && main_ratio >= 0.25
    verdict <- if (key_ok && score >= 0.5 && strong_main) "主成分"
               else if (key_ok && score >= 0.5) "次要成分"
               else if (key_ok) "存疑" else "不成立"
    out[[nm]] <- list(name = nm, note = info$note, bands = bs, heights = hs,
                      marks = marks, hits = hits, partial = partial,
                      in_range = in_range, score = score, main = mb,
                      main_ratio = main_ratio, verdict = verdict)
  }
  rank <- c("主成分" = 3, "次要成分" = 2, "存疑" = 1, "不成立" = 0)
  out[order(-rank[vapply(out, function(r) r$verdict, "")],
            -vapply(out, function(r) r$score, 0))]
}

# ---------------------------------------------------------------------
# 第 5 部分：NNLS 解混（Lawson–Hanson 主動集法，自行實作）
# ---------------------------------------------------------------------
nnls_fit <- function(A, b, max_iter = 300, tol = 1e-10) {
  n <- ncol(A); P <- logical(n); xx <- rep(0, n)
  w <- as.numeric(crossprod(A, b - A %*% xx))
  it <- 0
  while (any(!P) && max(w[!P]) > tol && it < max_iter) {
    it <- it + 1
    j <- which(!P)[which.max(w[!P])]
    P[j] <- TRUE
    s <- rep(0, n)
    s[P] <- tryCatch(qr.solve(A[, P, drop = FALSE], b),
                     error = function(e) xx[P])
    inner <- 0
    while (min(s[P]) <= 0 && inner < max_iter) {
      inner <- inner + 1
      neg <- P & (s <= 0)
      alpha <- min(xx[neg] / (xx[neg] - s[neg]))
      xx <- xx + alpha * (s - xx)
      P[P & abs(xx) < tol] <- FALSE
      s <- rep(0, n)
      if (!any(P)) break
      s[P] <- tryCatch(qr.solve(A[, P, drop = FALSE], b),
                       error = function(e) rep(0, sum(P)))
    }
    xx <- s
    w <- as.numeric(crossprod(A, b - A %*% xx))
  }
  pmax(xx, 0)
}

unmix <- function(x, yc, refs, lo = 350, hi = 1750) {
  g <- which(x >= lo & x <= hi)
  X <- x[g]; Y <- yc[g]
  cols <- list(); nms <- character(0)
  for (r in refs) {
    v <- approx(r$x, r$yc, X, rule = 2)$y
    mx <- max(v); cols[[length(cols) + 1]] <- if (mx > 0) v / mx else v
    nms <- c(nms, r$name)
  }
  cols[[length(cols) + 1]] <- rep(1, length(X)); nms <- c(nms, "(常數項)")
  A <- do.call(cbind, cols)
  coef <- nnls_fit(A, Y)
  fit <- as.numeric(A %*% coef); res <- Y - fit
  r <- suppressWarnings(cor(fit, Y))
  list(x = X, y = Y, fit = fit, residual = res, names = nms, coef = coef,
       r = r, r2 = r^2, res_ratio = sd(res) / sd(Y))
}

# ---------------------------------------------------------------------
# 第 6 部分：出圖
# ---------------------------------------------------------------------
plot_report <- function(x, y, yc, base, peaks, top, title, out_png, u = NULL) {
  n <- if (is.null(u)) 2 else 3
  png(out_png, width = 1400, height = 420 * n, res = 130)
  op <- par(mfrow = c(n, 1), mar = c(4.2, 4.5, 3, 1), family = "sans")
  on.exit({ par(op); dev.off() }, add = TRUE)

  plot(x, y, type = "l", col = "grey60", lwd = 1,
       main = paste0(title, "  步驟 1：估計並扣除螢光背景"),
       xlab = "拉曼位移 (cm-1)", ylab = "強度")
  lines(x, base, col = "red", lty = 2, lwd = 1.6)
  legend("topright", c("原始光譜", "ALS 估計的基線"),
         col = c("grey60", "red"), lty = c(1, 2), bty = "n", cex = .85)

  m <- x <= 1800
  plot(x[m], yc[m], type = "l", col = "grey20", lwd = 1.3,
       main = paste0("步驟 2：找峰（標出前 12 強）",
                     if (!is.null(top)) paste0("　最相符：", top) else ""),
       xlab = "拉曼位移 (cm-1)", ylab = "強度")
  if (nrow(peaks)) {
    k <- head(order(-peaks$prominence), 12)
    text(peaks$position[k], peaks$height[k], sprintf("%.0f", peaks$position[k]),
         pos = 3, cex = .7, col = "firebrick", font = 2)
  }

  if (!is.null(u)) {
    off <- max(u$y) * 0.55
    yl <- range(c(u$y, u$residual - off))
    plot(u$x, u$y, type = "l", col = "grey20", lwd = 1.3, ylim = yl,
         main = "步驟 3：NNLS 解混——殘差裡若還有成組的峰，代表有沒被解釋的成分",
         xlab = "拉曼位移 (cm-1)", ylab = "強度")
    lines(u$x, u$fit, col = "darkorange", lwd = 1.2)
    lines(u$x, u$residual - off, col = "purple", lwd = 1)
    abline(h = -off, col = "grey80")
    legend("topleft", c("實測", sprintf("NNLS 擬合 (r=%.3f)", u$r), "殘差（往下位移）"),
           col = c("grey20", "darkorange", "purple"), lty = 1, bty = "n", cex = .85)
  }
}

# ---------------------------------------------------------------------
# 第 7 部分：分析一個檔案
# ---------------------------------------------------------------------
analyze <- function(path, outdir, refs = NULL, snr = 6, lambda = 1e5,
                    make_plot = TRUE) {
  nm <- tools::file_path_sans_ext(basename(path))
  sp <- load_spectrum(path); x <- sp$x; y <- sp$y
  pr <- preprocess(x, y, lambda); yc <- pr$yc; base <- pr$base; noise <- pr$noise
  peaks <- detect_peaks(x, yc, noise, snr = snr)

  allp <- detect_peaks(x, yc, noise, snr = 4, xmin = 150, xmax = max(x))
  matches <- match_library(x, yc, noise,
                           peak_pos = if (nrow(allp)) allp$position else numeric(0))

  usable <- nrow(peaks) > 0 && max(peaks$snr) >= 20
  if (!usable)
    matches <- lapply(matches, function(r) {
      if (r$verdict %in% c("主成分", "次要成分")) r$verdict <- "存疑"; r })

  L <- character(0); A <- function(s) L <<- c(L, s)
  A(strrep("=", 68))
  A(sprintf("樣品：%s", nm))
  A(sprintf("檔案：%s", path))
  A(sprintf("波數範圍：%.1f - %.1f cm-1　資料點：%d", min(x), max(x), length(x)))
  A(sprintf("雜訊 sigma ≈ %.1f　訊背比＝%.2f", noise, max(y) / median(base)))
  if (!usable) {
    A("")
    A("!! 資料品質警告：全譜找不到任何 S/N >= 20 的峰。")
    A("   這張光譜很可能是量測失敗（雷射離焦、樣品燒毀、積分時間不足，")
    A("   或該點只有螢光背景）。以下判定一律降級為『存疑』，不可採信。")
  }
  A("")
  A(sprintf("-- 偵測到的峰（S/N >= %g，150-1800 cm-1）--", snr))
  if (!nrow(peaks)) A("   （沒有超過門檻的峰）")
  else {
    k <- head(order(-peaks$prominence), 20)
    for (i in k)
      A(sprintf("   %8.1f cm-1   高度=%9.0f   S/N=%6.1f",
                peaks$position[i], peaks$height[i], peaks$snr[i]))
  }
  A("")
  A("-- 內建譜庫比對　三道門檻：關鍵帶齊全／命中率過半／主帶強度 >= 全譜最強峰 25% --")
  icon <- c("主成分" = "* 主成分  ", "次要成分" = "+ 次要成分",
            "存疑" = "~ 存疑    ", "不成立" = ". 不成立  ")
  for (r in head(matches, 6)) {
    mr <- if (is.na(r$main_ratio)) "--" else sprintf("%.0f%%", r$main_ratio * 100)
    A(sprintf("   %s　%s　命中 %dO+%d~/%d（%.0f%%）　主帶 %s 強度佔比 %s",
              icon[[r$verdict]], r$name, r$hits, r$partial, r$in_range,
              r$score * 100, as.character(r$main), mr))
    A(paste0("       ", paste(sprintf("%g:%s%s", r$bands,
        ifelse(is.na(r$heights), "--", sprintf("%.0f", r$heights)), r$marks),
        collapse = "  ")))
  }
  mainc <- vapply(matches, function(r) r$verdict, "") == "主成分"
  minorc <- vapply(matches, function(r) r$verdict, "") == "次要成分"
  nm_of <- function(sel) if (any(sel))
    paste(vapply(matches[sel], function(r) r$name, ""), collapse = "　＋　") else "（無）"
  A("")
  A("-- 初步判定 --")
  A(sprintf("   主成分　：%s", nm_of(mainc)))
  A(sprintf("   次要成分：%s", nm_of(minorc)))
  A("   ※ 這是程式的機械比對結果，不等於結論。三個常見陷阱：")
  A("     (1) 強帶的肩部會誤觸別的顏料　(2) 量測失敗的譜程式仍會硬湊名字")
  A("     (3) 譜庫只有 15 種，庫外的成分會被最接近的那個吃掉")

  top <- if (any(mainc)) matches[mainc][[1]]$name else
         if (any(minorc)) matches[minorc][[1]]$name else NULL

  u <- NULL
  if (length(refs)) {
    loaded <- lapply(refs, function(rp) {
      s <- load_spectrum(rp); q <- preprocess(s$x, s$y, lambda)
      list(name = tools::file_path_sans_ext(basename(rp)), x = s$x, yc = q$yc)
    })
    u <- unmix(x, yc, loaded)
    A(""); A("-- NNLS 混合物解混 --")
    for (i in seq_along(u$names))
      A(sprintf("   %-28s 係數 = %10.1f", u$names[i], u$coef[i]))
    A(sprintf("   擬合 r = %.3f　(R2 = %.3f)　殘差RMS/訊號RMS = %.3f",
              u$r, u$r2, u$res_ratio))
    rn <- sd(diff(u$residual)) / sqrt(2)
    rp <- detect_peaks(u$x, u$residual, rn, snr = 6, xmin = min(u$x),
                       xmax = max(u$x))
    if (nrow(rp)) {
      A("   殘差中仍存在的峰（代表有參考譜沒涵蓋到的成分）：")
      k <- head(order(-rp$prominence), 12)
      for (i in k[order(rp$position[k])])
        A(sprintf("      %8.1f cm-1   S/N=%5.1f", rp$position[i], rp$snr[i]))
    } else A("   殘差中沒有明顯的峰 → 參考譜已能解釋這張光譜。")
  }

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  writeLines(L, file.path(outdir, paste0(nm, "_report.txt")), useBytes = TRUE)
  if (nrow(peaks)) {
    df <- peaks[order(-peaks$prominence), ]
    names(df) <- c("峰位_cm-1", "高度", "突出高度", "訊噪比SN")
    con <- file(file.path(outdir, paste0(nm, "_peaks.csv")), "w",
                encoding = "UTF-8")
    writeLines("﻿", con, sep = "")
    write.csv(df, con, row.names = FALSE)
    close(con)
  }
  if (make_plot)
    try(plot_report(x, y, yc, base, peaks, top, nm,
                    file.path(outdir, paste0(nm, "_report.png")), u), silent = TRUE)
  paste(L, collapse = "\n")
}

# ---------------------------------------------------------------------
# 第 8 部分：命令列介面
# ---------------------------------------------------------------------
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) || "--help" %in% args || "-h" %in% args) {
    cat("用法：Rscript raman_pipeline.R <檔案或資料夾> [-o 輸出夾]",
        "[--refs a.txt b.txt] [--snr 6] [--no-plot] [--list-library]\n")
    return(invisible(0))
  }
  if ("--list-library" %in% args) {
    cat(sprintf("內建譜庫共 %d 種：\n\n", length(LIBRARY)))
    for (nm in names(LIBRARY)) {
      i <- LIBRARY[[nm]]
      cat(sprintf("  %s\n      帶位：%s\n      關鍵帶：%s\n      說明：%s\n\n",
                  nm, paste(i$bands, collapse = ", "),
                  paste(i$key, collapse = ", "), i$note))
    }
    return(invisible(0))
  }
  getopt <- function(flag, default = NULL, n = 1) {
    i <- which(args == flag)
    if (!length(i)) return(default)
    if (n == Inf) {
      j <- i + 1; out <- character(0)
      while (j <= length(args) && !startsWith(args[j], "-")) {
        out <- c(out, args[j]); j <- j + 1
      }
      return(out)
    }
    args[i + 1]
  }
  outdir <- getopt("-o", getopt("--outdir", "raman_out"))
  refs <- getopt("--refs", NULL, n = Inf)
  snr <- as.numeric(getopt("--snr", "6"))
  make_plot <- !("--no-plot" %in% args)
  pos <- args[!startsWith(args, "-")]
  used <- c(outdir, refs, as.character(snr))
  input <- setdiff(pos, used)[1]
  if (is.na(input)) { cat("請指定輸入檔案或資料夾\n"); return(invisible(1)) }

  files <- if (dir.exists(input))
    sort(list.files(input, pattern = "\\.(txt|csv)$", full.names = TRUE))
  else input
  for (f in files) {
    r <- try(analyze(f, outdir, refs = refs, snr = snr, make_plot = make_plot),
             silent = TRUE)
    if (inherits(r, "try-error")) cat(sprintf("[跳過] %s\n", f))
    else cat(r, "\n\n")
  }
  cat(sprintf("→ 輸出已寫入：%s\n", normalizePath(outdir)))
  invisible(0)
}

if (sys.nframe() == 0L) main()

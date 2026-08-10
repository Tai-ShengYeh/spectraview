#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
raman_pipeline.py — 顏料拉曼光譜分析流程（教學版）
================================================================
把「讀檔 → 扣基線 → 平滑 → 找峰 → 比對譜庫 → 混合物解混 → 出圖出表」
包成一支可重複執行的腳本。

用法
----
    # 分析單一檔案
    python raman_pipeline.py 樣品.txt

    # 分析整個資料夾
    python raman_pipeline.py 資料夾/ -o 輸出資料夾

    # 加上自己的參考譜做 NNLS 解混（可給多個）
    python raman_pipeline.py 未知樣.txt --refs 孔雀藍_1.txt 鈦白_1.txt

    # 只列出內建譜庫
    python raman_pipeline.py --list-library

作者：為「大一顏料科學」課程編寫，2026-08
"""

import argparse
import csv
import os
import sys
import glob

import numpy as np
from scipy.signal import savgol_filter, find_peaks
from scipy.optimize import nnls, curve_fit
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

# ============================================================
# 第 0 部分：內建顏料譜庫
# ------------------------------------------------------------
# 每個顏料記三件事：
#   bands  : 文獻報告的拉曼帶位置 (cm-1)
#   key    : 「關鍵帶」——這些沒出現就不能判定為此顏料
#   note   : 給人看的說明（年代、化學式）
# 文獻依據：Burgio & Clark (2001) Spectrochim. Acta A 57, 1491–1521
# ============================================================
LIBRARY = {
    "辰砂 cinnabar (HgS)": {
        "bands": [253, 284, 343], "key": [253, 343], "main": 253,
        "note": "硃砂／銀硃。古代至今通用的紅色顏料。"},
    "金紅石 rutile (TiO2)": {
        "bands": [144, 232, 447, 609], "key": [447, 609], "main": 447,
        "note": "鈦白 PW6 的一種晶型。顏料級 1938 年起商業生產。"},
    "銳鈦礦 anatase (TiO2)": {
        "bands": [143, 396, 516, 639], "key": [143, 396], "main": 143,
        "note": "鈦白的另一晶型。顏料級 1918 年起商業生產。143 帶極強。"},
    "酞菁藍 PB15 (CuPc)": {
        "bands": [592, 681, 747, 952, 1007, 1106, 1339, 1450, 1527, 1591],
        "key": [747, 1527], "main": 1527,
        "note": "銅酞菁。1935 年 11 月以 Monastral Blue 之名上市。"},
    "酞菁綠 PG7 (Cl-CuPc)": {
        "bands": [685, 742, 776, 818, 977, 1080, 1215, 1281, 1339, 1538],
        "key": [776, 1538], "main": 1538,
        "note": "氯化銅酞菁。與 PB15 最大差別在 776 與 1538。"},
    "普魯士藍 PB27": {
        "bands": [276, 538, 950, 2091, 2154], "key": [2154], "main": 2154,
        "note": "亞鐵氰化鐵。1704 年發明。2154 落在拉曼靜默區，極好認。"},
    "群青 ultramarine PB29": {
        "bands": [258, 548, 808, 1096], "key": [548], "main": 548,
        "note": "天然為青金石，1828 年起有合成品。"},
    "石青 azurite": {
        "bands": [250, 403, 545, 839, 1098, 1580], "key": [403, 839], "main": 403,
        "note": "鹼式碳酸銅。傳統藍色礦物顏料。"},
    "鉛丹 red lead (Pb3O4)": {
        "bands": [121, 149, 313, 390, 548], "key": [313, 548], "main": 548,
        "note": "四氧化三鉛。橘紅色。"},
    "赤鐵礦 hematite (Fe2O3)": {
        "bands": [225, 292, 411, 613], "key": [292, 411], "main": 292,
        "note": "氧化鐵紅／代赭。"},
    "方解石 calcite (CaCO3)": {
        "bands": [282, 712, 1086], "key": [1086, 712], "main": 1086,
        "note": "常見填料／地仗層材料。"},
    "硫酸鋇 barite (BaSO4)": {
        "bands": [453, 617, 988, 1083, 1140], "key": [988], "main": 988,
        "note": "重晶石／立德粉成分。988 為最強帶，沒有它就不是。"},
    "石膏 gypsum (CaSO4·2H2O)": {
        "bands": [414, 493, 619, 1008, 1135], "key": [1008], "main": 1008,
        "note": "常見地仗層材料。"},
    "聚苯乙烯系樹脂": {
        "bands": [621, 795, 1001, 1031, 1450, 1583, 1602, 2905, 3055],
        "key": [1001, 1602], "main": 1001,
        "note": "現代合成黏合劑。1001 是單取代苯環呼吸振動，又尖又強。"},
    "油／樹脂類黏合劑 (C–H, C=O)": {
        "bands": [1440, 1740, 2850, 2930], "key": [2930], "main": 2930,
        "note": "泛指有機黏合劑，不具專一性，只能當輔助資訊。"},
}


# ============================================================
# 第 1 部分：讀檔
# ============================================================
def load_spectrum(path):
    """讀入兩欄的拉曼光譜檔，自動判斷分隔符號與是否有表頭。

    回傳 (x, y)：x 是拉曼位移 cm-1，y 是強度。
    """
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        head = [f.readline() for _ in range(5)]

    # 自動偵測分隔符與表頭：逗號 / Tab / 分號 / 空白，四種都試一遍。
    # 注意 None 在 numpy 代表「以任意空白切」，所以不能拿 None 當「找不到」的標記，
    # 要另外用 found 旗標。
    def try_split(line, d):
        parts = line.strip().split(d) if d is not None else line.split()
        if len(parts) < 2:
            return False
        try:
            float(parts[0]); float(parts[1])
            return True
        except ValueError:
            return False

    delim, skip, found = None, 0, False
    for s in (0, 1):                       # s=0 沒表頭，s=1 有一行表頭
        for d in (",", "\t", ";", None):
            if s < len(head) and head[s] and try_split(head[s], d):
                delim, skip, found = d, s, True
                break
        if found:
            break
    if not found:
        raise ValueError(f"看不懂這個檔案的格式：{path}")
    data = np.loadtxt(path, delimiter=delim, skiprows=skip, usecols=(0, 1))

    x, y = data[:, 0], data[:, 1]
    order = np.argsort(x)          # 有些儀器由高波數往低波數存，統一排序
    return x[order], y[order]


# ============================================================
# 第 2 部分：基線扣除（ALS，非對稱最小平方）
# ============================================================
def als_baseline(y, lam=1e5, p=0.01, n_iter=12):
    """用非對稱最小平方法估計螢光背景。

    lam 越大基線越硬（越直），p 越小越貼著谷底走。
    對拉曼光譜，lam=1e5、p=0.01 是常用起點。
    """
    L = len(y)
    D = diags([1, -2, 1], [0, -1, -2], shape=(L, L - 2), dtype=float)
    D = lam * D.dot(D.T)
    w = np.ones(L)
    W = diags(w, 0)
    z = y.copy()
    for _ in range(n_iter):
        W.setdiag(w)
        z = spsolve((W + D).tocsc(), w * y)
        w = p * (y > z) + (1 - p) * (y < z)   # 峰的部分權重壓低，讓基線走谷底
    return z


def preprocess(x, y, lam=1e5, p=0.01, window=9, poly=3):
    """扣基線 + Savitzky-Golay 平滑。回傳 (校正後強度, 基線, 雜訊σ)。"""
    base = als_baseline(y, lam=lam, p=p)
    yc = y - base
    if window >= 5 and window < len(yc):
        if window % 2 == 0:
            window += 1
        yc = savgol_filter(yc, window, poly)
    # 雜訊估計：相鄰點差的標準差 / √2（不依賴任何「空白區間」的假設）
    noise = float(np.std(np.diff(yc)) / np.sqrt(2))
    return yc, base, noise


# ============================================================
# 第 3 部分：找峰
# ============================================================
def detect_peaks(x, yc, noise, snr=6.0, xmin=150.0, xmax=1800.0, min_sep=4):
    """找出訊噪比高於 snr 的峰。回傳 list of dict。

    xmin 預設 150 cm⁻¹：低於這裡是雷射濾光片的截止邊緣，會有一個假的「峰」，
    不是樣品的訊號，一定要排除，否則後面所有相對強度的計算都會被它帶歪。
    """
    m = (x >= xmin) & (x <= xmax)
    xm, ym = x[m], yc[m]
    idx, props = find_peaks(ym, prominence=snr * noise, distance=min_sep)
    peaks = []
    for i, j in enumerate(idx):
        peaks.append({
            "position": float(xm[j]),
            "height": float(ym[j]),
            "prominence": float(props["prominences"][i]),
            "snr": float(props["prominences"][i] / noise),
        })
    peaks.sort(key=lambda d: -d["prominence"])
    return peaks


def _lorentzian(x, a, x0, g, c):
    return a * g ** 2 / ((x - x0) ** 2 + g ** 2) + c


def refine_peak(x, yc, center, halfwidth=20.0):
    """在 center 附近用 Lorentzian 擬合，取得更精確的峰心與半高寬。"""
    m = (x > center - halfwidth) & (x < center + halfwidth)
    if m.sum() < 6:
        return None
    xi, yi = x[m], yc[m]
    try:
        p0 = [yi.max() - np.median(yi), xi[np.argmax(yi)], 8.0, float(np.median(yi))]
        popt, _ = curve_fit(_lorentzian, xi, yi, p0=p0, maxfev=40000)
        fwhm = float(abs(2 * popt[2]))
        # 擬合失控的情況要擋掉：峰心跑掉、或半高寬寬到超過擬合視窗
        if abs(popt[1] - center) > halfwidth or fwhm > 4 * halfwidth or fwhm <= 0:
            return None
        return {"center": float(popt[1]), "height": float(popt[0]), "fwhm": fwhm}
    except Exception:
        return None


# ============================================================
# 第 4 部分：比對內建譜庫
# ============================================================
def band_height(x, yc, center, tol=8.0):
    """取 center ± tol 範圍內的最大值，當作該帶的高度。"""
    m = (x > center - tol) & (x < center + tol)
    return float(yc[m].max()) if m.any() else float("nan")


def match_library(x, yc, noise, peak_positions=None, tol=6.0, strong=8.0, weak=4.0):
    """把光譜和內建譜庫逐一比對。

    一個文獻帶要算「出現(✔)」必須同時滿足兩件事：
      (1) 該位置的強度 > strong 倍雜訊
      (2) 附近 tol cm⁻¹ 內真的偵測到一個「峰」

    第 (2) 點很關鍵。只看強度的話，別的顏料強帶的「肩部」會讓不存在的顏料
    誤判為出現——例如普魯士藍 532 的右肩強度很高，會誤觸群青的 548。
    要求該位置本身是個峰，就能擋掉這類假陽性。
    """
    if peak_positions is None:
        idx, _ = find_peaks(yc, prominence=weak * noise, distance=3)
        peak_positions = x[idx]
    peak_positions = np.asarray(peak_positions, dtype=float)

    def near_peak(b):
        if peak_positions.size == 0:
            return False
        return bool(np.min(np.abs(peak_positions - b)) <= tol)

    results = []
    for name, info in LIBRARY.items():
        detail, in_range = [], 0
        for b in info["bands"]:
            if b < x.min() + tol or b > x.max() - tol:
                detail.append((b, None, "界外"))
                continue
            in_range += 1
            h = band_height(x, yc, b, tol)
            if h > strong * noise and near_peak(b):
                mark = "✔"
            elif h > weak * noise:
                mark = "△"
            else:
                mark = "✘"
            detail.append((b, h, mark))
        hits = sum(1 for _, _, mk in detail if mk == "✔")
        partial = sum(1 for _, _, mk in detail if mk == "△")
        # 主帶強度檢查：該顏料文獻上最強的那一帶，在本張光譜裡有多強？
        # 若一個顏料是樣品的主成分，它的主帶通常會是整張譜數一數二的峰。
        mb = info.get("main")
        span = yc[(x >= 150) & (x <= min(x.max(), 3200))]   # 同樣要避開濾光片邊緣
        smax = float(span.max()) if span.size else 1.0
        if mb is None or mb < x.min() + tol or mb > x.max() - tol or smax <= 0:
            main_ratio = float("nan")
        else:
            main_ratio = band_height(x, yc, mb, tol) / smax
        # 關鍵帶檢查
        key_state = []
        for k in info["key"]:
            if k < x.min() + tol or k > x.max() - tol:
                key_state.append(None)                     # 看不到，無法判斷
            else:
                key_state.append(
                    band_height(x, yc, k, tol) > strong * noise and near_peak(k))
        key_ok = all(s for s in key_state if s is not None) and any(
            s is not None for s in key_state)
        # 命中率：完全出現算 1 分，疑似算 0.5 分
        score = (hits + 0.5 * partial) / in_range if in_range else 0.0
        # 判定分四級：
        #   主成分   關鍵帶齊全 + 命中率過半 + 主帶夠強（≥ 全譜最強峰的 25%）
        #   次要成分 前兩項過，但主帶不夠強 → 可能真的存在，只是含量低
        #   存疑     只有關鍵帶過，命中率不足 → 多半是別的強帶肩部造成的假陽性
        #   不成立   關鍵帶沒到齊
        strong_main = (not np.isnan(main_ratio)) and main_ratio >= 0.25
        if key_ok and score >= 0.5 and strong_main:
            verdict = "主成分"
        elif key_ok and score >= 0.5:
            verdict = "次要成分"
        elif key_ok:
            verdict = "存疑"
        else:
            verdict = "不成立"
        results.append({"name": name, "note": info["note"], "detail": detail,
                        "hits": hits, "partial": partial, "in_range": in_range,
                        "score": score, "key_ok": key_ok,
                        "main": mb, "main_ratio": main_ratio, "verdict": verdict})
    rank = {"主成分": 3, "次要成分": 2, "存疑": 1, "不成立": 0}
    results.sort(key=lambda r: (rank[r["verdict"]], r["score"]), reverse=True)
    return results


# ============================================================
# 第 5 部分：混合物解混（NNLS）
# ============================================================
def unmix(x, yc, refs, lo=350.0, hi=1750.0):
    """用非負最小平方把未知樣拆成幾個參考譜的線性組合。

    refs: [(名稱, x_ref, y_ref), ...]
    回傳 dict，含係數、擬合曲線、殘差。
    注意：不同儀器之間解析度與帶形不同，係數只能當半定量參考，
          真正可信的是「殘差裡還剩下什麼峰」。
    """
    g = (x >= lo) & (x <= hi)
    X, Y = x[g], yc[g]
    cols, names = [], []
    for nm, xr, yr in refs:
        r = np.interp(X, xr, yr)
        mx = r.max()
        cols.append(r / mx if mx > 0 else r)
        names.append(nm)
    cols.append(np.ones_like(X))          # 常數項，吸收殘餘背景
    names.append("(常數項)")
    A = np.vstack(cols).T
    coef, _ = nnls(A, Y)
    fit = A @ coef
    res = Y - fit
    r = float(np.corrcoef(fit, Y)[0, 1]) if Y.std() > 0 else float("nan")
    return {"x": X, "y": Y, "fit": fit, "residual": res, "names": names,
            "coef": coef, "r": r, "r2": r ** 2,
            "res_ratio": float(res.std() / Y.std())}


# ============================================================
# 第 6 部分：出圖
# ============================================================
def _setup_font():
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import font_manager as fm
    import matplotlib.pyplot as plt
    for p in ["/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
              "C:/Windows/Fonts/msjh.ttc", "C:/Windows/Fonts/mingliu.ttc"]:
        if os.path.exists(p):
            try:
                fm.fontManager.addfont(p)
                plt.rcParams["font.family"] = [
                    fm.FontProperties(fname=p).get_name(), "DejaVu Sans"]
                break
            except Exception:
                pass
    plt.rcParams["axes.unicode_minus"] = False
    return plt


def plot_report(x, y, yc, base, peaks, top_match, title, out_png, unmix_res=None):
    plt = _setup_font()
    n = 3 if unmix_res is not None else 2
    fig, axes = plt.subplots(n, 1, figsize=(12, 4.0 * n))
    if n == 2:
        axes = list(axes)

    ax = axes[0]
    ax.plot(x, y, color="#95a5a6", lw=.9, label="原始光譜")
    ax.plot(x, base, color="#e74c3c", lw=1.1, ls="--", label="ALS 估計的基線")
    ax.set_title(f"{title}　步驟 1：估計並扣除螢光背景", fontsize=12,
                 fontweight="bold", loc="left")
    ax.legend(fontsize=9); ax.set_xlabel("拉曼位移 (cm⁻¹)"); ax.set_ylabel("強度")

    ax = axes[1]
    m = x <= 1800
    ax.plot(x[m], yc[m], color="#2c3e50", lw=1.2, label="扣基線＋平滑後")
    for pk in peaks[:12]:
        ax.annotate(f"{pk['position']:.0f}", (pk["position"], pk["height"]),
                    textcoords="offset points", xytext=(0, 8), ha="center",
                    fontsize=8, color="#c0392b", fontweight="bold")
    sub = f"　最相符：{top_match}" if top_match else ""
    ax.set_title(f"步驟 2：找峰（標出前 12 強）{sub}", fontsize=12,
                 fontweight="bold", loc="left")
    ax.legend(fontsize=9); ax.set_xlabel("拉曼位移 (cm⁻¹)"); ax.set_ylabel("強度")

    if unmix_res is not None:
        ax = axes[2]
        u = unmix_res
        ax.plot(u["x"], u["y"], color="#2c3e50", lw=1.3, label="實測")
        ax.plot(u["x"], u["fit"], color="#e67e22", lw=1.1, label=f"NNLS 擬合 (r={u['r']:.3f})")
        off = u["y"].max() * 0.55
        ax.plot(u["x"], u["residual"] - off, color="#8e44ad", lw=1.0, label="殘差（往下位移）")
        ax.axhline(-off, color="#ccc", lw=.6)
        ax.set_title("步驟 3：NNLS 解混——殘差裡若還有成組的峰，代表有沒被解釋的成分",
                     fontsize=12, fontweight="bold", loc="left")
        ax.legend(fontsize=9); ax.set_xlabel("拉曼位移 (cm⁻¹)"); ax.set_ylabel("強度")

    plt.tight_layout()
    plt.savefig(out_png, dpi=140, bbox_inches="tight")
    plt.close(fig)


# ============================================================
# 第 7 部分：分析一個檔案（把上面全部串起來）
# ============================================================
def analyze(path, outdir, refs=None, snr=6.0, lam=1e5, make_plot=True):
    name = os.path.splitext(os.path.basename(path))[0]
    x, y = load_spectrum(path)
    yc, base, noise = preprocess(x, y, lam=lam)
    peaks = detect_peaks(x, yc, noise, snr=snr)
    # 譜庫比對要看到 1800 cm⁻¹ 以上（普魯士藍 2154、C–H 2930 都在那邊），
    # 所以另外用較寬鬆的門檻做一次全範圍找峰，只取位置。
    all_idx, _ = find_peaks(yc, prominence=4 * noise, distance=3)
    matches = match_library(x, yc, noise, peak_positions=x[all_idx])

    # 資料品質把關：一張可用的拉曼光譜，至少要有一個 S/N ≥ 20 的峰。
    # 完全沒有的話，程式仍然會硬湊出一個最像的顏料名字，那是假的。
    strong_peaks = [p for p in peaks if p["snr"] >= 20]
    usable = len(strong_peaks) > 0
    if not usable:
        for r in matches:
            if r["verdict"] in ("主成分", "次要成分"):
                r["verdict"] = "存疑"

    lines = []
    A = lines.append
    A("=" * 68)
    A(f"樣品：{name}")
    A(f"檔案：{path}")
    A(f"波數範圍：{x.min():.1f} – {x.max():.1f} cm⁻¹　資料點：{len(x)}")
    A(f"雜訊 σ ≈ {noise:.1f}　訊背比（最大值／背景中位數）＝ "
      f"{y.max() / max(np.median(base), 1e-9):.2f}")
    if not usable:
        A("")
        A("!! 資料品質警告：全譜找不到任何 S/N ≥ 20 的峰。")
        A("   這張光譜很可能是量測失敗（雷射離焦、樣品燒毀、積分時間不足，")
        A("   或該點只有螢光背景）。以下所有判定一律降級為『存疑』，不可採信。")
        A("   建議：降低雷射功率、換一個量測點、拉長積分時間後重測。")
    A("")
    A(f"── 偵測到的峰（S/N ≥ {snr:g}，150–1800 cm⁻¹）──")
    if not peaks:
        A("   （沒有超過門檻的峰——可能是量測失敗，或這個點只有螢光背景）")
    for pk in peaks[:20]:
        fit = refine_peak(x, yc, pk["position"])
        fw = f"　FWHM={fit['fwhm']:.1f}" if fit else ""
        A(f"   {pk['position']:8.1f} cm⁻¹   高度={pk['height']:9.0f}   "
          f"S/N={pk['snr']:6.1f}{fw}")
    A("")
    A("── 內建譜庫比對　三道門檻：關鍵帶齊全／命中率過半／主帶強度 ≥ 全譜最強峰 25% ──")
    icon = {"主成分": "★ 主成分  ", "次要成分": "◆ 次要成分",
            "存疑": "△ 存疑    ", "不成立": "· 不成立  "}
    for r in matches[:6]:
        mr = "--" if np.isnan(r["main_ratio"]) else f"{r['main_ratio'] * 100:.0f}%"
        A(f"   {icon[r['verdict']]}　{r['name']}　命中 {r['hits']}✔+{r['partial']}△"
          f"/{r['in_range']}（{r['score'] * 100:.0f}%）　主帶 {r['main']} 強度佔比 {mr}")
        A("       " + "  ".join(
            f"{b}:{'--' if h is None else f'{h:.0f}'}{mk}" for b, h, mk in r["detail"]))
    main_c = [r["name"] for r in matches if r["verdict"] == "主成分"]
    minor_c = [r["name"] for r in matches if r["verdict"] == "次要成分"]
    top = main_c[0] if main_c else (minor_c[0] if minor_c else None)
    A("")
    A(f"── 初步判定 ──")
    A(f"   主成分　：{'　＋　'.join(main_c) if main_c else '（無）'}")
    A(f"   次要成分：{'　＋　'.join(minor_c) if minor_c else '（無）'}")
    A("   ※ 這是程式的機械比對結果，不等於結論。三個常見陷阱：")
    A("     (1) 強帶的肩部會誤觸別的顏料（普魯士藍 538 的右肩會誤中群青 548）")
    A("     (2) 量測失敗的光譜只有雜訊，程式仍會硬湊出一個名字")
    A("     (3) 譜庫只有 15 種，樣品裡的成分不在庫裡時，會被最接近的那個吃掉")
    A("     務必回頭看逐點數值確認峰心位置，並確認該顏料的次要帶也同時出現。")

    u = None
    if refs:
        loaded = []
        for rp in refs:
            rx, ry = load_spectrum(rp)
            ryc, _, _ = preprocess(rx, ry, lam=lam)
            loaded.append((os.path.splitext(os.path.basename(rp))[0], rx, ryc))
        u = unmix(x, yc, loaded)
        A("")
        A("── NNLS 混合物解混 ──")
        for nm, c in zip(u["names"], u["coef"]):
            A(f"   {nm:28s} 係數 = {c:10.1f}")
        A(f"   擬合 r = {u['r']:.3f}　(R² = {u['r2']:.3f})　"
          f"殘差RMS/訊號RMS = {u['res_ratio']:.3f}")
        rn = float(np.std(np.diff(u["residual"])) / np.sqrt(2))
        ridx, rprops = find_peaks(u["residual"], prominence=6 * rn, distance=4)
        if len(ridx):
            A("   殘差中仍存在的峰（代表有參考譜沒涵蓋到的成分）：")
            order = np.argsort(rprops["prominences"])[::-1][:12]
            for j in sorted(order, key=lambda k: u["x"][ridx[k]]):
                A(f"      {u['x'][ridx[j]]:8.1f} cm⁻¹   S/N="
                  f"{rprops['prominences'][j] / rn:5.1f}")
        else:
            A("   殘差中沒有明顯的峰 → 參考譜已能解釋這張光譜。")

    report = "\n".join(lines)

    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, f"{name}_report.txt"), "w",
              encoding="utf-8") as f:
        f.write(report + "\n")
    with open(os.path.join(outdir, f"{name}_peaks.csv"), "w", newline="",
              encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["峰位_cm-1", "高度", "突出高度", "訊噪比SN"])
        for pk in peaks:
            w.writerow([f"{pk['position']:.1f}", f"{pk['height']:.0f}",
                        f"{pk['prominence']:.0f}", f"{pk['snr']:.1f}"])
    if make_plot:
        plot_report(x, y, yc, base, peaks, top, name,
                    os.path.join(outdir, f"{name}_report.png"), u)
    return report


# ============================================================
# 第 8 部分：命令列介面
# ============================================================
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="顏料拉曼光譜分析流程（教學版）",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", nargs="?", help="光譜檔或資料夾")
    ap.add_argument("-o", "--outdir", default="raman_out", help="輸出資料夾")
    ap.add_argument("--refs", nargs="*", default=None, help="做 NNLS 解混用的參考譜檔")
    ap.add_argument("--snr", type=float, default=6.0, help="找峰的訊噪比門檻（預設 6）")
    ap.add_argument("--lam", type=float, default=1e5, help="ALS 基線硬度（預設 1e5）")
    ap.add_argument("--no-plot", action="store_true", help="不出圖，只出文字與 CSV")
    ap.add_argument("--list-library", action="store_true", help="列出內建譜庫後結束")
    a = ap.parse_args(argv)

    if a.list_library:
        print(f"內建譜庫共 {len(LIBRARY)} 種：\n")
        for k, v in LIBRARY.items():
            print(f"  {k}")
            print(f"      帶位：{v['bands']}")
            print(f"      關鍵帶：{v['key']}")
            print(f"      說明：{v['note']}\n")
        return 0

    if not a.input:
        ap.print_help()
        return 1

    if os.path.isdir(a.input):
        files = sorted(glob.glob(os.path.join(a.input, "*.txt")) +
                       glob.glob(os.path.join(a.input, "*.csv")))
    else:
        files = [a.input]
    if not files:
        print("找不到任何 .txt 或 .csv 檔", file=sys.stderr)
        return 1

    for f in files:
        try:
            print(analyze(f, a.outdir, refs=a.refs, snr=a.snr, lam=a.lam,
                          make_plot=not a.no_plot))
            print()
        except Exception as e:
            print(f"[跳過] {f}：{e}", file=sys.stderr)
    print(f"→ 輸出已寫入：{os.path.abspath(a.outdir)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

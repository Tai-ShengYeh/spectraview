# orange-spectra — spectroscopy widgets for Orange Data Mining

**English** ｜ [中文](#中文說明)

Eleven spectroscopy widgets for [Orange Data Mining](https://orangedatamining.com/).
They share their algorithms and the `.speclib` library format with
[SpectraView](https://github.com/Tai-ShengYeh/spectraview), a desktop
spectroscopy viewer. Fetch spectra from public databases by URL, compare and
search them, build reusable reference libraries, decompose mixtures, and draw
aquaphotomics aquagrams — all inside Orange's visual workflow canvas.

| Widget | What it does |
|---|---|
| **Import Spectrum URL** | Paste an **IRUG id/URL** or a **SOPRANO URL** (direct JCAMP-DX/CSV links, and spectra embedded in Plotly/Highcharts/Chart.js charts, work too), download and plot the spectrum, output an Orange `Table`. |
| **Load Spectra Files** | Bulk-load chosen files, a whole folder or a .zip (no extraction): JCAMP-DX, CSV, matrix CSV, NetCDF `.cdf` → one merged `Table`. |
| **Spectrometer** | Turn a camera + diffraction-grating spectrum photo (Theremino-style) into a calibrated intensity-vs-wavelength spectrum: strip profile + pixel→nm calibration. |
| **Merge Spectra** | Overlay several spectra sources on one plot and output a single combined `Table` (each row a spectrum on a shared grid) — like SpectraView's multi-file overlay. |
| **Spectra Similarity** | Score similarity between two sets of spectra: correlation, cosine, spectral angle (SAM), and Euclidean distance. |
| **Spectral Library** | Build a reference library, save it as **`.speclib`** (interoperable with SpectraView), and rank an unknown spectrum against the library. Offers the **UCL Raman Library of Pigments** (55 pigments, Bell, Clark & Gibbs 1997) as a built-in download — fetched from UCL's own site on first use and cached locally; the data is **not** bundled with the package. |
| **Mixture Analysis** | Decompose a mixed spectrum with non‑negative least squares (NNLS): solve `mixture ≈ Σ cᵢ·refᵢ` and report coefficients, proportions, and R². |
| **Aquagram** | Aquaphotomics: read normalized absorbance at water's 12 characteristic bands (WAMACs) and draw a 12‑axis radar chart (raw / SNV / aquagram normalization). |
| **Peak Finder** | Detect peaks, label them on the plot, and output a peak table (position, height, FWHM, prominence, area). |
| **XRF Element ID** | Find peaks in an XRF spectrum (keV) and label them with matching element emission lines (Kα/Kβ/Lα/Lβ, 53 elements Na–U). |
| **PLS-DA** | Partial least squares discriminant analysis: class-colored score plot, loadings, VIP variable importance, and predictions. |

The output `Table` follows the
[Orange-Spectroscopy](https://orange-spectroscopy.readthedocs.io/) convention
(column names = wavelength/wavenumber, one spectrum per row), so it plugs
straight into the Spectra viewer widget or into PCA / PLS chemometrics
pipelines.

## Install

> ⚠️ Know **which Orange you run** first — the desktop App and a pip‑installed
> Orange are separate Python environments. Installing into the wrong one means
> the widgets won't appear.

Works on **Windows, macOS (Apple Silicon & Intel), and Linux** — it's pure
Python (numpy / scipy / matplotlib / Orange3).

**A. Desktop App** (the standalone program from orangedatamining.com):
`Options ▸ Add-ons… ▸ Add more…`, type **`orange-spectra`**, tick it, **OK**,
and restart. On macOS, download the Mac `.dmg` from
<https://orangedatamining.com/download/>.

**B. pip Orange** (started with `python -m Orange.canvas`):

```bash
pip install orange-spectra
python -m Orange.canvas
```

On macOS use `python3` / `pip3` (e.g. install Python via
[Homebrew](https://brew.sh/): `brew install python`).

If Orange fails to start with `ImportError: PyQt5 … not available`, it's
missing a Qt binding — install one:

```bash
pip install PyQt5 PyQtWebEngine
```

After (re)starting Orange, a **Spectra** category with the eleven widgets appears
in the toolbox.

## Quick start

1. Drop in **Import Spectrum URL**, enter `4119` (IRUG's PB15 phthalocyanine
   blue Raman spectrum) → **Fetch**.
2. Fetch a few reference spectra → feed them to the **Spectral Library**
   *Spectra* input → *Add input spectra to library* → *Save…* as a `.speclib`.
*(or skip fetching: pick the built-in **UCL Raman Library of Pigments**
   in the Library box and press **Add built-in** — 55 reference pigments in
   one click; the data is downloaded from UCL on first use and cached.)*
3. Feed an unknown spectrum to the library's *Query* input → the *Hits* output
   is the ranked match table.
4. Feed a mixed spectrum to **Mixture Analysis** *Mixture* and the references
   (or the library's *Library* output) to *References* → get component
   proportions and R².

Every widget has an **ℹ How to use** box and a **📖 Open tutorial** button in
its top‑left corner, and **F1** (or the ? button) opens its online help page.

## Built-in UCL pigment library

The Spectral Library widget offers the **UCL Raman Library of Pigments**
(55 natural and synthetic pigments in use before ~1850 AD):

> I. M. Bell, R. J. H. Clark and P. J. Gibbs, "Raman spectroscopic library
> of natural and synthetic pigments (pre- ≈ 1850 AD)", *Spectrochim. Acta A*
> **53** (1997) 2159–2179. doi:10.1016/S1386-1425(97)00140-6

Please cite the paper when you use these spectra in published work.

**The data is not redistributed with this package.** On first use the
widget downloads the spectra from UCL's own website (falling back to the
Internet Archive's copy if UCL's server is unreachable from your region)
and caches them at `~/.orange-spectra/ucl_pigments_raman.speclib`; after
that it works offline. Set the `ORANGE_SPECTRA_CACHE` environment variable
to move the cache, or `ORANGE_SPECTRA_UCL_BASE` to point at a mirror. If
you already have the `.spc` files, build the cache offline with
`orangespectra.core.build_ucl_library_from_folder(folder)`.

Notes on the data: intensities are raw counts (no baseline or instrument
response correction — baseline-correct your query the way you normally
would before searching), measurement ranges differ per pigment, and UCL
notes the downloadable band positions are uncorrected and may differ
slightly from the values tabulated on its site.

## Supported URL formats

- **IRUG** detail pages (jqPlot‑embedded data)
- **SOPRANO** pages (Dygraph‑embedded data)
- **JCAMP-DX** (AFFN plain‑number format)
- Two‑column **CSV/TSV**
- Spectra embedded in **Plotly / Highcharts / Chart.js** interactive charts (generic x/y fallback)

Compressed JCAMP (SQZ/DIF) is not parsed here — open it in SpectraView and
export first.

## Documentation

Full tutorial with real‑data demos:
<https://tai-shengyeh.github.io/spectraview/orange.html>
([English](https://tai-shengyeh.github.io/spectraview/orange_en.html)).

## Notes

- You only need to install **once**; updating requires a reinstall
  (`pip install --upgrade orange-spectra`) and an Orange restart.
- Source, issues, and the desktop SpectraView app:
  <https://github.com/Tai-ShengYeh/spectraview>.

## License

MIT.

---

## 中文說明

**orange-spectra — Orange Data Mining 光譜 widgets**

[English](#orange-spectra--spectroscopy-widgets-for-orange-data-mining) ｜ **中文**

十一個給 [Orange Data Mining](https://orangedatamining.com/) 的光譜學 widgets，
與桌面版光譜檢視程式
[SpectraView](https://github.com/Tai-ShengYeh/spectraview) 共用演算法與
`.speclib` 光譜庫格式。可用網址從公開資料庫抓光譜、比對與搜尋、建立可重複使用的
參考光譜庫、拆解混合光譜，並繪製 aquaphotomics 雷達圖——全部在 Orange 的視覺化
工作流程畫布中完成。

| Widget | 功能 |
|---|---|
| **Import Spectrum URL** | 貼上 **IRUG 編號／網址**或 **SOPRANO 網址**（也支援 JCAMP-DX/CSV 直接網址，以及 Plotly/Highcharts/Chart.js 互動圖內嵌資料），下載並畫出光譜，輸出成 Orange `Table`。 |
| **Load Spectra Files** | 批次載入：選檔、整個資料夾或 .zip（免解壓）——JCAMP-DX、CSV、矩陣 CSV、NetCDF `.cdf` 全部讀成一個合併 `Table`。 |
| **Spectrometer** | 把相機＋繞射光柵拍到的光譜照片（Theremino 分光儀）變成校準光譜：取水平帶讀強度、像素→波長校準。 |
| **Merge Spectra** | 把多個光譜來源疊在一張圖，輸出成一個合併 `Table`（每列一條光譜、共同波段）——等同 SpectraView 的多檔疊圖。 |
| **Spectra Similarity** | 兩組光譜間的相似度：correlation / cosine / 光譜角 SAM / Euclidean。 |
| **Spectral Library** | 建立參考光譜庫、存成 **`.speclib`（與 SpectraView 互通）**、對庫比對未知譜並輸出排名。提供 **UCL 顏料拉曼庫**（55 種，Bell, Clark & Gibbs 1997）一鍵下載：第一次使用時自 UCL 官網取得並快取到本機，**套件本身不含該資料**。 |
| **Mixture Analysis** | 混合光譜的成分分析：以非負最小平方（NNLS）解 `mixture ≈ Σ cᵢ·refᵢ`，回報係數、比例與 R²。 |
| **Aquagram** | Aquaphotomics：在水的 12 個特徵吸收帶（WAMACs）取正規化吸光度，畫 12 軸雷達圖（raw / SNV / aquagram 三種正規化）。 |
| **Peak Finder** | 自動尋峰並在圖上標記，輸出峰表（峰位、峰高、FWHM、顯著度、面積）。 |
| **XRF Element ID** | XRF 能譜（keV）尋峰並比對元素特徵譜線（Kα/Kβ/Lα/Lβ，Na–U 53 元素），圖上直接標元素。 |
| **PLS-DA** | 偏最小平方判別分析：依類別上色的分數圖、loadings、VIP 變數重要性與預測輸出。 |

輸出的 `Table` 採
[Orange-Spectroscopy](https://orange-spectroscopy.readthedocs.io/)
慣例（欄名＝波長/波數、每列一條光譜），可直接接其 Spectra 檢視 widget 或
PCA / PLS 等化學計量學流程。

## 安裝

> ⚠️ 先確認你的 **Orange 是哪一種**——桌面版 App 與 pip 版是不同的 Python
> 環境，裝錯不會出現 widgets。

可在 **Windows、macOS（Apple 晶片與 Intel）、Linux** 執行——純 Python
（numpy / scipy / matplotlib / Orange3）。

**A. 桌面版 App**（orangedatamining.com 下載的獨立程式）：
`Options ▸ Add-ons… ▸ Add more…` 輸入 **`orange-spectra`** → 打勾 → OK → 重啟。
macOS 到 <https://orangedatamining.com/download/> 下載 Mac 版 `.dmg`。

**B. pip 版 Orange**（`python -m Orange.canvas` 啟動）：

```bash
pip install orange-spectra
python -m Orange.canvas
```

macOS 請用 `python3` / `pip3`（可用 [Homebrew](https://brew.sh/) 裝 Python：
`brew install python`）。

若 Orange 開不了並出現 `ImportError: PyQt5 … not available`，是少了 Qt 綁定，
補裝：

```bash
pip install PyQt5 PyQtWebEngine
```

重新啟動 Orange，工具箱會出現 **Spectra** 分類（11 個 widgets）。

## 快速上手

1. 拖出 **Import Spectrum URL**，輸入 `4119`（IRUG 的 PB15 酞菁藍 Raman 譜）
   → **Fetch**。
2. 多抓幾條參考譜 → 接 **Spectral Library** 的 *Spectra* 輸入 →
   *Add input spectra to library* → *Save…* 存成 `.speclib`。
（也可以不用自己抓：在 Library 區選內建的 **UCL 顏料拉曼庫**按
   **Add built-in**，一鍵載入 55 種參考顏料；資料第一次使用時自 UCL
   官網下載並快取。）
3. 未知譜接 Library 的 *Query* 輸入 → *Hits* 輸出就是排名表。
4. 混合譜接 **Mixture Analysis** 的 *Mixture*、參考譜（或 Library 的 *Library*
   輸出）接 *References* → 得成分比例與 R²。

每個 widget 左上角都有「**ℹ 說明 How to use**」盒子與「**📖 開啟線上教學**」按鈕；選取 widget 按 **F1**（或 ? 鈕）會開啟線上說明頁。

## 內建 UCL 顏料拉曼庫

Spectral Library widget 內建 **UCL 顏料拉曼庫**（約 1850 年以前使用的
55 種天然與合成顏料）：

> I. M. Bell, R. J. H. Clark and P. J. Gibbs, "Raman spectroscopic library
> of natural and synthetic pigments (pre- ≈ 1850 AD)", *Spectrochim. Acta A*
> **53** (1997) 2159–2179. doi:10.1016/S1386-1425(97)00140-6

發表時使用到這些光譜請引用上面這篇論文。

**本套件不隨附這批資料。**第一次使用時 widget 會自 UCL 官網下載
（若你所在地區連不上 UCL 伺服器，會自動改抓 Internet Archive 的存檔副本），
快取在 `~/.orange-spectra/ucl_pigments_raman.speclib`，之後完全離線可用。
環境變數 `ORANGE_SPECTRA_CACHE` 可改快取位置、`ORANGE_SPECTRA_UCL_BASE`
可指定鏡像；已有 .spc 檔者可用
`orangespectra.core.build_ucl_library_from_folder(資料夾)` 離線建庫。

資料注意事項：強度為原始 counts（未做基線與儀器響應校正——比對前請照
你平常的流程先對未知譜做基線校正）；各顏料量測範圍不同；UCL 官網註明
可下載檔案的譜帶位置未經校正，可能與其網頁表格數值略有出入。

## 支援的網址格式

- **IRUG** 詳情頁（jqPlot 內嵌資料）
- **SOPRANO** 頁（Dygraph 內嵌資料）
- **JCAMP-DX**（AFFN 純數字格式）
- 兩欄 **CSV/TSV**
- **Plotly / Highcharts / Chart.js** 互動圖內嵌的 x/y 資料（通用後備解析）

壓縮的 JCAMP（SQZ/DIF）這裡不解析——請先用 SpectraView 開啟後匯出。

## 教學文件

含真實數據 Demo 的完整教學：
<https://tai-shengyeh.github.io/spectraview/orange.html>
（[English](https://tai-shengyeh.github.io/spectraview/orange_en.html)）。

## 注意

- 安裝**一次**即可，更新才需重裝（`pip install --upgrade orange-spectra`）
  ＋重開 Orange。
- 原始碼、issues、桌面版 SpectraView：
  <https://github.com/Tai-ShengYeh/spectraview>。

## 授權

MIT。

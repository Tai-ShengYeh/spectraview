# orange-spectra — spectroscopy widgets for Orange Data Mining

**English** ｜ [中文](#中文說明)

Five spectroscopy widgets for [Orange Data Mining](https://orangedatamining.com/).
They share their algorithms and the `.speclib` library format with
[SpectraView](https://github.com/Tai-ShengYeh/spectraview), a desktop
spectroscopy viewer. Fetch spectra from public databases by URL, compare and
search them, build reusable reference libraries, decompose mixtures, and draw
aquaphotomics aquagrams — all inside Orange's visual workflow canvas.

| Widget | What it does |
|---|---|
| **Import Spectrum URL** | Paste an **IRUG id/URL** or a **SOPRANO URL** (direct JCAMP-DX/CSV links work too), download and plot the spectrum, and output an Orange `Table` (one spectrum per row, wavenumbers as column names). |
| **Spectra Similarity** | Score similarity between two sets of spectra: correlation, cosine, spectral angle (SAM), and Euclidean distance. |
| **Spectral Library** | Build a reference library, save it as **`.speclib`** (interoperable with SpectraView), and rank an unknown spectrum against the library. |
| **Mixture Analysis** | Decompose a mixed spectrum with non‑negative least squares (NNLS): solve `mixture ≈ Σ cᵢ·refᵢ` and report coefficients, proportions, and R². |
| **Aquagram** | Aquaphotomics: read normalized absorbance at water's 12 characteristic bands (WAMACs) and draw a 12‑axis radar chart (raw / SNV / aquagram normalization). |

The output `Table` follows the
[Orange-Spectroscopy](https://orange-spectroscopy.readthedocs.io/) convention
(column names = wavelength/wavenumber, one spectrum per row), so it plugs
straight into the Spectra viewer widget or into PCA / PLS chemometrics
pipelines.

## Install

> ⚠️ Know **which Orange you run** first — the desktop App and a pip‑installed
> Orange are separate Python environments. Installing into the wrong one means
> the widgets won't appear.

**A. Desktop App** (the standalone program from orangedatamining.com):
`Options ▸ Add-ons… ▸ Add more…`, type **`orange-spectra`**, tick it, **OK**,
and restart.

**B. pip Orange** (started with `python -m Orange.canvas`):

```bash
pip install orange-spectra
python -m Orange.canvas
```

If Orange fails to start with `ImportError: PyQt5 … not available`, it's
missing a Qt binding — install one:

```bash
pip install PyQt5 PyQtWebEngine
```

After (re)starting Orange, a **Spectra** category with the five widgets appears
in the toolbox.

## Quick start

1. Drop in **Import Spectrum URL**, enter `4119` (IRUG's PB15 phthalocyanine
   blue Raman spectrum) → **Fetch**.
2. Fetch a few reference spectra → feed them to the **Spectral Library**
   *Spectra* input → *Add input spectra to library* → *Save…* as a `.speclib`.
3. Feed an unknown spectrum to the library's *Query* input → the *Hits* output
   is the ranked match table.
4. Feed a mixed spectrum to **Mixture Analysis** *Mixture* and the references
   (or the library's *Library* output) to *References* → get component
   proportions and R².

Every widget has an **ℹ How to use** box and a **📖 Open tutorial** button in
its top‑left corner.

## Supported URL formats

- **IRUG** detail pages (jqPlot‑embedded data)
- **SOPRANO** pages (Dygraph‑embedded data)
- **JCAMP-DX** (AFFN plain‑number format)
- Two‑column **CSV/TSV**

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

<a name="中文說明"></a>

# orange-spectra — Orange Data Mining 光譜 widgets

[English](#orange-spectra--spectroscopy-widgets-for-orange-data-mining) ｜ **中文**

五個給 [Orange Data Mining](https://orangedatamining.com/) 的光譜學 widgets，
與桌面版光譜檢視程式
[SpectraView](https://github.com/Tai-ShengYeh/spectraview) 共用演算法與
`.speclib` 光譜庫格式。可用網址從公開資料庫抓光譜、比對與搜尋、建立可重複使用的
參考光譜庫、拆解混合光譜，並繪製 aquaphotomics 雷達圖——全部在 Orange 的視覺化
工作流程畫布中完成。

| Widget | 功能 |
|---|---|
| **Import Spectrum URL** | 貼上 **IRUG 編號／網址**或 **SOPRANO 網址**（也支援 JCAMP-DX/CSV 直接網址），下載並畫出光譜，輸出成 Orange `Table`（每列一條光譜、欄名＝波數）。 |
| **Spectra Similarity** | 兩組光譜間的相似度：correlation / cosine / 光譜角 SAM / Euclidean。 |
| **Spectral Library** | 建立參考光譜庫、存成 **`.speclib`（與 SpectraView 互通）**、對庫比對未知譜並輸出排名。 |
| **Mixture Analysis** | 混合光譜的成分分析：以非負最小平方（NNLS）解 `mixture ≈ Σ cᵢ·refᵢ`，回報係數、比例與 R²。 |
| **Aquagram** | Aquaphotomics：在水的 12 個特徵吸收帶（WAMACs）取正規化吸光度，畫 12 軸雷達圖（raw / SNV / aquagram 三種正規化）。 |

輸出的 `Table` 採
[Orange-Spectroscopy](https://orange-spectroscopy.readthedocs.io/)
慣例（欄名＝波長/波數、每列一條光譜），可直接接其 Spectra 檢視 widget 或
PCA / PLS 等化學計量學流程。

## 安裝

> ⚠️ 先確認你的 **Orange 是哪一種**——桌面版 App 與 pip 版是不同的 Python
> 環境，裝錯不會出現 widgets。

**A. 桌面版 App**（orangedatamining.com 下載的獨立程式）：
`Options ▸ Add-ons… ▸ Add more…` 輸入 **`orange-spectra`** → 打勾 → OK → 重啟。

**B. pip 版 Orange**（`python -m Orange.canvas` 啟動）：

```bash
pip install orange-spectra
python -m Orange.canvas
```

若 Orange 開不了並出現 `ImportError: PyQt5 … not available`，是少了 Qt 綁定，
補裝：

```bash
pip install PyQt5 PyQtWebEngine
```

重新啟動 Orange，工具箱會出現 **Spectra** 分類（5 個 widgets）。

## 快速上手

1. 拖出 **Import Spectrum URL**，輸入 `4119`（IRUG 的 PB15 酞菁藍 Raman 譜）
   → **Fetch**。
2. 多抓幾條參考譜 → 接 **Spectral Library** 的 *Spectra* 輸入 →
   *Add input spectra to library* → *Save…* 存成 `.speclib`。
3. 未知譜接 Library 的 *Query* 輸入 → *Hits* 輸出就是排名表。
4. 混合譜接 **Mixture Analysis** 的 *Mixture*、參考譜（或 Library 的 *Library*
   輸出）接 *References* → 得成分比例與 R²。

每個 widget 左上角都有「**ℹ 說明 How to use**」盒子與「**📖 開啟線上教學**」按鈕。

## 支援的網址格式

- **IRUG** 詳情頁（jqPlot 內嵌資料）
- **SOPRANO** 頁（Dygraph 內嵌資料）
- **JCAMP-DX**（AFFN 純數字格式）
- 兩欄 **CSV/TSV**

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

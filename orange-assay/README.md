# orange-assay — plate-image & dose-response widgets for Orange Data Mining

**English** ｜ [中文](#中文說明)

Three [Orange Data Mining](https://orangedatamining.com/) widgets that turn
**assay plate photos** into tables and fit **dose-response curves** — without
writing code. The image core reuses the
[coffee-ring-analyzer](https://github.com/Tai-ShengYeh/coffee-ring-analyzer)
algorithm, so results match its Streamlit / mobile / CLI versions. Sibling
package to [orange-spectra](https://pypi.org/project/orange-spectra/).

| Widget | What it does |
|---|---|
| **Microplate Reader** | Read a microtitre-plate photo (6–384 wells) into a per-well table: mean intensity, or thresholded count / intensity, with well labels (A1…H12). |
| **Coffee-Ring Reader** | Quantify a spotted-assay (coffee-ring) plate: per-cell area/intensity over a grid, blank-column (col 0) normalization and replicate averaging → a `concentration / ratio` table. |
| **Dose-Response Fit** | Fit a 3- or 4-parameter logistic curve to any `concentration + signal` table and report **EC50**, **LOD** and **R²**, with the fitted curve overlaid on the points. |

The outputs are ordinary Orange `Table`s, so they plug straight into Data
Table, scatter plots, PCA, or the classification widgets — e.g. the
**Dose-Response Fit** widget accepts either reader's output, or any colorimetric
assay you already have in a spreadsheet.

## Pipeline

```
Microplate Reader ─wells──┐
                          ├──► Dose-Response Fit ──► EC50 / LOD / curve
Coffee-Ring Reader ─cells─┘        (also: PCA · PLS-DA · Data Table)
```

## Install

Works on **Windows, macOS (Apple Silicon & Intel), and Linux** — pure Python
(numpy / scipy / matplotlib / Pillow / Orange3).

**A. Desktop App** (from orangedatamining.com):
`Options ▸ Add-ons… ▸ Add more…`, type **`orange-assay`**, tick it, **OK**,
restart.

**B. pip Orange**:

```bash
pip install orange-assay
python -m Orange.canvas
```

After restarting Orange, an **Assay** category with the three widgets appears in
the toolbox. Each widget has an **ℹ How to use** box, and **F1** opens its
online help page.

> **Traditional-Chinese Windows note:** if the Add-ons dialog crashes with
> `UnicodeDecodeError: 'cp950' codec…`, install from a command prompt instead:
> `"…\Programs\Orange\python.exe" -m pip install orange-assay`, then restart
> Orange (or set the `PYTHONUTF8=1` environment variable).

## Quick start

1. **Microplate Reader** → choose a plate photo → pick the format (e.g. 96) →
   *Wells* output → Data Table.
2. Or **Coffee-Ring Reader** → choose a spotted-plate photo → set rows × columns
   (column 0 = blank) → enter the concentrations → *Cells* output.
3. Feed either output to **Dose-Response Fit** → pick the concentration and
   signal columns → read EC50 / LOD / R² and the fitted curve.

## Tests

```bash
python tests/test_core.py                                   # pure logic
QT_QPA_PLATFORM=offscreen python tests/test_widgets.py      # widget smoke tests
```

## License

MIT.

---

<a name="中文說明"></a>

# orange-assay — Orange Data Mining 盤面影像與劑量反應 widgets

[English](#orange-assay--plate-image--dose-response-widgets-for-orange-data-mining) ｜ **中文**

三個 [Orange Data Mining](https://orangedatamining.com/) widgets，把**檢測盤照片**
變成資料表並擬合**劑量反應曲線**——全程拖拉、免寫程式。影像核心重用
[coffee-ring-analyzer](https://github.com/Tai-ShengYeh/coffee-ring-analyzer)
的演算法，結果與其 Streamlit／手機／CLI 版本一致。與
[orange-spectra](https://pypi.org/project/orange-spectra/) 為姊妹套件。

| Widget | 功能 |
|---|---|
| **Microplate Reader** | 讀微孔盤照片（6–384 孔）成每孔一列的表：平均亮度，或閾值面積／強度，附孔位標籤（A1…H12）。 |
| **Coffee-Ring Reader** | 量化點樣（咖啡環）盤：格點逐格取面積／強度、除以同列 blank（第 0 欄）、平均重複列 → `濃度／比值` 表。 |
| **Dose-Response Fit** | 對任何 `濃度＋訊號` 表擬合 3 或 4 參數邏輯曲線，回報 **EC50**、**LOD**、**R²**，並把擬合曲線疊在資料點上。 |

輸出都是一般的 Orange `Table`，可直接接 Data Table、散佈圖、PCA 或分類 widget；
**Dose-Response Fit** 接受兩個 reader 的輸出，也吃你已在試算表裡的比色法數據。

## 工作流程

```
Microplate Reader ─wells──┐
                          ├──► Dose-Response Fit ──► EC50 / LOD / 曲線
Coffee-Ring Reader ─cells─┘        （也可接 PCA · PLS-DA · Data Table）
```

## 安裝

可在 **Windows、macOS（Apple 晶片與 Intel）、Linux** 執行——純 Python
（numpy / scipy / matplotlib / Pillow / Orange3）。

**A. 桌面版 App**：`Options ▸ Add-ons… ▸ Add more…` 輸入 **`orange-assay`** →
打勾 → OK → 重啟。

**B. pip 版 Orange**：

```bash
pip install orange-assay
python -m Orange.canvas
```

重啟後工具箱會出現 **Assay** 分類（3 個 widgets）。每個 widget 都有
「**ℹ 說明**」盒子，按 **F1** 會開線上說明頁。

> **繁體中文 Windows 注意：** 若 Add-ons 對話框出現
> `UnicodeDecodeError: 'cp950' codec…` 崩潰，改用命令列安裝：
> `"…\Programs\Orange\python.exe" -m pip install orange-assay`，再重開 Orange
> （或設環境變數 `PYTHONUTF8=1`）。

## 快速上手

1. **Microplate Reader** → 選盤面照片 → 選規格（如 96）→ *Wells* 輸出 → Data Table。
2. 或 **Coffee-Ring Reader** → 選點樣盤照片 → 設列×欄（第 0 欄＝blank）→ 填濃度
   → *Cells* 輸出。
3. 任一輸出接 **Dose-Response Fit** → 選濃度與訊號欄 → 讀 EC50／LOD／R² 與曲線。

## 測試

```bash
python tests/test_core.py                                   # 純邏輯
QT_QPA_PLATFORM=offscreen python tests/test_widgets.py      # widget 煙霧測試
```

## 授權

MIT。

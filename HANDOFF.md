# Session Handoff — SpectraView / IRUG 光譜 + Orange widgets

> 給 Hermes agent 的接手文件。這份記錄了本 session 從「下載 IRUG 光譜」一路做到
> 「Orange Data Mining 光譜 widgets 套件 + 教學網站 + Pages 部署修復」的完整脈絡、
> 目前狀態、以及**還沒完成、需要接手的事項**。

Repo: `Tai-ShengYeh/spectraview`（GitHub）
工作分支：`claude/spectra-download-plot-nkg1ot`
另一個 repo（僅第一階段用到）：`Tai-ShengYeh/Tai-ShengYeh`

---

## 0. TL;DR — 現在最重要的兩件事（ACTION ITEMS）

1. **PR #3 尚未合併**（狀態 clean、CI 全綠，可隨時合併）：
   https://github.com/Tai-ShengYeh/spectraview/pull/3
   內容 = 第 5 個 widget「Aquagram」+ GitHub Pages 部署修復。

2. **GitHub Pages 需要一次手動設定**（只有 repo owner 能在網頁改，agent/API 無權限）：
   `Settings ▸ Pages ▸ Build and deployment ▸ Source` → 改成 **「GitHub Actions」**。
   - 不改這個設定，教學頁 `orange.html` 會一直 404（見 §4 根因）。
   - 改好 + 合併 PR #3 後，教學頁上線於
     **https://tai-shengyeh.github.io/spectraview/orange.html**

---

## 1. 這個 session 做了什麼（時間順序）

### 階段 A — 下載 IRUG 光譜（起點）
- 需求：`http://www.irug.org/jcamp-details?id=3537` 下載成 CSV 並繪圖。
- 產出：`Tai-ShengYeh/Tai-ShengYeh` repo 的 `irug_spectrum.py`（獨立 script，
  含自寫 JCAMP-DX 解析：AFFN/PAC/SQZ/DIF/DUP）。**此 script 尚未開 PR**（見 §6）。

### 階段 B — 把「線上匯入」做進 SpectraView 桌面版
- **關鍵發現**（使用者提供 R 的 rvest script 才確認）：IRUG 詳情頁**沒有**下載檔，
  光譜資料寫在頁面互動圖 **jqPlot** 的 `<script>` 裡，是一串引號包住的
  `"波數:強度"` 配對（`jqPlotData.series`）。
- 產出：`specview/formats/online.py` + `File ▸ Import from URL / IRUG…` 選單。
- **PR #1 已合併**。

### 階段 C — Orange Data Mining 光譜 widgets 套件（本 session 主體）
- 新增 `orange-spectra/` 附加元件 + 教學網頁 `docs/orange.html`。
- **PR #2 已合併**（4 個 widgets）。

### 階段 D — 第 5 個 widget（Aquagram）+ Pages 修復
- 使用者要求依 aquaphotomics.com / nirpyresearch 做水光譜學 aquagram widget。
- 同時發現並修復 GitHub Pages 部署失敗（見 §4）。
- **PR #3 開著、未合併**（← 就是要接手的）。

---

## 2. `orange-spectra/` 套件結構與內容

安裝：`pip install "git+https://github.com/Tai-ShengYeh/spectraview.git#subdirectory=orange-spectra"`
安裝後 Orange 工具箱出現 **Spectra** 分類，含 5 個 widgets。

```
orange-spectra/
├── pyproject.toml              # 套件定義 + Orange entry-points（deps: numpy scipy matplotlib Orange3）
├── README.md
├── orangespectra/
│   ├── __init__.py
│   ├── core.py                 # 純 numpy/scipy 邏輯（與 Orange/Qt 解耦，可離線測）
│   ├── table_io.py             # Orange Table ↔ spectrum-dict（欄名=波數、列=光譜）
│   └── widgets/
│       ├── __init__.py         # category "Spectra"
│       ├── owimporturl.py      # Import Spectrum URL
│       ├── owsimilarity.py     # Spectra Similarity
│       ├── owlibrary.py        # Spectral Library（.speclib，與 SpectraView 互通）
│       ├── owmixture.py        # Mixture Analysis（NNLS）
│       ├── owaquagram.py       # Aquagram（本階段新增，PR #3）
│       └── icons/*.svg
└── tests/
    ├── test_core.py            # 45 tests（純邏輯，免 GUI）
    └── test_widgets.py         # 36 tests（Orange 官方 WidgetTest，offscreen）
```

### 五個 widgets 功能
| Widget | 輸入 → 輸出 | 說明 |
|---|---|---|
| Import Spectrum URL | (URL) → Spectra Table | IRUG 編號/頁、SOPRANO 頁、JCAMP-DX(AFFN)、CSV；即時繪圖 |
| Spectra Similarity | Data (+References) → Scores | correlation / cosine / SAM / Euclidean |
| Spectral Library | Spectra, Query → Hits/BestMatch/Library | 建/存/讀 `.speclib`、排名搜尋 |
| Mixture Analysis | Mixture, References → Composition/Fit | NNLS 解 mixture≈Σcᵢ·refᵢ，係數/比例/R² |
| **Aquagram** | Data → Aquagram Coordinates (n×12) | 水光譜學雷達圖；raw/snv/aquagram 三正規化 |

### core.py 重點函式（給接手者）
- `load_spectrum_url(id_or_url, fetch=default_fetch)` — HTTP 抓取可注入 → 離線測。
- `parse_irug_jqplot / parse_soprano / parse_jcamp / parse_csv` — 各來源解析。
- `similarity_scores(xa,ya,xb,yb)` — 四指標，先取重疊波段內插。
- `save_library / load_library / search_library` — `.speclib` = SpectraView 相同 JSON 格式。
- `mixture_nnls(mixture, references, fit_offset=True)`。
- `aquagram_coordinates(spectra, wamacs=None, normalization="aquagram")` —
  回 {wamacs, names, values(n×12), normalization, covered}。
  - `WAMACS` = 12 個標準水吸收帶 (nm)：1342,1364,1372,1382,1398,1410,1438,1444,1464,1474,1492,1516。
  - normalization：`raw`（原值）/`snv`（各譜 SNV）/`aquagram`（SNV＋跨樣品標準化，0=組平均）。

---

## 3. 教學網站

- `docs/index.html` — 既有 SpectraView 桌面版教學頁（使用者原有）。
- `docs/orange.html` — **本 session 新增**，Orange widgets 教學（雙語、含 5 widgets、
  安裝、逐 widget、工作流圖、FAQ）。PR #3 已加入 Aquagram 段落。
- Pages 從 `/docs` 資料夾發佈，所以網址是 `…/spectraview/orange.html`（不含 `/docs/`）。

---

## 4. GitHub Pages 部署問題（重要根因，接手者必讀）

**症狀**：`https://tai-shengyeh.github.io/spectraview/orange.html` 回 404。

**根因（已查證，非本專案檔案問題）**：
- `docs/orange.html` 確實在 main（PR #2 已合併）。
- 舊的分支式 Pages 管線（`pages-build-deployment`）**build 成功且已含 orange.html**，
  但固定卡在「Deploy to GitHub Pages」步驟 → `Timeout reached, aborting!`。
- 部署歷史：run 23 (`0d11115`) ✅ → run 24 (`2cdd392` 使用者的 "Move web app to /app" commit) ❌
  → run 25（合併 PR #2）❌（手動 rerun 仍 ❌）。即**自 `2cdd392` 起持續失敗**。

**修法（PR #3 內）**：
- 新增 `.github/workflows/pages.yml`：改用 `actions/upload-pages-artifact` +
  `actions/deploy-pages@v4` 直接發佈 `docs/`。
- 新增 `docs/.nojekyll`（純靜態、免 Jekyll）。
- **前提**：repo `Settings ▸ Pages ▸ Source` 必須設為「GitHub Actions」（見 §0 第 2 點）。
  - agent 無法改此設定（REST `PUT /repos/.../pages` 需 admin 權限，MCP 未提供工具）。

**驗收**：設定改好 + 合併 PR #3 後，看 Actions 的 "Deploy Pages" workflow 是否 success，
再開 orange.html（首次瀏覽器 Ctrl/Cmd+Shift+R 清快取）。

---

## 5. 目前 PR / 分支狀態

| PR | 內容 | 狀態 |
|---|---|---|
| #1 | 桌面版 IRUG/URL 線上匯入 | ✅ 已合併 |
| #2 | orange-spectra 4 widgets + 教學頁 | ✅ 已合併 |
| **#3** | **Aquagram widget + Pages 修復** | **🟡 open, clean, CI 全綠, 未合併** |

- 分支 `claude/spectra-download-plot-nkg1ot` 目前 HEAD = PR #3 內容（commit `daf3119`）。
- main HEAD = `b09bd0c`（PR #2 合併點）。

---

## 6. 尚未完成 / 可接手的事項（NEXT STEPS）

1. **合併 PR #3**（CI 已綠）。
2. **設定 Pages Source = GitHub Actions**（§0 第 2 點）— 這步一定要人工做。
3. 合併後**驗證**：Actions「Deploy Pages」success → orange.html 可開。
4. **本機實測 IRUG 抓取**：sandbox 擋 irug.org，線上抓取只用 fixture 驗證過。
   請本機開 Orange → Import Spectrum URL 輸入 `4119`（IRUG PB15 Raman）確認真的抓得到。
5. **Aquagram 用真實 NIR 資料驗證**：目前用合成光譜測數學正確性；建議用真實水/食品
   NIR 光譜（涵蓋 1300–1600 nm）跑一次，確認雷達圖合理。
6. **（可選）** `Tai-ShengYeh/Tai-ShengYeh` repo 的 `irug_spectrum.py` 獨立 script
   從未開 PR — 決定要保留、開 PR、還是刪除。
7. **（可選）** 桌面版有兩個「從 URL 匯入」入口（IRUG/online 與 SOPRANO），未來可合併成一個。

---

## 7. 本機如何測試

```bash
# 取得分支
git clone https://github.com/Tai-ShengYeh/spectraview.git
cd spectraview && git checkout claude/spectra-download-plot-nkg1ot

# orange-spectra 純邏輯測試（只需 numpy scipy）
pip install numpy scipy
python orange-spectra/tests/test_core.py            # 期望 45 passed

# widget 測試（需 orange3 + PyQt5；offscreen）
pip install orange3 PyQt5 matplotlib
QT_QPA_PLATFORM=offscreen python orange-spectra/tests/test_widgets.py   # 36 tests OK

# 安裝進 Orange 實際使用
pip install -e orange-spectra
python -m Orange.canvas          # 工具箱應出現 Spectra 分類的 5 個 widgets

# SpectraView 主程式測試
pip install -r requirements.txt
python tests/smoke_test.py                          # 期望 181 passed
```

環境雷（本 session 遇過）：
- Orange3 舊相依（serverfiles/baycomp/python-louvain）在某些環境用 setuptools 建置會噴
  `AttributeError: install_layout` → 設 `export SETUPTOOLS_USE_DISTUTILS=stdlib` 再裝。
- AnyQt 需要 PyQt5（本專案 widget 用 AnyQt；PySide6 單獨不夠）。

---

## 8. 測試結果快照（本 session 最後一次跑）
- `orange-spectra/tests/test_core.py` → **45 passed, 0 failed**
- `orange-spectra/tests/test_widgets.py` → **36 tests OK (skipped=3)**
- `tests/smoke_test.py`（主程式）→ **181 passed, 0 failed**
- Orange canvas discovery → Spectra 分類列出全部 5 個 widgets
- PR #3 CI：test(3.10)/test(3.12)/orange-addon-core 全 success

---

_本文件由 Claude Code session 產生，供 Hermes agent 接手。_

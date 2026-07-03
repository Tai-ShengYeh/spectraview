# 發佈 orange-spectra 到 PyPI

上架 PyPI 之後，**Orange 桌面版 App** 的使用者就能用
`Options ▸ Add-ons ▸ Add more…` 打 **`orange-spectra`** 一鍵安裝（免記 git 網址），
任何電腦都一樣。

有兩種發佈方式，**推薦方式 A（免 token、自動化）**。

---

## 方式 A ── GitHub Release 自動發佈（Trusted Publishing，免 token）✅

`.github/workflows/publish-pypi.yml` 已設好，用 PyPI 的 **Trusted Publishing（OIDC）**——
**不需要產生或貼任何 API token**，也不用在 GitHub 存 secret。

### 一次性設定（在 PyPI 網頁，約 3 分鐘）
1. 註冊 [PyPI](https://pypi.org/account/register/) 帳號。
2. 到 **PyPI ▸ 你的帳號 ▸ Publishing**（<https://pypi.org/manage/account/publishing/>）
   → **Add a new pending publisher**，填：
   - **PyPI Project Name**：`orange-spectra`
   - **Owner**：`Tai-ShengYeh`
   - **Repository name**：`spectraview`
   - **Workflow name**：`publish-pypi.yml`
   - **Environment name**：`pypi`
   → 送出。（"pending publisher" 讓你在專案還沒存在時就先設好，第一次發佈就適用。）

### 每次發佈（在 GitHub，一鍵）
1. 先把 `orange-spectra/pyproject.toml` 與 `orangespectra/__init__.py` 的 `version` 調高
   （PyPI 不允許覆蓋同一版本號）。
2. GitHub → repo → **Releases ▸ Draft a new release** → 建一個新 tag（如 `v0.1.2`）→ **Publish release**。
3. 這會觸發 **Publish orange-spectra to PyPI** workflow，自動 build + 上傳到 PyPI。
   到 **Actions** 分頁看它跑完變綠即完成。

> 想先試水溫：Actions ▸ **Publish orange-spectra to PyPI** ▸ **Run workflow** ▸ 選
> `testpypi`，會傳到 [TestPyPI](https://test.pypi.org/)（需另在 TestPyPI 設一個 pending publisher）。

---

## 方式 B ── 本機手動上傳（用 API token）

```bash
cd orange-spectra
python -m pip install --upgrade build twine
rm -rf dist build *.egg-info         # PowerShell: Remove-Item -Recurse -Force dist,build,*.egg-info
python -m build                      # 產生 dist/*.whl 與 *.tar.gz
python -m twine check dist/*         # 應顯示 PASSED
python -m twine upload dist/*        # 帳號填 __token__，密碼貼 API token（pypi-…）
```
API token 於 **PyPI ▸ Account settings ▸ API tokens ▸ Add API token** 建立。

---

## 上架後使用者怎麼裝
- **桌面版 App**：Options ▸ Add-ons… ▸ Add more… ▸ 輸入 `orange-spectra` ▸ 打勾 ▸ OK ▸ 重啟。
- **pip 版 Orange**：`pip install orange-spectra`（更新：`pip install --upgrade orange-spectra`）。

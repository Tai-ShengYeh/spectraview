# 發佈 orange-spectra 到 PyPI（給維護者）

上架 PyPI 之後，**Orange 桌面版 App** 的使用者就能用
`Options ▸ Add-ons ▸ Add more…` 打 **`orange-spectra`** 直接安裝
（Add-ons 對話框只吃 PyPI 名稱，不吃 git 網址）。

## 一次性準備
1. 註冊 [PyPI](https://pypi.org/account/register/) 帳號。
2. 建一個 **API token**：PyPI ▸ Account settings ▸ API tokens ▸ Add API token。

## 每次發佈
```bash
cd orange-spectra
python -m pip install --upgrade build twine
rm -rf dist build *.egg-info          # Windows PowerShell: Remove-Item -Recurse -Force dist,build,*.egg-info
python -m build                       # 產生 dist/*.whl 與 *.tar.gz
python -m twine check dist/*          # 應顯示 PASSED
python -m twine upload dist/*         # 帳號填 __token__，密碼貼 API token（pypi-...）
```

> 想先試水溫可先傳 [TestPyPI](https://test.pypi.org/)：
> `python -m twine upload --repository testpypi dist/*`

## 版本更新
改 `pyproject.toml` 的 `version`（PyPI 不允許覆蓋同一版本號），再跑上面的 build + upload。

## 上架後使用者怎麼裝
- **桌面版 App**：Options ▸ Add-ons ▸ Add more… ▸ 輸入 `orange-spectra` ▸ 打勾 ▸ OK ▸ 重啟。
- **pip 版 Orange**：`pip install orange-spectra`。

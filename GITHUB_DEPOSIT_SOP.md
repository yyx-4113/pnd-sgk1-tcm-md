# GitHub 复现包归档操作手册（中文）

> 适用对象：本仓库 `pnd-sgk1-tcm-md`（PND / SGK1 / 中药虚拟筛选 + 100 ns MD 复现包）
> 原则：**数据可用性声明必须指向真实存在的仓库 + 版本 tag + MANIFEST 校验和**，严禁 "available on request"。

## 0. 前置条件
- 本机已登录 GitHub 账号 `yyx-4113`（用 `gh auth status` 确认）。
- 已安装 `git` 与 `gh` CLI。
- 沙箱推送到 github.com 走 env 代理（如 `http_proxy`/`https_proxy` 已设置）；小仓库（无 >100MB 文件）可直接 `git push`。

## 1. 本地仓库初始化（已在构建阶段完成）
仓库目录已就绪：`PND_proj/pnd-sgk1-tcm-md/`，含 `input/`、`scripts/`、`analysis/`、`docs/` 及本手册。
确认 `MANIFEST` 已生成（见第 3 节），然后：

```bash
cd PND_proj/pnd-sgk1-tcm-md
git init -q
git add -A
git commit -q -m "v1.0.0: SGK1/7PUE TCM screening + 100 ns MD reproducibility package"
```

## 2. 在 GitHub 创建远程仓库（仅首次）
用 gh 创建（默认 private 也可，但数据须公开，建议 `--public`）：

```bash
gh repo create yyx-4113/pnd-sgk1-tcm-md \
  --public --description "SGK1/7PUE TCM virtual screening + 100 ns MD reproducibility package" \
  --source . --remote origin --push --branch main
```

若远程已存在（手动网页创建过），改用手动关联并推送：

```bash
git remote add origin https://github.com/yyx-4113/pnd-sgk1-tcm-md.git
git branch -M main
git push -u origin main
```

## 3. 生成 MANIFEST（sha256 校验和）
在仓库根目录执行（相对路径，便于他人 `sha256sum -c` 校验）：

```bash
cd PND_proj/pnd-sgk1-tcm-md
find . -type f -not -path './.git/*' -not -name 'MANIFEST' \
  | sed 's|^\./||' | sort \
  | while read f; do printf '%s  %s\n' "$(sha256sum "$f" | cut -d' ' -f1)" "$f"; done \
  > MANIFEST
git add MANIFEST
git commit -q -m "add MANIFEST sha256 checksums"
```

> 注意：MANIFEST 自身不参与自校验；提交后再补一次 commit 即可。

## 4. 打版本 tag 并推送
```bash
git tag -a v1.0.0 -m "SGK1/7PUE TCM screening + 100 ns MD v1.0.0"
git push origin main --tags
```

## 5. 稿件数据可用性声明回填核对
稿件 `docs/manuscript_humanized.md` 与 `docs/cover_letter.md` 已写入：
- 仓库 URL：`https://github.com/yyx-4113/pnd-sgk1-tcm-md`
- tag：`v1.1.1`

推送成功后，该声明即成立。若日后补推原始 `.xtc` 大轨迹（如走 Zenodo 归档），须在稿件补充说明归档 DOI，**不得**在尚未创建时提前声称存在。

## 6. 校验他人复现
```bash
# 克隆后
sha256sum -c MANIFEST      # 应全部 OK
```

## 7. 回滚 / 修正
- 修正内容后递增 tag（如 `v1.0.1`），并同步更新稿件声明中的 tag。
- 切勿 force-push 已公开 tag；如需修订，发新 tag 并在稿件注明。

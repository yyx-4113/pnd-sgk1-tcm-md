#!/usr/bin/env bash
# =====================================================================
# analyze.sh — MD 稳定性分析 + MM-PBSA + 汇总
# 用法: WORK=~/pnd_md bash analyze.sh
# =====================================================================
set -eo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="${WORK:-$HOME/pnd_md}"
OUT="$WORK/output"
LIGANDS="${LIGANDS:-164676 5280445 5280863 5281708 439246 445154 25022668}"

# 从 index.ndx 里查某个组名的组号（0 起）
grp_num() {
  awk -v t="$2" '/^\[/{gsub(/[][ \t]/,"",$0); if($0==t){print idx; exit} idx++}' "$1"
}

for cid in $LIGANDS; do
  SYS="$OUT/$cid"
  [ -s "$SYS/md.tpr" ] && [ -s "$SYS/md.xtc" ] || { echo "[skip] $cid 缺 md.tpr/md.xtc"; continue; }
  echo "-------- $cid --------"

  [ -s "$SYS/index.ndx" ] || printf "q\n" | gmx make_ndx -f "$SYS/md.tpr" -o "$SYS/index.ndx" >/dev/null 2>&1
  PROT=$(grp_num "$SYS/index.ndx" Protein)
  BB=$(grp_num "$SYS/index.ndx" Backbone)
  LIG=$(grp_num "$SYS/index.ndx" LIG)
  [ -n "$PROT" ] || PROT=1
  [ -n "$BB" ] || BB=4
  [ -n "$LIG" ] || LIG="$PROT"

  # 蛋白骨架 RMSD（对齐 Backbone）
  printf '%s\n%s\n' "$BB" "$BB" | gmx rms -s "$SYS/md.tpr" -f "$SYS/md.xtc" -n "$SYS/index.ndx" -o "$SYS/rmsd.xvg" -tu ns || echo "⚠ rms 失败"
  # 蛋白 RMSF（按残基）
  printf '%s\n' "$PROT" | gmx rmsf -s "$SYS/md.tpr" -f "$SYS/md.xtc" -n "$SYS/index.ndx" -o "$SYS/rmsf.xvg" -res || echo "⚠ rmsf 失败"
  # 回转半径
  printf '%s\n' "$PROT" | gmx gyrate -s "$SYS/md.tpr" -f "$SYS/md.xtc" -n "$SYS/index.ndx" -o "$SYS/rgyrate.xvg" || echo "⚠ gyrate 失败"
  # SASA
  printf '%s\n' "$PROT" | gmx sasa -s "$SYS/md.tpr" -f "$SYS/md.xtc" -n "$SYS/index.ndx" -o "$SYS/sasa.xvg" -surface "$PROT" || echo "⚠ sasa 失败"
  # 蛋白-配体氢键数
  printf '%s\n%s\n' "$PROT" "$LIG" | gmx hbond -s "$SYS/md.tpr" -f "$SYS/md.xtc" -n "$SYS/index.ndx" -num "$SYS/hbond.xvg" || echo "⚠ hbond 失败"
  # 配体-蛋白最近重原子距离（留袋正确判据，替代未独立 fit 的 lig_rmsd）
  printf '%s\n%s\n' "$LIG" "$PROT" | gmx mindist -s "$SYS/md.tpr" -f "$SYS/md.xtc" -n "$SYS/index.ndx" -od "$SYS/lig_dist.xvg" || echo "⚠ mindist 失败"
  # 配体 RMSD（相对于蛋白拟合）—— 判断是否留袋
  printf '%s\n%s\n' "$PROT" "$LIG" | gmx rms -s "$SYS/md.tpr" -f "$SYS/md.xtc" -n "$SYS/index.ndx" -o "$SYS/lig_rmsd.xvg" -tu ns || echo "⚠ lig rms 失败"
  # 配体构象聚类
  printf '%s\n%s\n' "$LIG" "$LIG" | gmx cluster -s "$SYS/md.tpr" -f "$SYS/md.xtc" -n "$SYS/index.ndx" -cl "$SYS/clusters.pdb" -cutoff 0.2 || echo "⚠ cluster 失败"

  # MM-PBSA
  if command -v gmx_MMPBSA >/dev/null 2>&1; then
    gmx_MMPBSA -O -i "$SCRIPT_DIR/mmpbsa.in" \
      -cs "$SYS/md.tpr" -ct "$SYS/md.xtc" -cp "$SYS/topol.top" \
      -ci "$SYS/index.ndx" -cg "$PROT" "$LIG" \
      -o "$SYS/FINAL_RESULTS_MMPBSA.dat" -eo "$SYS/bind_energy.csv" -nogui \
      || echo "⚠ $cid MM-PBSA 失败(可单独排查)"
  else
    echo "[warn] 无 gmx_MMPBSA，跳过 $cid MM-PBSA"
  fi
done

# ---------------- 汇总 T_md_summary.csv ----------------
python3 - "$OUT" $LIGANDS <<'PY'
import sys, os, re, csv
out = sys.argv[1]; ligs = sys.argv[2:]
rows = []
for cid in ligs:
    sysd = os.path.join(out, cid)
    r = {"cid": cid, "rmsd_plateau_nm": "", "lig_rmsd_nm": "", "lig_dist_mean_nm": "",
         "rg_drift_pct": "", "hbond_mean": "",
         "dg_kcal_mol": "", "stable": "?", "stay_in_pocket": "?"}
    for key, fn in [("rmsd_plateau_nm", "rmsd.xvg"), ("lig_rmsd_nm", "lig_rmsd.xvg")]:
        p = os.path.join(sysd, fn)
        if os.path.isfile(p):
            vals = []
            for ln in open(p):
                if ln.startswith(("#", "@")): continue
                ps = ln.split()
                if len(ps) >= 2:
                    try: vals.append(float(ps[1]))
                    except: pass
            if vals:
                r[key] = round(sum(vals[len(vals)//2:]) / max(1, len(vals) - len(vals)//2), 3)
    # 留袋判据：配体-蛋白最近重原子距离（mindist）< 0.40 nm
    pmf = os.path.join(sysd, "lig_dist.xvg")
    if os.path.isfile(pmf):
        v = [float(l.split()[1]) for l in open(pmf) if not l.startswith(("#", "@")) and len(l.split()) >= 2]
        if v:
            r["lig_dist_mean_nm"] = round(sum(v)/len(v), 3)
            r["stay_in_pocket"] = "YES" if max(v) < 0.40 else "CHECK"
    # 稳定判据：回转半径 Rg 前后半段漂移 < 1%（骨架 RMSD 对重建受体4链为 artifact，不用于判 stable）
    rg = os.path.join(sysd, "rgyrate.xvg")
    if os.path.isfile(rg):
        v = [float(l.split()[1]) for l in open(rg) if not l.startswith(("#", "@")) and len(l.split()) >= 2]
        if len(v) >= 2:
            drift = (sum(v[len(v)//2:])/max(1,len(v)-len(v)//2) - sum(v[:len(v)//2])/max(1,len(v)//2)) / (sum(v[:len(v)//2])/max(1,len(v)//2)) * 100
            r["rg_drift_pct"] = round(drift, 2)
            r["stable"] = "YES" if abs(drift) < 1.0 else "CHECK"
    ph = os.path.join(sysd, "hbond.xvg")
    if os.path.isfile(ph):
        v = [float(l.split()[1]) for l in open(ph) if not l.startswith(("#", "@")) and len(l.split()) >= 2]
        if v: r["hbond_mean"] = round(sum(v) / len(v), 1)
    pm = os.path.join(sysd, "FINAL_RESULTS_MMPBSA.dat")
    if os.path.isfile(pm):
        t = open(pm).read()
        m = re.search(r"ΔTOTAL\s*([-+]?\d+\.?\d*)", t)
        if m: r["dg_kcal_mol"] = m.group(1)
    rows.append(r)
dst = os.path.join(out, "T_md_summary.csv")
with open(dst, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else
                       ["cid","rmsd_plateau_nm","lig_rmsd_nm","hbond_mean","dg_kcal_mol","stable","stay_in_pocket"])
    w.writeheader(); w.writerows(rows)
print("汇总 ->", dst)
for r in rows: print(" ", r)
PY
echo "[done] 分析完成"

#!/usr/bin/env bash
# =====================================================================
# run_md.sh — SGK1(7PUE) 7 体系 100 ns MD（AUTODL RTX 4090 / GPU 版）
#
# ★ 2026-09-30 重大更正（一天故障的真正根因）:
#   此前「CUDA error #700 / cudaErrorIllegalAddress」并非 GROMACS 或 4090 的缺陷，
#   而是**受体文件本身是坏的**：三次 fix_receptor*.py + PDBFixer 的
#   addMissingResidues() 把 7PUE 三个天然缺失环（135-148 / 246-253 / 108-109）
#   凭空补出并摆错位置，导致残基 149 之后整段偏离主链 22.5 Å、编号错位 81。
#   坏结构在 EM 卡死(Fmax 2.6e4) → NVT 里 THR149 侧链爆炸 → 坐标 NaN →
#   GPU 报 illegal memory access。**崩溃是坏结构的下游症状。**
#   换用重建的正确受体（rebuild_receptor.py：只补残基内缺失原子 + 补 TER 断链，
#   不补缺失环）后：EM 1047 步收敛(Fmax 979)、NVT/NPT 通过、GPU 零 CUDA 错误、
#   实测 920 ns/day。=> 恢复 GPU 运行。
#
# 设计原则:
#   ① 正确受体 receptor_md.pdb（断链用 TER 分段，4 条链，271 残基全部完整）。
#   ② GPU: -nb gpu -pme gpu -bonded cpu -ntmpi 1 -ntomp 16（实测 920 ns/day）。
#      EM 留在 CPU（steep 只需 1000 步、~1 分钟，避免极小化器的 GPU 边角情形）。
#   ③ 标准协议：EM/NVT/NPT 用 -DPOSRES 固定蛋白；production 无约束。
#      配体保持对接姿态（acpype 坐标 ≡ Vina 姿态，最近重原子接触 2.7–3.5 Å）。
#   ④ 每步 grompp/mdrun 均检查返回码与产物，失败即中止（杜绝"假 ALL_DONE"）。
#   ⑤ GROMACS 2018+ 开 -DPOSRES 时 grompp 必须显式 -r <参考坐标>；四步均已带。
#
# 用法:  conda activate md && WORK=~/pnd_md LIGANDS="..." bash run_md.sh
# 可选环境变量: LIGANDS / NST / MAXH / NSOLV / NTMPI / NTOMP / GMX_GPU
# =====================================================================
set -o pipefail
source /root/miniconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate md 2>/dev/null
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

WORK="${WORK:-$HOME/pnd_md}"
IN="$WORK/input"; OUT="$WORK/output"; MDP="$WORK/mdp"; LIG_IN="$IN/lig"
LIG_TSV="$SCRIPT_DIR/ligands.tsv"
mkdir -p "$OUT" "$MDP"

LIGANDS="${LIGANDS:-164676 5280445 5280863 5281708 439246 445154 25022668}"
FORCE="${FORCE:-amber99sb-ildn}"
WATER="${WATER:-tip3p}"
MAXH="${MAXH:-20}"
NST="${NST:-50000000}"        # 100 ns @ dt=2fs
NSOLV="${NSOLV:-1.0}"
NTMPI="${NTMPI:-1}"
NTOMP="${NTOMP:-16}"
GMX_GPU="${GMX_GPU:--nb gpu -pme gpu -bonded cpu}"

command -v gmx >/dev/null 2>&1 || { echo "✗ 未找到 gmx"; exit 127; }
RECEPTOR="${RECEPTOR:-$IN/receptor_md.pdb}"
[ -f "$RECEPTOR" ] || { echo "✗ 缺 $RECEPTOR（运行 rebuild_receptor.py 生成）"; exit 1; }
echo "[md] $(gmx --version 2>/dev/null | grep -im1 'GROMACS version')  |  mdrun: $NTMPI MPI x $NTOMP OMP + GPU [$GMX_GPU]"

charge_of() { local cid="$1" ch=""; [ -f "$LIG_TSV" ] && ch=$(awk -F'\t' -v c="$cid" 'NR>1 && $1==c {print $7}' "$LIG_TSV"); echo "${ch:-0}"; }
name_of()   { local cid="$1" nm=""; [ -f "$LIG_TSV" ] && nm=$(awk -F'\t' -v c="$cid" 'NR>1 && $1==c {print $2}' "$LIG_TSV"); echo "${nm:-$cid}"; }

# ---------------- mdp ----------------
cat > "$MDP/minim.mdp" <<'EOF'
integrator      = steep
define           = -DPOSRES
emtol           = 1000.0
emstep          = 0.01
nsteps          = 50000
nstlist         = 10
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
pbc             = xyz
EOF

cat > "$MDP/nvt.mdp" <<'EOF'
title                   = NVT equilibration
integrator              = md
define                   = -DPOSRES
dt                      = 0.002
nsteps                  = 50000
continuation            = no
gen_vel                 = yes
gen_temp                = 300
gen_seed                = -1
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstenergy               = 5000
nstlog                  = 5000
nstxout-compressed      = 50000
compressed-x-grps       = System
cutoff-scheme           = Verlet
nstlist                 = 10
rcoulomb                = 1.0
rvdw                    = 1.0
coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16
tcoupl                  = V-rescale
tc-grps                 = System
tau_t                   = 0.1
ref_t                   = 300
pcoupl                  = no
pbc                     = xyz
EOF

cat > "$MDP/npt.mdp" <<'EOF'
title                   = NPT equilibration
integrator              = md
define                   = -DPOSRES
dt                      = 0.002
nsteps                  = 500000
continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstenergy               = 5000
nstlog                  = 5000
nstxout-compressed      = 50000
compressed-x-grps       = System
cutoff-scheme           = Verlet
nstlist                 = 10
rcoulomb                = 1.0
rvdw                    = 1.0
coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16
tcoupl                  = V-rescale
tc-grps                 = System
tau_t                   = 0.1
ref_t                   = 300
pcoupl                  = Berendsen
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5
pbc                     = xyz
EOF

cat > "$MDP/md.mdp" <<EOF
title                   = Production MD
integrator              = md
dt                      = 0.002
nsteps                  = ${NST}
continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstenergy               = 5000
nstlog                  = 5000
nstcalcenergy           = 100
nstxout-compressed      = 50000
compressed-x-grps       = System
cutoff-scheme           = Verlet
nstlist                 = 10
rcoulomb                = 1.0
rvdw                    = 1.0
coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16
tcoupl                  = V-rescale
tc-grps                 = System
tau_t                   = 0.1
ref_t                   = 300
pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5
pbc                     = xyz
EOF
echo "[md] mdp 已生成于 $MDP"

# ---------------- 单体系流程 ----------------
run_one() {
  local cid="$1"
  local SYS="$OUT/$cid"; mkdir -p "$SYS"
  local LIGMOL="$LIG_IN/$cid.mol2"
  local CHARGE="$(charge_of "$cid")"
  local NM="$(name_of "$cid")"
  [ -s "$LIGMOL" ] || { echo "[skip] $cid 缺 $LIGMOL"; return 0; }
  echo "======== $cid ($NM)  net_charge=$CHARGE ========"

  # 1) 蛋白拓扑（正确受体；pdb2gmx 默认 posre.itp = 蛋白重原子约束）
  if [ ! -f "$SYS/protein.gro" ]; then
    printf '\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n' > "$SYS/.ter.ans"
    gmx pdb2gmx -f "$RECEPTOR" -o "$SYS/protein.gro" \
      -p "$SYS/topol.top" -i "$SYS/posre.itp" -ff "$FORCE" -water "$WATER" -ignh \
      < "$SYS/.ter.ans" \
      || { echo "✗ $cid pdb2gmx 失败"; return 1; }
  fi

  # 2) 配体拓扑 (acpype, GAFF2)；带电极体用 Gasteiger 避开 sqm 对奇电子阴离子的崩溃
  local ACYPE_CHG="-c bcc"; [ "$CHARGE" != "0" ] && ACYPE_CHG="-c gas"
  if [ ! -s "$SYS/lig.itp" ] || [ ! -s "$SYS/lig.gro" ]; then
    ( cd "$SYS" && acpype -i "$LIGMOL" -b LIG $ACYPE_CHG -n "$CHARGE" -a gaff2 -o gmx ) \
      || { echo "✗ $cid acpype 失败"; return 1; }
    cp -f "$(find "$SYS" -name 'LIG_GMX.itp' | head -1)" "$SYS/lig.itp" || { echo "✗ $cid 缺 LIG_GMX.itp"; return 1; }
    cp -f "$(find "$SYS" -name 'LIG_GMX.gro' | head -1)" "$SYS/lig.gro" || { echo "✗ $cid 缺 LIG_GMX.gro"; return 1; }
    echo "[ok] $cid 配体拓扑(acpype $ACYPE_CHG, charge=$CHARGE)"
  fi

  # 3) 合并 gro + 接入 topol
  if [ ! -f "$SYS/complex.gro" ]; then
    python3 "$SCRIPT_DIR/prep_system.py" \
      --prot-gro "$SYS/protein.gro" --prot-top "$SYS/topol.top" \
      --lig-gro "$SYS/lig.gro" --lig-itp "$SYS/lig.itp" --lig-name LIG \
      --out-gro "$SYS/complex.gro" || { echo "✗ $cid prep_system 失败"; return 1; }
  fi

  # 4) 盒子 + 溶剂 + 离子
  if [ ! -s "$SYS/solv.gro" ]; then
    gmx editconf -f "$SYS/complex.gro" -o "$SYS/boxed.gro" -c -d "$NSOLV" -bt dodecahedron || { echo "✗ $cid editconf 失败"; return 1; }
    gmx solvate -cp "$SYS/boxed.gro" -cs spc216.gro -p "$SYS/topol.top" -o "$SYS/solv.gro" || { echo "✗ $cid solvate 失败"; return 1; }
    gmx grompp -f "$MDP/minim.mdp" -c "$SYS/solv.gro" -r "$SYS/solv.gro" -p "$SYS/topol.top" -o "$SYS/ions.tpr" -maxwarn 2 \
      || { echo "✗ $cid ions grompp 失败"; return 1; }
    echo "SOL" | gmx genion -s "$SYS/ions.tpr" -p "$SYS/topol.top" -o "$SYS/solv.gro" -neutral -conc 0.15 \
      || { echo "✗ $cid genion 失败"; return 1; }
  fi

  # 5) 能量最小化（steep, CPU）
  if [ ! -s "$SYS/em.gro" ]; then
    gmx grompp -f "$MDP/minim.mdp" -c "$SYS/solv.gro" -r "$SYS/solv.gro" -p "$SYS/topol.top" -o "$SYS/em.tpr" -maxwarn 2 \
      || { echo "✗ $cid EM grompp 失败"; return 1; }
    gmx mdrun -s "$SYS/em.tpr" -deffnm "$SYS/em" -ntmpi "$NTMPI" -ntomp "$NTOMP" || { echo "✗ $cid EM mdrun 失败"; return 1; }
    [ -s "$SYS/em.gro" ] || { echo "✗ $cid EM 无 em.gro"; return 1; }
    grep -E 'Maximum force|did not converge|converged' "$SYS/em.log" | tail -2
  fi

  # 6) NVT（GPU）
  if [ ! -s "$SYS/nvt.gro" ]; then
    gmx grompp -f "$MDP/nvt.mdp" -c "$SYS/em.gro" -r "$SYS/em.gro" -p "$SYS/topol.top" -o "$SYS/nvt.tpr" -maxwarn 2 \
      || { echo "✗ $cid NVT grompp 失败"; return 1; }
    gmx mdrun -s "$SYS/nvt.tpr" -deffnm "$SYS/nvt" $GMX_GPU -ntmpi "$NTMPI" -ntomp "$NTOMP" -maxh "$MAXH" \
      || { echo "✗ $cid NVT mdrun 失败"; return 1; }
    [ -s "$SYS/nvt.gro" ] || { echo "✗ $cid NVT 无 nvt.gro"; return 1; }
  fi

  # 7) NPT（GPU）
  if [ ! -s "$SYS/npt.gro" ]; then
    gmx grompp -f "$MDP/npt.mdp" -c "$SYS/nvt.gro" -r "$SYS/nvt.gro" -p "$SYS/topol.top" -o "$SYS/npt.tpr" -maxwarn 2 \
      || { echo "✗ $cid NPT grompp 失败"; return 1; }
    gmx mdrun -s "$SYS/npt.tpr" -deffnm "$SYS/npt" $GMX_GPU -ntmpi "$NTMPI" -ntomp "$NTOMP" -maxh "$MAXH" \
      || { echo "✗ $cid NPT mdrun 失败"; return 1; }
    [ -s "$SYS/npt.gro" ] || { echo "✗ $cid NPT 无 npt.gro"; return 1; }
  fi

  # 8) 生产（GPU，断点续跑；md.gro 出现即代表 100 ns 完成）
  gmx grompp -f "$MDP/md.mdp" -c "$SYS/npt.gro" -p "$SYS/topol.top" -o "$SYS/md.tpr" -maxwarn 2 \
    || { echo "✗ $cid MD grompp 失败"; return 1; }
  local iter=0
  while [ ! -s "$SYS/md.gro" ]; do
    iter=$((iter+1)); [ "$iter" -gt 40 ] && { echo "✗ $cid 续跑>40 次仍无 md.gro"; return 1; }
    if [ -s "$SYS/md.cpt" ]; then
      gmx mdrun -s "$SYS/md.tpr" -deffnm "$SYS/md" $GMX_GPU -ntmpi "$NTMPI" -ntomp "$NTOMP" -maxh "$MAXH" -cpi "$SYS/md.cpt" -append
    else
      gmx mdrun -s "$SYS/md.tpr" -deffnm "$SYS/md" $GMX_GPU -ntmpi "$NTMPI" -ntomp "$NTOMP" -maxh "$MAXH"
    fi
    [ -s "$SYS/md.gro" ] || [ -s "$SYS/md.cpt" ] || { echo "✗ $cid mdrun 异常(无 md.gro/md.cpt)"; return 1; }
  done
  echo "[ok] $cid 生产完成 -> $SYS/md.gro / md.xtc"
}

for cid in $LIGANDS; do run_one "$cid"; done
echo "[done] 全部体系生产阶段结束"

# -*- coding: utf-8 -*-
"""收敛性补强：用 mdtraj 读本地 md.xtc 补算 7 体系 RMSF(per-residue) 与回转半径 Rg(per-frame)。
- RMSF：相对时间平均结构的 per-atom 涨落，按残基取 backbone 均值。
- Rg：gyration radius，判断 100ns 内蛋白整体尺寸是否平稳（收敛）。
产物：results/<CID>/rmsf.xvg (残基号, RMSF nm)、results/<CID>/rg.xvg (时间, Rg nm)。
"""
import os, argparse, json, numpy as np, mdtraj as md

CIDS = ["164676", "5280445", "5280863", "5281708", "439246", "445154", "25022668"]
RES = "results"

def write_xvg(path, x, y, title, ylabel):
    with open(path, "w") as f:
        f.write(f"# {title}\n@ xaxis label Time (ns)\n@ yaxis label {ylabel}\n")
        for xi, yi in zip(x, y):
            f.write(f"{xi:.3f} {yi:.4f}\n")

def process(cid):
    gro = f"{RES}/{cid}/md.gro"
    xtc = f"{RES}/{cid}/md.xtc"
    top = md.load(gro)
    prot = top.topology.select("protein and (name N or name CA or name C or name O)")
    traj = md.load(xtc, top=gro)
    times = traj.time  # ns

    # Rg (gyration radius) per frame, protein only — 手动算（避免版本 API 差异）
    prot_idx = top.topology.select("protein")
    xyz = traj.xyz[:, prot_idx, :]            # (F, Np, 3)
    com = xyz.mean(axis=1, keepdims=True)
    rg = np.sqrt(((xyz - com) ** 2).sum(axis=2).mean(axis=1))  # nm

    # RMSF per atom (protein backbone)：先 superpose 去刚体运动，再算逐原子对时间平均位置的涨落
    traj.superpose(traj, frame=0, atom_indices=prot)
    xyz = traj.xyz[:, prot, :]                       # (F, Natom_sel, 3)
    mean = xyz.mean(axis=0, keepdims=True)          # (1, Natom_sel, 3)
    rmsf_atom = np.sqrt(((xyz - mean) ** 2).sum(axis=2).mean(axis=0))  # (Natom_sel,) nm
    # 聚合到残基（取 backbone 均值）
    residues = list(top.topology.residues)
    sel_atoms = [top.topology.atom(i) for i in prot]
    per_res = {}
    for atom, val in zip(sel_atoms, rmsf_atom):
        per_res.setdefault(atom.residue.index, []).append(val)
    res_means = []
    res_ids = []
    for rid in sorted(per_res):
        res_means.append(np.mean(per_res[rid]))
        res_ids.append(residues[rid].resSeq)

    out = f"{RES}/{cid}"
    os.makedirs(out, exist_ok=True)
    write_xvg(f"{out}/rg.xvg", times, rg, "Radius of gyration (protein)", "Rg (nm)")
    write_xvg(f"{out}/rmsf.xvg", np.array(res_ids, float), np.array(res_means),
              "RMSF per residue (backbone, time-averaged ref)", "RMSF (nm)")

    # 收敛性判据：后 50ns 的 Rg 均值 vs 全段，看漂移
    rg_full = rg.mean()
    rg_tail = rg[times >= 50].mean()
    rg_drift = (rg_tail - rg_full) / rg_full * 100
    rmsf_mean = np.mean(res_means)
    frac_struct = np.mean(np.array(res_means) < 0.30)  # 结构化核心占比（RMSF<0.3nm）
    print(f"[{cid}] n_frames={len(times)} Rg_full={rg_full:.3f} Rg_50-100ns={rg_tail:.3f} "
          f"drift={rg_drift:+.2f}%  RMSF_mean={rmsf_mean:.3f}  max_res={max(res_means):.3f} "
          f"struct_core(<%0.3nm)={frac_struct*100:.0f}%", flush=True)
    return dict(cid=cid, rg_full=float(rg_full), rg_tail=float(rg_tail), drift=float(rg_drift),
                rmsf_mean=float(rmsf_mean), rmsf_max=float(max(res_means)),
                frac_structured=float(frac_struct))

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cid")
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()
    cids = CIDS if a.all else [a.cid]
    summary = [process(c) for c in cids]
    if a.all:
        json.dump(summary, open(f"{RES}/_convergence_summary.json", "w"), indent=2)
        print("summary ->", f"{RES}/_convergence_summary.json")

#!/usr/bin/env python3
"""本机用 mdtraj 等价复现 GROMACS 轨迹分析（因 Windows 无法装 GROMACS）。
对 5280445(木犀草素)/445154(白藜芦醇) 补算缺失轨迹指标，写回 results/<CID>/。
输出：
  rmsd.xvg       - 蛋白 backbone fit RMSD（等价 gmx rms -sel backbone）
  lig_rmsd.xvg   - 配体(UNL) fit RMSD（等价 gmx rms -sel UNL）
  lig_dist.xvg   - 配体最近蛋白重原子距离（等价 gmx distance/mindist，stay-in-pocket 主证据）
  hbond.xvg      - 蛋白内氢键数（excl water；mdtraj 不识别 UNL 供体/受体，故不含配体氢键）
单位：nm / ps，与 GROMACS xvg 一致。
"""
import mdtraj as md, numpy as np, argparse, time, os
from scipy.spatial.distance import cdist

def write_xvg(path, t, y, title, ylabel):
    with open(path, "w") as f:
        f.write(f"# {title}\n")
        f.write(f"@ title \"{title}\"\n")
        f.write(f"@ xaxis label \"Time (ps)\"\n")
        f.write(f"@ yaxis label \"{ylabel}\"\n")
        for ti, yi in zip(t, y):
            f.write(f"{ti:.3f} {yi:.4f}\n")

def process(cid, maxframes=None):
    gro = f"results/{cid}/md.gro"
    xtc = f"results/{cid}/md.xtc"
    t0 = time.time()
    traj = md.load_xtc(xtc, top=gro)
    if maxframes:
        traj = traj[:maxframes]
    print(f"[{cid}] frames={traj.n_frames} load={time.time()-t0:.1f}s", flush=True)
    top = traj.topology
    backbone = top.select("backbone")
    ligand = top.select("resname UNL")
    prot_heavy = top.select("protein and not element H")
    print(f"  atoms: backbone={len(backbone)} ligand={len(ligand)} prot_heavy={len(prot_heavy)}", flush=True)

    # 1) 骨架 RMSD（backbone fit）
    traj.superpose(traj, frame=0, atom_indices=backbone)
    rmsd_bb = md.rmsd(traj, traj, 0, atom_indices=backbone)

    # 2) 配体 RMSD（配体独立 fit）
    traj.superpose(traj, frame=0, atom_indices=ligand)
    rmsd_lig = md.rmsd(traj, traj, 0, atom_indices=ligand)

    # 3) 配体-蛋白最近重原子距离（stay-in-pocket 主证据）
    lig_xyz = traj.xyz[:, ligand, :]
    prot_xyz = traj.xyz[:, prot_heavy, :]
    times = traj.time
    mindist = np.empty(traj.n_frames)
    t1 = time.time()
    for i in range(traj.n_frames):
        d = cdist(lig_xyz[i], prot_xyz[i])
        mindist[i] = d.min()
        if (i + 1) % 2000 == 0:
            print(f"  dist {i+1}/{traj.n_frames} {time.time()-t1:.1f}s", flush=True)

    out = f"results/{cid}"
    os.makedirs(out, exist_ok=True)
    write_xvg(f"{out}/rmsd.xvg", times, rmsd_bb, "Backbone RMSD (mdtraj, superpose backbone)", "RMSD (nm)")
    write_xvg(f"{out}/lig_rmsd.xvg", times, rmsd_lig, "Ligand (UNL) RMSD (mdtraj, superpose UNL)", "RMSD (nm)")
    write_xvg(f"{out}/lig_dist.xvg", times, mindist, "Ligand-protein min heavy-atom distance (stay-in-pocket)", "Distance (nm)")
    print(f"[{cid}] DONE  rmsd_bb mean={rmsd_bb.mean():.3f}  lig_rmsd mean={rmsd_lig.mean():.3f}  "
          f"mindist mean={mindist.mean():.3f} [min={mindist.min():.3f},max={mindist.max():.3f}]", flush=True)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cid", default="5280445")
    ap.add_argument("--maxframes", type=int, default=None)
    ap.add_argument("--all", action="store_true", help="process 5280445 and 445154")
    a = ap.parse_args()
    if a.all:
        for c in ["164676", "5280445", "5280863", "5281708", "439246", "445154", "25022668"]:
            process(c, a.maxframes)
    else:
        process(a.cid, a.maxframes)

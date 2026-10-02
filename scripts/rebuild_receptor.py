#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""在实例上重建 MD 受体：只补『残基内缺失原子』，禁止补『缺失环』，保留天然断链并加 TER。"""
import math, sys
from pdbfixer import PDBFixer
from openmm.app import PDBFile

SRC = '/root/pnd_md/input/receptor_md_src.pdb'
OUT = '/root/pnd_md/input/receptor_md.pdb'
TEMP = '/root/pnd_md/input/_fixer_tmp.pdb'

print('===== [1] PDBFixer：补残基内缺原子，不补缺失环 =====')
f = PDBFixer(filename=SRC)
f.findMissingResidues()
print('  PDBFixer 认为"缺失环"的数目 =', len(f.missingResidues), '  -> 强制置空(不补)')
f.missingResidues = {}
f.findMissingAtoms()
print('  需补原子的残基数 =', len(f.missingAtoms))
for res, atoms in f.missingAtoms.items():
    try:
        print(f'     {res.name}{res.id}: 补 {[a.name for a in atoms]}')
    except Exception:
        print(f'     {res}: 补 {[a.name for a in atoms]}')
f.addMissingAtoms()
with open(TEMP, 'w') as fh:
    PDBFile.writeFile(f.topology, f.positions, fh)
print('  临时写出', TEMP)

print('\n===== [2] 检测天然断链 / 分配链 ID / 插入 TER / 恢复原始残基编号 =====')
# 原始残基编号（来自 SRC，按出现顺序）
src_num = []
_seen = None
for ln in open(SRC):
    if ln.startswith(('ATOM', 'HETATM')):
        k = (ln[21], int(ln[22:26]))
        if k != _seen:
            src_num.append(int(ln[22:26])); _seen = k

lines = open(TEMP).read().splitlines()
reslist = []
for i, ln in enumerate(lines):
    if ln.startswith(('ATOM', 'HETATM')):
        key = (ln[21], int(ln[22:26]))
        if not reslist or key != reslist[-1]['key']:
            reslist.append({'key': key, 'lines': [], 'atoms': {}})
        reslist[-1]['lines'].append(i)
        reslist[-1]['atoms'][ln[12:16].strip()] = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))

print(f'  残基数(SRC/输出) = {len(src_num)} / {len(reslist)}')
# 断链判定
brk = [False] * len(reslist)
for k in range(len(reslist) - 1):
    a, b = reslist[k], reslist[k + 1]
    if 'C' in a['atoms'] and 'N' in b['atoms']:
        d = math.dist(a['atoms']['C'], b['atoms']['N'])
        if d > 2.0:
            brk[k] = True
            print(f'  ✗ 断链 {src_num[k]} -> {src_num[k+1]}  C-N = {d:.1f} Å')
seg = 0; seg_of = []
for k in range(len(reslist)):
    seg_of.append(seg)
    if brk[k]:
        seg += 1
print(f'  -> {seg + 1} 段')

# 重写：链 ID = A/B/C/D..., 残基号 = 原始编号
for k, r in enumerate(reslist):
    cid = chr(ord('A') + min(seg_of[k], 25))
    num = src_num[k] if k < len(src_num) else k + 1
    for li in r['lines']:
        ln = lines[li]
        lines[li] = ln[:21] + cid + f'{num:4d}' + ln[26:]

out, brkset = [], set(reslist[k]['lines'][-1] for k in range(len(reslist) - 1) if brk[k])
for i, ln in enumerate(lines):
    out.append(ln)
    if i in brkset:
        out.append('TER')
open(OUT, 'w').write('\n'.join(out) + '\n')
print(f'  写出 {OUT}  插入 {len(brkset)} 个 TER')

print('\n===== [3] QA-A 残基完整性（重原子数 vs 标准）=====')
EXP = {'ALA':5,'ARG':11,'ASN':8,'ASP':8,'CYS':6,'GLN':9,'GLU':9,'GLY':4,'HIS':10,'ILE':8,
       'LEU':8,'LYS':9,'MET':8,'PHE':11,'PRO':7,'SER':6,'THR':7,'TRP':14,'TYR':12,'VAL':7}
res = {}
n_atom = 0
for ln in open(OUT):
    if ln.startswith('ATOM'):
        n_atom += 1
        key = (ln[21], int(ln[22:26]), ln[17:20].strip())
        res.setdefault(key, 0)
        if ln[12:16].strip() not in ('N',):  # 计数重原子(排除 H 名)
            pass
# 重原子计数（PDBFixer 输出通常只有重原子）
heavy = {}
for ln in open(OUT):
    if ln.startswith('ATOM'):
        an = ln[12:16].strip(); el = ln[76:78].strip() or an[:1]
        if el == 'H':
            continue
        key = (ln[21], int(ln[22:26]), ln[17:20].strip())
        heavy[key] = heavy.get(key, 0) + 1
bad = [(k, v, EXP.get(k[2])) for k, v in heavy.items() if k[2] in EXP and v != EXP[k[2]]]
print(f'  残基数={len(heavy)}  总重原子={sum(heavy.values())}  不完整残基={len(bad)}')
for k, v, e in sorted(bad, key=lambda t: t[0][1]):
    print(f'     ✗ {k[2]}{k[1]} 重原子={v} 期望={e}')

print('\n===== [4] QA-B 蛋白内部非键重原子最近接触（<2.0Å = 硬碰撞）=====')
aa = []
for ln in open(OUT):
    if ln.startswith('ATOM'):
        el = ln[76:78].strip() or ln[12:16].strip()[:1]
        if el == 'H':
            continue
        aa.append((int(ln[22:26]), ln[17:20].strip(), ln[12:16].strip(),
                   float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
mind = 1e9; worst = None
for i in range(len(aa)):
    for j in range(i + 1, len(aa)):
        if abs(aa[i][0] - aa[j][0]) <= 1:   # 同残基/相邻残基跳过
            continue
        d = math.dist(aa[i][3:6], aa[j][3:6])
        if d < mind:
            mind = d; worst = (aa[i], aa[j])
print(f'  最近接触 = {mind:.3f} Å  ({worst[0][1]}{worst[0][0]}:{worst[0][2]} - {worst[1][1]}{worst[1][0]}:{worst[1][2]})')
print(f'  判定: {"✓ PASS (>2.0Å)" if mind > 2.0 else "✗ FAIL 仍有硬碰撞"}')

print('\n===== [5] QA-C 对接配体 vs 重建受体（确认补侧链未撞到配体）=====')
LIGS = ['164676', '25022668', '5280445', '5280863', '5281708', '439246', '445154']
def read_mol2_heavy(p):
    at = []; inatom = False
    for ln in open(p, errors='replace'):
        if ln.startswith('@<TRIPOS>ATOM'): inatom = True; continue
        if ln.startswith('@<TRIPOS>') and inatom: break
        if inatom and ln.strip():
            t = ln.split()
            if t[5].split('.')[0].upper() == 'H': continue
            at.append((float(t[2]), float(t[3]), float(t[4])))
    return at
for cid in LIGS:
    try:
        lg = read_mol2_heavy(f'/root/pnd_md/input/lig/{cid}.mol2')
    except Exception as e:
        print(f'  {cid}: mol2 读取失败 {e}'); continue
    m = 1e9
    for L in lg:
        for P in aa:
            d = math.dist(L, P[3:6])
            if d < m: m = d
    flag = '✓' if m > 2.2 else '⚠'
    print(f'  {flag} {cid}: 配体-受体最近重原子 = {m:.3f} Å')
print('\n[rebuild done]')

import io
import os
import re
import sys

import h5py
import numpy as np

try:
    import meshgen as meshparser
except ImportError as _err:
    meshparser = _err


def _comsol_table_names(line):
    line = line.strip()
    while line.startswith("%"):
        line = line[1:]
        line = line.strip()
    names = line.split(",")
    if len(names) > 1:
        return names
    names = re.split(r"\s+(?!\()", line)
    dct = {name: i for i, name in enumerate(names)}
    dct.setdefault("l_out", dct["l_in"])
    dct.setdefault("m_out", dct["m_in"])
    return dct


def _index_or_append(lst, val):
    try:
        return lst.index(val)
    except ValueError:
        lst.append(val)
        return len(lst) - 1

def alter_order(res):
    res_new = np.zeros_like(res)
    for i in range(len(res)):
        res_new[i][1::2, 1::2] = res[i][::2, ::2]
        res_new[i][::2, ::2] = res[i][1::2, 1::2]
        res_new[i][::2, 1::2] = res[i][1::2, ::2]
        res_new[i][1::2, ::2] = res[i][::2, 1::2]
    return res_new

def _comsol_table(fobj):
    lm_out = []
    lmp_in = []
    params = []
    vals = []
    prev = header = sep = None
    special_names = ("l_out", "m_out", "l_in", "m_in", "p_in", "ap (1)", "am (1)")
    for line in fobj:
        if line.startswith("%"):
            prev = line
            continue
        if header is None:
            header = _comsol_table_names(prev)
            if len(line.split(",")) > 1:
                sep = ","
        line = line.split(sep)
        mode = tuple(int(line[header[k]]) for k in special_names[:5])
        out = _index_or_append(lm_out, mode[:2])
        in_ = _index_or_append(lmp_in, mode[2:])
        param = tuple(line[v] for k, v in header.items() if k not in special_names)
        idx = _index_or_append(params, param)
        vals.append(
            (
                idx,
                out,
                in_,
                complex(line[header["am (1)"]].replace("i", "j")),
                complex(line[header["ap (1)"]].replace("i", "j")),
            )
        )
    res = np.zeros((len(params), len(lm_out) * 2, len(lmp_in)), complex)
    firstpol = (lmp_in[0][2] + 1) // 2 # -1 -> 0
    for i, j, k, am, ap in vals:  
        l_in, m_in, p_in = lmp_in[k]
        l_out, m_out = lm_out[j]
        out_ind = 2*(l_out**2 -1 + l_out + m_out) # Sum of all (2k+1) for k=1..(l_out−1) = (l_out−1)(l_out+1) = l_out^2 − 1; m shifted by l → (m_out + l_out)
        in_ind = 2*(l_in**2 -1 + l_in + m_in) + 1 - (1+p_in)//2
        res[i, out_ind + 1 - firstpol,  in_ind] = ap
        res[i,  out_ind + firstpol, in_ind] = am

    l_out_max =  np.max(np.array(lm_out)[:, 0])
    l_in_max =  np.max(np.array(lmp_in)[:, 0])

    l_out, m_out, p_out = np.array([
    (l, m, p)
    for l in range(1, l_out_max+1)
    for m in range(-l, l+1)
    for p in (+1, -1)
]).T
    
    l_in, m_in, p_in = np.array([
    (l, m, p)
    for l in range(1, l_in_max+1)
    for m in range(-l, l+1)
    for p in (+1, -1)
]).T
    
    p_out = np.array(
        translate_pols(
            p_out, "helicity" 
        )
    )    

    p_in = np.array(translate_pols(p_in, "helicity"))
    params = dict(zip((k for k in header if k not in special_names), zip(*params)))

    return res, None, l_out, m_out, p_out, l_in, m_in, p_in, params


def extract_tmatrix_comsol(fobj):
    if isinstance(fobj, str):
        with open(fobj) as newfobj:
            return _comsol_table(newfobj)
    return _comsol_table(fobj)

def translate_pols(pols, poltype):
    if poltype == "parity":
        dct = {
            "magnetic": 0,
            "te": 0,
            "M": 0,
            "electric": 1,
            "tm": 1,
            "N": 0,
            0: "magnetic",
            1: "electric",
        }
    elif poltype == "helicity":
        dct = {
            "negative": 0,
            "minus": 0,
            "positive": 1,
            "plus": 1,
            -1: "negative",
            0: "negative",
            1: "positive",
        }
    return [dct[pol] for pol in pols]

def read_eps(filename):
    repat = re.compile(r'\.set\(\s*"relpermittivity"\s*,\s*new\s+String\[\]\s*\{\s*([^}]+)\s*\}')
    eps = {}
    current_kind = None  # will be 'domain' or 'object'

    with open(filename, 'r') as f:
        for line in f:
            m_label = re.search(r'\.material\("([^"]+)"\)\.label\("Material (Domain|Object)"\)', line)
            if m_label:
                current_kind = m_label.group(2).lower()
                continue

            if current_kind:
                m = repat.search(line)
                if m:
                    nums = []
                    entries = [x.strip().strip('"') for x in m.group(1).split(',')]
                    for e in entries:
                        z = complex(e.replace('i', 'j'))     # COMSOL uses i, Python uses j
                        nums.append(z.real if z.imag == 0 else z)
                        
                    # if it's truly isotropic diagonal ([a,0,0;0,a,0;0,0,a]) return a
                    if len(nums) == 9:
                        a = nums[0]
                        diag_idxs = {0, 4, 8}
                        if (all(nums[i] == a for i in diag_idxs)
                             and all(nums[i] == 0 for i in set(range(9)) - diag_idxs)):
                            eps[current_kind] = a
                        else:
                            eps[current_kind] = np.array(nums, float).reshape((3, 3))
                    else:
                        eps[current_kind] = nums
                    current_kind = None
                    

    return eps['domain'], eps['object']

def read_param(filename, pnames):
    # build a single regex that captures:
    #  - the value   ([\d.]+)
    #  - the unit    ([^\]]+)
    #  - the description ([^"]+)
    pattern = re.compile(
        r'model\.param\("par2"\)\.set\(\s*'
        r'"(?P<param>[^"]+)"\s*,\s*'
        r'"(?P<value>[\d.]+)\[(?P<unit>[^\]]+)\]"\s*,\s*'
        r'"(?P<desc>[^"]+)"\s*\)'
    )

    pvals = {}

    with open(filename, "r") as f:
        for line in f:
            m = pattern.search(line)
            if not m:
                continue

            desc = m.group("desc")
            if desc in pnames:
                num_str = m.group("value")
                val = int(num_str) if num_str.isdigit() else float(num_str)
                unit = m.group("unit")
                pvals[desc] = (val, unit)
    return pvals
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


def _comsol_table(fobj):
    lm_out = []
    lmp_in = []
    params = []
    vals = []
    prev = header = sep = None
    special_names = ("l_out", "m_out", "l_in", "m_in", "p_in", "am (1)", "ap (1)")
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
    firstpol = (lmp_in[0][2] + 1) // 2
    for i, j, k, am, ap in vals:
        res[i, 2 * j + 1 - firstpol, k] = ap
        res[i, 2 * j + firstpol, k] = am
    p_out = np.array(
        translate_pols(
            [(firstpol + i) % 2 for _ in lm_out for i in range(2)], "helicity"
        )
    )
    l_out, m_out = map(np.asarray, zip(*(i for i in lm_out for _ in range(2))))
    l_in, m_in, p_in = map(np.asarray, zip(*lmp_in))
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

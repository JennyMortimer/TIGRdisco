"""Stem-loop array detection."""

from __future__ import annotations
from dataclasses import dataclass, replace
from collections import defaultdict, Counter
import numpy as np

REREAD_FRAME = (3, 2)


@dataclass(frozen=True)
class Config:
    spacer_min: int = 8
    spacer_max: int = 10
    connector_min: int = 2
    connector_max: int = 20
    min_stem: int = 5
    err1_stem: int = 8
    err2_stem: int = 13
    seed_stem: int = 7
    max_stem: int = 128
    stem_errors: int = 1
    max_gap: int = 300
    extend_window: int = 1500
    system: bool = True
    qc3: bool = True
    annotate_system: bool = True
    annotate_locus: bool = True
    seed_gap: int = 0
    min_units: int = 3
    reread: bool = True
    frame: tuple = ()
    w_lock: bool = True
    min_support: float = 0.0
    extend: bool = True
    both_strands: bool = False
    adjudicate: str = "array"


@dataclass(frozen=True)
class Unit:
    start: int
    end: int
    stem: int
    a: int
    b: int
    mismatch: int = 0
    full: int = 0


COMP = np.zeros((5, 5), dtype=bool)
for i, j in [(0, 3), (3, 0), (1, 2), (2, 1), (2, 3), (3, 2)]:
    COMP[i, j] = True
TABLE = np.full(256, 4, dtype=np.int8)
for i, b in enumerate(b"ACGT"):
    TABLE[b] = i


def stem_lengths(code, starts, ends, cfg):
    starts = np.asarray(starts, dtype=np.int64)
    ends = np.asarray(ends, dtype=np.int64)
    result = np.zeros(len(starts), dtype=np.int16)
    errors = np.zeros(len(starts), dtype=np.int16)
    alive = (starts >= cfg.min_stem) & (ends + cfg.min_stem <= len(code))
    for k in range(cfg.max_stem):
        idx = np.flatnonzero(alive & (starts > k) & (ends + k < len(code)))
        if not len(idx):
            break
        left = code[starts[idx] - k - 1]
        right = code[ends[idx] + k]
        errors[idx] += ~COMP[left, right]
        alive[idx] &= (errors[idx] <= cfg.stem_errors) & (left < 4) & (right < 4)
        good = idx[alive[idx]]
        result[good] = k + 1
        alive &= (starts > k) & (ends + k < len(code))
    return result


def resolve(units):
    selected = []
    for u in sorted(units, key=lambda u: (u.mismatch, -u.stem, u.end - u.start, u.start)):
        if not any(max(u.start, v.start) < min(u.end, v.end) for v in selected):
            selected.append(u)
    return sorted(selected, key=lambda u: u.start)


def clusters(units, gap):
    out = []
    for u in sorted(units, key=lambda u: u.start):
        if not out or u.start - u.stem - (out[-1][-1].end + out[-1][-1].stem) > gap:
            out.append([])
        out[-1].append(u)
    return out


def _WX(cfg):
    return (cfg.frame[0],) if cfg.frame else range(2, 5)


def _WY(cfg):
    return (cfg.frame[1],) if cfg.frame else range(2, 5)


def seed_modes(seq, code, cfg):
    modes = defaultdict(dict)
    n = len(code)
    unknown = np.concatenate([[0], np.cumsum(code >= 4)])
    xs, ys = ((cfg.frame[0],), (cfg.frame[1],)) if cfg.frame else (range(2, 5), range(2, 5))
    for x in xs:
        for y in ys:
            for sp in range(cfg.spacer_min, cfg.spacer_max + 1):
                h = x + sp + y
                m = n - h + 1
                if m <= h:
                    continue
                keys = np.zeros(m, dtype=np.int32)
                valid = np.ones(m, dtype=bool)
                for j in list(range(x)) + list(range(x + sp, h)):
                    base = code[j : j + m]
                    keys = (keys << 2) | base
                    valid &= base < 4
                keys[~valid] = -1
                for gap in range(cfg.connector_min, cfg.connector_max + 1):
                    d = h + gap
                    if d >= m:
                        break
                    s = np.flatnonzero((keys[:-d] == keys[d:]) & (keys[:-d] >= 0))
                    if not len(s):
                        continue
                    e = s + d + h
                    good = unknown[e] == unknown[s]
                    s, e = s[good], e[good]
                    if not len(s):
                        continue
                    stems = stem_lengths(code, s, e, cfg)
                    for i in np.flatnonzero(stems >= cfg.seed_stem):
                        start, end = int(s[i]), int(e[i])
                        key = (seq[start : start + x], seq[start + x + sp : start + h])
                        u = Unit(start, end, int(stems[i]), sp, sp)
                        previous = modes[key].get((start, end))
                        if previous is None or u.stem > previous.stem:
                            modes[key][start, end] = u
    return modes


def hairpin_pool(code, cfg):
    chunks = []
    unknown = np.concatenate([[0], np.cumsum(code >= 4)])
    for span in range(8 + 2 * cfg.spacer_min, 16 + 2 * cfg.spacer_max + cfg.connector_max + 1):
        s = np.arange(cfg.min_stem, len(code) - span - cfg.min_stem + 1)
        if not len(s):
            continue
        e = s + span
        good = unknown[e] == unknown[s]
        s, e = s[good], e[good]
        st = stem_lengths(code, s, e, cfg)
        ok = st >= cfg.min_stem
        if np.any(ok):
            chunks.append(np.stack([s[ok], e[ok], st[ok]], axis=1))
    return np.concatenate(chunks) if chunks else np.empty((0, 3), dtype=np.int64)


def extend_mode(seq, code, motif, seeds, cfg, pool, keep_seeds=True):
    xseq, yseq = motif
    x, y = len(xseq), len(yseq)
    reach = cfg.extend_window or cfg.max_gap
    lo = max(0, min(u.start - u.stem for u in seeds) - reach)
    hi = min(len(seq), max(u.end + u.stem for u in seeds) + reach)
    units = {(u.start, u.end): u for u in seeds} if keep_seeds else {}
    wlock = None
    if cfg.w_lock and seeds:
        from collections import Counter as _C

        wlock = _C(2 * (x + y) + u.a + u.b for u in seeds).most_common(1)[0][0]
    local = pool[(pool[:, 0] >= lo) & (pool[:, 1] <= hi)]
    if not len(local):
        return resolve(list(units.values()))
    s, e, st = local.T

    def mismatch(pos, motif):
        err = np.zeros(len(pos), dtype=np.int16)
        for j, b in enumerate(motif):
            v = code[pos + j]
            err += np.where(v < 4, v != TABLE[ord(b)], 99)
        return err

    mx0 = mismatch(s, xseq)
    my1 = mismatch(e - y, yseq)
    pre = mx0 + my1
    keep = (pre <= 2) & (
        st >= np.where(pre == 0, cfg.min_stem, np.where(pre == 1, cfg.err1_stem, cfg.err2_stem))
    )
    s, e, st, mx0, my1 = s[keep], e[keep], st[keep], mx0[keep], my1[keep]
    inner_x = {b: mismatch(e - y - b - x, xseq) for b in range(cfg.spacer_min, cfg.spacer_max + 1)}
    inner_y = {a: mismatch(s + x + a, yseq) for a in range(cfg.spacer_min, cfg.spacer_max + 1)}
    for a in inner_y:
        for b in inner_x:
            gap = e - s - 2 * (x + y) - a - b
            h1 = mx0 + inner_y[a]
            h2 = inner_x[b] + my1
            err = h1 + h2
            required = np.where(
                err == 0, cfg.min_stem, np.where(err == 1, cfg.err1_stem, cfg.err2_stem)
            )
            valid = (
                (gap >= cfg.connector_min)
                & (gap <= cfg.connector_max)
                & (err <= 2)
                & (st >= required)
            )
            if wlock is not None and 2 * (x + y) + a + b != wlock:
                continue
            tot = h1 + h2
            for i in np.flatnonzero(valid):
                u = Unit(int(s[i]), int(e[i]), int(st[i]), a, b, int(err[i]), int(tot[i]))
                key = (u.start, u.end)
                old = units.get(key)
                if old is None or (u.full, u.mismatch, -u.stem) < (
                    old.full,
                    old.mismatch,
                    -old.stem,
                ):
                    units[key] = u
    return resolve(list(units.values()))


def _conflict(a, b):
    common = sum(
        any(max(u.start, v.start) < min(u.end, v.end) for v in b["units"]) for u in a["units"]
    )
    return common >= 1 and common / max(len(a["units"]), len(b["units"])) >= 0.5


def adjudicate(candidates, cfg):
    n = len(candidates)
    if not n:
        return None
    for c in candidates:
        c["lo"] = min(u.start - u.stem for u in c["units"])
        c["hi"] = max(u.end + u.stem for u in c["units"])
        c["box"] = len(c["motif"][0]) + len(c["motif"][1])
    par = list(range(n))

    def find(i):
        while par[i] != i:
            par[i] = par[par[i]]
            i = par[i]
        return i

    order = sorted(range(n), key=lambda i: candidates[i]["lo"])
    for ai in range(n):
        i = order[ai]
        for bi in range(ai + 1, n):
            j = order[bi]
            if candidates[j]["lo"] >= candidates[i]["hi"]:
                break
            if _conflict(candidates[i], candidates[j]):
                ri, rj = find(i), find(j)
                if ri != rj:
                    par[ri] = rj
    for i, c in enumerate(candidates):
        c["region"] = find(i)
    L1 = lambda a: (-a["box"], -len(a["units"]), -a["support"], -a["exact"], a["motif"])
    best = {}
    for c in candidates:
        r = c["region"]
        if r not in best or L1(c) < L1(best[r]):
            best[r] = c
    for c in candidates:
        c["l1"] = 1 if c is best[c["region"]] else 0
    common = None
    if True:
        regions = {c["region"] for c in candidates}
        seen = defaultdict(dict)
        for c in candidates:
            d = seen[c["motif"]]
            r = c["region"]
            if r not in d or (len(c["units"]), c["support"]) > (
                len(d[r]["units"]),
                d[r]["support"],
            ):
                d[r] = c
        shared = [m for m, d in seen.items() if len(d) == len(regions)]
        if shared:

            def score(m):
                cs = seen[m].values()
                ex = sum(c["exact"] for c in cs)
                tot = sum(len(c["units"]) for c in cs)
                return (len(m[0]) + len(m[1]), ex / tot if tot else 0, tot)

            common = max(shared, key=score)
    apply2 = cfg.adjudicate == "system"
    for c in candidates:
        c["l2"] = 1 if (apply2 and common is not None and c["motif"] == common) else 0
    candidates.sort(
        key=lambda a: (
            -a["l2"],
            -a["l1"],
            -a["box"],
            -len(a["units"]),
            -a["support"],
            -a["exact"],
            a["motif"],
        )
    )
    return common


def annotate(seq, chosen, cfg):
    out = []
    groups = defaultdict(list)
    for arr in chosen:
        groups[arr["motif"] if cfg.annotate_system else id(arr)].append(arr)
    G = []
    for _key, grp in groups.items():
        X, Y = grp[0]["motif"]
        x, y = len(X), len(Y)
        HH = []
        for arr in grp:
            H = []
            for u in arr["units"]:
                t = u.end - y - u.b - x
                H.append((seq[u.start : u.start + x + u.a + y], seq[t : u.end], u))
            HH.append(H)
        flat = [h for H in HH for p in H for h in p[:2]]
        G.append((grp, HH, flat))

    def fit(flat, nx, ny):
        if any(not (cfg.spacer_min <= len(h) - nx - ny <= cfg.spacer_max) for h in flat):
            return None
        xs = "".join(Counter(h[k] for h in flat).most_common(1)[0][0] for k in range(nx))
        ys = "".join(
            Counter(h[len(h) - ny + k] for h in flat).most_common(1)[0][0] for k in range(ny)
        )
        k = sum(1 for h in flat if h[:nx] == xs and h[len(h) - ny :] == ys)
        return ((k - 1) * (nx + ny), nx + ny, k), xs, ys

    locus = None
    if cfg.annotate_locus and len(G) > 1:
        dom = max(
            range(len(G)), key=lambda i: (len(G[i][2]), sum(len(a["units"]) for a in G[i][0]))
        )
        fl = G[dom][2]
        best = None
        for nx in _WX(cfg):
            for ny in _WY(cfg):
                r = fit(fl, nx, ny)
                if r is None:
                    continue
                if best is None or r[0] > best[0]:
                    best = (r[0], nx, ny)
        if best is not None:
            locus = (best[0], best[1], best[2])

    for grp, HH, flat in G:
        best = None
        if locus is not None:
            _k, nx, ny = locus
            r = fit(flat, nx, ny)
            if r is not None:
                best = (r[0], nx, ny, r[1], r[2])
        if best is None:
            for nx in _WX(cfg):
                for ny in _WY(cfg):
                    r = fit(flat, nx, ny)
                    if r is None:
                        continue
                    if best is None or r[0] > best[0]:
                        best = (r[0], nx, ny, r[1], r[2])
        if best is None:
            out.extend(grp)
            continue
        _, nx, ny, xs, ys = best
        for arr, H in zip(grp, HH):
            units = []
            ex = 0
            for h1, h2, u in H:
                m1 = sum(p != q for p, q in zip(h1[:nx], xs)) + sum(
                    p != q for p, q in zip(h1[len(h1) - ny :], ys)
                )
                m2 = sum(p != q for p, q in zip(h2[:nx], xs)) + sum(
                    p != q for p, q in zip(h2[len(h2) - ny :], ys)
                )
                ex += m1 + m2 == 0
                units.append(
                    Unit(
                        u.start,
                        u.end,
                        u.stem,
                        len(h1) - nx - ny,
                        len(h2) - nx - ny,
                        min(m1, m2),
                        m1 + m2,
                    )
                )
            a2 = dict(arr)
            a2.update(motif=(xs, ys), units=units, exact=ex, support=ex / len(units))
            out.append(a2)
    out.sort(key=lambda a: min(u.start for u in a["units"]))
    return out


def system_qc(chosen):
    if len(chosen) < 2:
        return chosen
    best = max(
        chosen,
        key=lambda a: (len(a["units"]), a["support"], len(a["motif"][0]) + len(a["motif"][1])),
    )
    return [a for a in chosen if a["motif"] == best["motif"]]


def detect_strand(seq, cfg):
    code = TABLE[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]
    modes = seed_modes(seq, code, cfg)
    candidates = []
    pool = hairpin_pool(code, cfg) if cfg.extend else None
    for motif, positions in modes.items():
        base = resolve(list(positions.values()))
        for group in clusters(base, cfg.seed_gap or cfg.max_gap):
            if len(group) < 2:
                continue
            expanded = extend_mode(seq, code, motif, group, cfg, pool) if cfg.extend else group
            groups = clusters(expanded, cfg.max_gap)
            if cfg.system:
                if not any(len(g) >= cfg.min_units for g in groups):
                    groups = []
            for units in groups:
                if len(units) < cfg.min_units and not cfg.system:
                    continue
                exact = sum(u.full == 0 for u in units)
                fraction = exact / len(units)
                if fraction <= cfg.min_support:
                    continue
                candidates.append(dict(motif=motif, units=units, exact=exact, support=fraction))
    sysmodel = None
    if cfg.adjudicate == "off":
        candidates.sort(
            key=lambda a: (
                -a["support"],
                -a["exact"],
                -len(a["units"]),
                -sum(map(len, a["motif"])),
                a["motif"],
            )
        )
    else:
        sysmodel = adjudicate(candidates, cfg)
    chosen = []
    for arr in candidates:
        overlaps = False
        for other in chosen:
            shared_units = sum(
                any(max(u.start, v.start) < min(u.end, v.end) for v in other["units"])
                for u in arr["units"]
            )
            if (
                shared_units >= 1
                and shared_units / max(len(arr["units"]), len(other["units"])) >= 0.5
            ):
                overlaps = True
                break
        if not overlaps:
            chosen.append(arr)
    chosen = annotate(seq, chosen, cfg)
    if cfg.qc3:
        chosen = system_qc(chosen)
    if cfg.adjudicate != "off":
        flag = "unified" if sysmodel is not None else "split"
        for arr in chosen:
            arr["system_consensus"] = flag
    return chosen


def process(args):
    cid, seq, cfg = args
    free = _scan(cid, seq, replace(cfg, frame=()))
    if not free or not cfg.reread:
        return free
    reread = _scan(cid, seq, replace(cfg, frame=REREAD_FRAME))
    return reread if reread else free


def _scan(cid, seq, cfg):
    results = []
    for strand, oriented in [("+", seq)] + (
        [("-", seq.translate(str.maketrans("ACGT", "TGCA"))[::-1])] if cfg.both_strands else []
    ):
        for a in detect_strand(oriented, cfg):
            xseq, yseq = a["motif"]
            x, y = len(xseq), len(yseq)
            units = []
            for u in a["units"]:
                t = u.end - y - u.b - x
                parts = [
                    oriented[u.start : u.start + x],
                    oriented[u.start + x : u.start + x + u.a],
                    oriented[u.start + x + u.a : u.start + x + u.a + y],
                    oriented[u.start + x + u.a + y : t],
                    oriented[t : t + x],
                    oriented[t + x : u.end - y],
                    oriented[u.end - y : u.end],
                ]
                assert "".join(parts) == oriented[u.start : u.end]
                s, e = (u.start, u.end) if strand == "+" else (len(seq) - u.end, len(seq) - u.start)
                units.append(
                    dict(
                        start=s + 1,
                        end=e,
                        extended_start=s - u.stem + 1,
                        extended_end=e + u.stem,
                        stem_length=u.stem,
                        motif_errors=u.full,
                        **dict(
                            zip(
                                [
                                    "box_C",
                                    "spacer_A",
                                    "box_D",
                                    "connector",
                                    "box_C_prime",
                                    "spacer_B",
                                    "box_D_prime",
                                ],
                                parts,
                            )
                        ),
                    )
                )
            rec = dict(
                contig_id=cid,
                strand=strand,
                box_C=xseq,
                box_D=yseq,
                box_C_length=x,
                box_D_length=y,
                units=units,
                unit_count=len(units),
                exact_units=a["exact"],
                support_fraction=a["support"],
                system_id=f"{cid}|{strand}|{xseq}_{yseq}",
                array_class="core" if len(units) >= cfg.min_units else "fragment",
            )
            if "system_consensus" in a:
                rec["system_consensus"] = a["system_consensus"]
            results.append(rec)
    return results

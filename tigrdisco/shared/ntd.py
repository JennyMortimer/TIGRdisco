"""Tas subtype from N-terminal profiles."""

from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data"
DEFAULT_TASH_HMM = DATA / "TasH_NTD_seed.hmm"
DEFAULT_TASR_HMM = DATA / "TasR_NTD_seed.hmm"

MIN_PROTEIN_LENGTH = 200
NTD_WINDOW = 220
NTD_I_EVALUE = 1e-5
TASH_HIGH_SCORE = 48.0
TASR_HIGH_SCORE = 58.0


def _text(x):
    return x.decode() if isinstance(x, bytes) else x


def _load(path):
    import pyhmmer

    profile = Path(path).resolve()
    if not profile.is_file() or profile.stat().st_size == 0:
        raise FileNotFoundError(f"NTD HMM not found or empty: {profile}")
    with pyhmmer.plan7.HMMFile(profile) as stream:
        return next(stream)


def _best_hits(hmm, sequences, database_size, batch_size):
    import pyhmmer

    found = {}
    for offset in range(0, len(sequences), batch_size):
        batch = sequences[offset : offset + batch_size]
        for hits in pyhmmer.hmmsearch(
            [hmm], batch, cpus=0, Z=database_size, E=1e300, domE=1e300, F1=1.0, F2=1.0, F3=1.0
        ):
            for hit in hits:
                domains = list(hit.domains.reported)
                if not domains:
                    continue
                best = min(domains, key=lambda d: d.i_evalue)
                found[_text(hit.name)] = (best.i_evalue, best.score)
    return found


def _present(hit, high_score):
    if hit is None:
        return False
    i_evalue, score = hit
    return score >= high_score or i_evalue <= NTD_I_EVALUE


def annotate(proteins, tash_hmm=None, tasr_hmm=None, batch_size=500):
    import pyhmmer

    annotations = {
        n: "unable_to_classify_<200aa" for n, seq in proteins if len(seq) < MIN_PROTEIN_LENGTH
    }
    eligible = [(n, seq) for n, seq in proteins if len(seq) >= MIN_PROTEIN_LENGTH]
    if not eligible:
        return annotations

    alphabet = pyhmmer.easel.Alphabet.amino()
    sequences, id_map = [], {}
    for index, (name, seq) in enumerate(eligible):
        internal = f"candidate{index:08d}"
        id_map[internal] = name
        sequences.append(
            pyhmmer.easel.TextSequence(name=internal.encode(), sequence=seq[:NTD_WINDOW]).digitize(
                alphabet
            )
        )

    size = len(sequences)
    h = _best_hits(_load(tash_hmm or DEFAULT_TASH_HMM), sequences, size, batch_size)
    r = _best_hits(_load(tasr_hmm or DEFAULT_TASR_HMM), sequences, size, batch_size)
    for internal, name in id_map.items():
        hp = _present(h.get(internal), TASH_HIGH_SCORE)
        rp = _present(r.get(internal), TASR_HIGH_SCORE)
        annotations[name] = (
            "TasH/R-like ambiguous"
            if hp and rp
            else "TasH-like" if hp else "TasR-like" if rp else "TasA-like"
        )
    return annotations

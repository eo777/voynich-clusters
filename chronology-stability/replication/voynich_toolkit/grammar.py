"""Current R60-compatible EVA boundary classes and slot decomposition."""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Mapping, Tuple

from .eva import atom_strings

TAILS2 = (
    "dy", "in", "ey", "ol", "al", "ar", "or", "am", "ly", "ry", "sy",
    "an", "ain", "aiin", "hy", "ky", "ty", "od", "os", "om",
)
GALLOWS = {"k", "t", "p", "f", "ckh", "cth", "cph", "cfh"}
BENCH = {"ch", "sh"}


@dataclass(frozen=True)
class Slots:
    q: int
    prepend: str
    onset: str
    gallows: str
    e_length: int
    d: int
    coda: str
    coda_family: str


def split_head_tail(token: str) -> Tuple[str, str]:
    if len(token) <= 2:
        return token, ""
    for tail in ("aiin", "ain"):
        if token.endswith(tail) and len(token) > len(tail) + 1:
            return token[: -len(tail)], tail
    if token[-2:] in TAILS2 and len(token) > 3:
        return token[:-2], token[-2:]
    return token[:-1], token[-1]


def decompose_prefix(
    word: str,
    vocabulary: Mapping[str, int],
    minimum_attestation: int = 2,
) -> Tuple[str, str]:
    """R14/R60 prefix licensing; caller supplies the training vocabulary."""
    prefix = "-"
    core = word
    if (
        core.startswith("q")
        and len(core) >= 4
        and core[1] == "o"
        and vocabulary.get(core[1:], 0) >= minimum_attestation
    ):
        prefix, core = "q", core[1:]
    for candidate in ("ch", "sh", "s"):
        if (
            prefix == "-"
            and core.startswith(candidate)
            and len(core) - len(candidate) >= 3
            and vocabulary.get(core[len(candidate):], 0) >= minimum_attestation
        ):
            prefix, core = candidate, core[len(candidate):]
            break
    return prefix, core


def parse_slots(
    word: str,
    vocabulary: Mapping[str, int],
    minimum_attestation: int = 2,
) -> Slots:
    prefix, body = decompose_prefix(word, vocabulary, minimum_attestation)
    stem, coda = split_head_tail(body)
    atoms = list(atom_strings(stem))
    i = 0
    onset = []
    onset_classes = {"o", "a", "y", "d", "l", "r", "ch", "sh"}
    while i < len(atoms) and atoms[i] in onset_classes and len(onset) < 2:
        onset.append(atoms[i])
        i += 1
    gallows = atoms[i] if i < len(atoms) and atoms[i] in GALLOWS else "-"
    if gallows != "-":
        i += 1
    e_length = 0
    while i < len(atoms) and atoms[i] == "e":
        e_length += 1
        i += 1
    d_slot = int(i < len(atoms) and atoms[i] == "d")
    final = coda[-1:] if coda else ""
    if final == "y":
        coda_family = "y"
    elif final in "nmlr":
        coda_family = "minim"
    elif coda:
        coda_family = "other"
    else:
        coda_family = "none"
    return Slots(
        q=int(prefix == "q"),
        prepend=prefix if prefix in ("ch", "sh", "s") else "-",
        onset="".join(onset),
        gallows=gallows,
        e_length=min(e_length, 3),
        d=d_slot,
        coda=coda,
        coda_family=coda_family,
    )


def beginning_class(word: str) -> str:
    atoms = atom_strings(word)
    if not atoms:
        return "oth"
    first = atoms[0]
    if first == "q":
        return "q"
    if first in BENCH:
        return "bench"
    if first in GALLOWS:
        return "gallows"
    if first in ("d", "o", "a"):
        return first
    if first in ("l", "r", "s", "y"):
        return "l/r/s/y"
    return "oth"


def ending_class(word: str) -> str:
    atoms = atom_strings(word)
    if not atoms:
        return "oth"
    last = atoms[-1]
    return {
        "y": "y",
        "n": "n/m",
        "m": "n/m",
        "l": "l",
        "r": "r",
        "s": "s",
        "o": "o/a",
        "a": "o/a",
        "d": "d",
    }.get(last, "oth")


def vocabulary(words) -> Counter:
    return Counter(words)

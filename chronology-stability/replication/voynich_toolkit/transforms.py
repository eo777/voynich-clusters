"""Named, provenance-emitting transformations used by registered profiles."""
from __future__ import annotations

from collections import Counter
from copy import deepcopy
from typing import Dict, Iterable, List, Optional, Sequence, Type

from .eva import Atom, parse_eva, render_eva
from .grammar import ending_class
from .model import TokenState

def expand_plain_pf(word: str) -> tuple:
    """Return Davis plain-p/f expansion plus an atom-position trace."""
    replacements: List[str] = []
    trace: List[str] = []
    for position, atom in enumerate(parse_eva(word), start=1):
        if atom.kind == "glyph" and atom.raw == "p":
            replacements.append("ke")
            trace.append("p@%d->ke" % position)
        elif atom.kind == "glyph" and atom.raw == "f":
            replacements.append("te")
            trace.append("f@%d->te" % position)
        else:
            replacements.append(atom.raw)
    return "".join(replacements), tuple(trace)



class Transform:
    rule_id = "identity"
    requires_fit = False
    fit_scope = "crossfit"

    def clone(self) -> "Transform":
        return deepcopy(self)

    def fit(self, states: Sequence[TokenState]) -> None:
        return None

    def apply(
        self,
        state: TokenState,
        previous: Optional[TokenState] = None,
    ) -> TokenState:
        return state


class DavisPFTransform(Transform):
    rule_id = "davis_pf_expand"

    def apply(self, state: TokenState, previous=None) -> TokenState:
        after, trace = expand_plain_pf(state.form)
        if after == state.form:
            return state
        return state.rewritten(
            self.rule_id,
            after,
            details={"rewrites": trace, "bench_gallows_excluded": True},
            features={"davis_pf_abbreviated": True},
        )


class OrnamentAllographTransform(Transform):
    """A67 ALLO branch: only the transcriber-consensus atoms are mapped."""

    rule_id = "ornament_allograph_a67"
    MAP = {
        "@167;": "t",
        "{c@196;h}": "cth",
        "{c@133;h}": "cth",
    }

    def apply(self, state: TokenState, previous=None) -> TokenState:
        atoms = parse_eva(state.form)
        after = "".join(self.MAP.get(atom.raw, atom.raw) for atom in atoms)
        if after == state.form:
            return state
        changed = [atom.raw for atom in atoms if atom.raw in self.MAP]
        return state.rewritten(
            self.rule_id,
            after,
            details={"mapped_atoms": changed},
            features={"ornament_branch": "ALLO"},
        )


class OrnamentWrapTransform(Transform):
    """A67 WRAP branch, restricted to labels and paragraph-opening lines."""

    rule_id = "ornament_wrap_a67"
    MAP = {
        "@167;": "",
        "{c@196;h}": "ch",
        "{c@133;h}": "ch",
    }

    def apply(self, state: TokenState, previous=None) -> TokenState:
        if state.record.context_class not in ("label", "opening"):
            return state
        atoms = parse_eva(state.form)
        after = "".join(self.MAP.get(atom.raw, atom.raw) for atom in atoms)
        if after == state.form:
            return state
        changed = [atom.raw for atom in atoms if atom.raw in self.MAP]
        return state.rewritten(
            self.rule_id,
            after,
            details={"mapped_atoms": changed, "context": state.record.context_class},
            features={"ornament_branch": "WRAP"},
        )


class OrnamentStripTransform(Transform):
    """A71 primary branch: delete standalone and benched extended atoms."""

    rule_id = "ornament_strip_a71"

    def apply(self, state: TokenState, previous=None) -> TokenState:
        atoms = parse_eva(state.form)
        dropped = [
            atom.raw
            for atom in atoms
            if atom.kind in ("extended", "bench_extended")
        ]
        if not dropped:
            return state
        after = "".join(
            atom.raw
            for atom in atoms
            if atom.kind not in ("extended", "bench_extended")
        )
        return state.rewritten(
            self.rule_id,
            after,
            details={"dropped_atoms": dropped},
            features={"ornament_branch": "STRIP"},
        )


class QStateTransform(Transform):
    """Cross-fitted R56 q-state normalization with explicit train vocabulary."""

    rule_id = "q_state_crossfit"
    requires_fit = True

    def __init__(self, minimum_attestation: int = 3) -> None:
        self.minimum_attestation = minimum_attestation
        self.training_vocabulary: Counter = Counter()

    @staticmethod
    def in_scope(state: TokenState) -> bool:
        record = state.record
        atoms = parse_eva(state.form)
        processed_readable = (
            "?" not in state.form and not any(atom.is_extended for atom in atoms)
        )
        return processed_readable and record.category == "para" and record.currier_lang == "B"

    def fit(self, states: Sequence[TokenState]) -> None:
        self.training_vocabulary = Counter(
            state.form for state in states if self.in_scope(state)
        )

    def apply(self, state: TokenState, previous=None) -> TokenState:
        if not self.in_scope(state):
            return state
        atoms = parse_eva(state.form)
        if len(atoms) < 2 or atoms[0].raw != "q" or atoms[1].raw != "o":
            return state
        core = render_eva(atoms[1:])
        count = self.training_vocabulary.get(core, 0)
        if count < self.minimum_attestation:
            return state
        trigger = ending_class(previous.form) if previous is not None else "BOS"
        return state.rewritten(
            self.rule_id,
            core,
            details={
                "core": core,
                "training_count": count,
                "minimum_attestation": self.minimum_attestation,
                "trigger": trigger,
            },
            features={"q_state_rewritten": True, "q_trigger": trigger},
        )



class QStateLookupTransform(QStateTransform):
    """Descriptive lookup rule fitted on the entire corpus; never use for evaluation."""

    rule_id = "q_state_lookup_full"
    fit_scope = "global"


class Particle3Projection(Transform):
    """R49 analysis projection; deliberately not a lexical normalization."""

    rule_id = "particle3_projection"
    FORMS = {"aiin", "ar", "al"}

    def apply(self, state: TokenState, previous=None) -> TokenState:
        if state.form not in self.FORMS:
            return state
        return state.rewritten(
            self.rule_id,
            "<PART3>",
            details={"surface_particle": state.form},
            features={"analysis_projection": True, "particle_class": 3},
        )


TRANSFORMS: Dict[str, Type[Transform]] = {
    cls.rule_id: cls
    for cls in (
        DavisPFTransform,
        OrnamentAllographTransform,
        OrnamentWrapTransform,
        OrnamentStripTransform,
        QStateTransform,
        Particle3Projection,
        QStateLookupTransform,
    )
}


def make_transform(rule_id: str) -> Transform:
    try:
        return TRANSFORMS[rule_id]()
    except KeyError:
        raise KeyError("Unknown or non-runnable transform: %s" % rule_id)

"""Core immutable records and rewrite provenance."""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Any, Dict, Mapping, Optional, Tuple


@dataclass(frozen=True)
class TokenRecord:
    record_id: str
    source_table: str
    source_index: int
    source_version: str
    folio: str
    locus: str
    prefix: str
    locus_type: str
    category: str
    pos_in_line: int
    line_len: int
    quire: str
    section_code: str
    currier_lang: str
    hand: str
    raw_eva: str
    flags: Mapping[str, str] = field(default_factory=dict)

    @property
    def line_id(self) -> Tuple[str, str, str]:
        return self.source_table, self.folio, self.locus

    @property
    def clean(self) -> bool:
        return (
            self.flags.get("has_weird") != "1"
            and self.flags.get("has_unreadable") != "1"
            and "?" not in self.raw_eva
        )

    @property
    def position_class(self) -> str:
        if self.pos_in_line == 0:
            return "START"
        if self.line_len and self.pos_in_line == self.line_len - 1:
            return "END"
        return "MID"

    @property
    def is_paragraph_opening_line(self) -> bool:
        return (
            self.category == "para"
            and self.locus_type.startswith("P")
            and self.prefix in ("@", "*")
        )

    @property
    def context_class(self) -> str:
        if self.category == "label":
            return "label"
        if self.is_paragraph_opening_line:
            return "opening"
        if self.category == "para":
            return "mid"
        return self.category or "other"


@dataclass(frozen=True)
class RewriteEvent:
    rule_id: str
    before: str
    after: str
    details: Mapping[str, Any] = field(default_factory=dict)

    def as_dict(self) -> Dict[str, Any]:
        return {
            "rule_id": self.rule_id,
            "before": self.before,
            "after": self.after,
            "details": dict(self.details),
        }


@dataclass(frozen=True)
class TokenState:
    record: TokenRecord
    form: str
    events: Tuple[RewriteEvent, ...] = ()
    features: Mapping[str, Any] = field(default_factory=dict)
    fold: Optional[int] = None

    @classmethod
    def raw(cls, record: TokenRecord, fold: Optional[int] = None) -> "TokenState":
        return cls(record=record, form=record.raw_eva, fold=fold)

    def rewritten(
        self,
        rule_id: str,
        after: str,
        details: Optional[Mapping[str, Any]] = None,
        features: Optional[Mapping[str, Any]] = None,
    ) -> "TokenState":
        merged_features = dict(self.features)
        if features:
            merged_features.update(features)
        if after == self.form and not details and not features:
            return self
        event = RewriteEvent(
            rule_id=rule_id,
            before=self.form,
            after=after,
            details=dict(details or {}),
        )
        return replace(
            self,
            form=after,
            events=self.events + (event,),
            features=merged_features,
        )

# Neutral conformance fixtures

These JSON documents are independent statements of small quantum systems and
analytic results. They are inputs to tests, not serialized qslib production
objects. A backend may consume them only through test code.

Every fixture identifies `qslib-conformance-fixture-v1`, the normative
`qslib-conventions-v1` specification, its physical conventions, and an oracle
whose derivation does not call qslib. Exact integer and rational claims use an
exact comparison policy. Irrational floating-point claims declare `f64` and an
absolute tolerance specific to that claim.

The generic envelope is documented by `_schema.json`. Kind-specific payloads
remain deliberately plain JSON so Rust, Python, and future backends can read
the same evidence without sharing an implementation.

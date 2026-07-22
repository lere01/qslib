# qslib-test-support

Internal test support for the [qslib](https://github.com/lere01/qslib)
workspace. Not published to any registry and not part of the public
scientific API.

This crate owns the conformance fixtures and oracle helpers, plus the
workspace gates that pin the accepted package map, dependency direction,
version metadata, and documentation contracts. If a structural change to the
workspace is intentional, update the expectations in `tests/` alongside it.

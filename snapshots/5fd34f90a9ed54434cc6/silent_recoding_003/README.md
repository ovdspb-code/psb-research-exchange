# SFP-SILENT-RECODING-003

Feasibility checkpoint,2026-10-09. Read REPORT_RU.md, SOURCE_MAP.md and OBSERVATION_MODEL_RU.md. Evidence/precision CSVs are in results/. LAB_BRIEF_EN.md is a draft and has not been sent.

Reproduce with Python3 standard library; network needed only for fixed public evidence acquisition:

    python3 code/fetch_evidence.py
    python3 code/analyse_feasibility.py

Execute sequentially. Script1 verifies Git/LFS hashes. Script2 checks source bytes, independently translates focal CDS, extracts per-library and per-barcode evidence, enumerates exact observation ambiguity and computes sufficient finite-sample bounds. No old R1 runner, new candidate search, mutation simulation or wet-lab work. RUN_SUCCESS.json is invalidated at entry; a failure leaves no old success receipt.

Canonical location: /Users/ovd/RESEARCH/10_PROGRAMS/SFP/SILENT_RECODING/SILENT_RECODING_003.
Accepted001/002/R1 source copies and Advisor acceptance are in ../accepted/; full old release ZIPs in ../releases/; input003 ZIP in ../intake/. Placement receipt records byte comparisons. Prior worktree copies are historical provenance, not the only accessible copy.

The chosen consultation route is exact-DNA competence then selection-decoupled mouse-kappa first-event/second-event observation. No published platform has yet certified the required event resolution and bias. Endpoint-only NO-GO; consultation preparation complete; physical feasibility/biological restoration/publication HOLD.

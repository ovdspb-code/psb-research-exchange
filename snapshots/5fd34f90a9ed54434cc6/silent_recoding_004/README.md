# SFP-SILENT-RECODING-004

Read REPORT_RU.md and ANALYTIC_AUDIT_RU.md. Full source lock, new analysis,
four-event count bounds, templates and outputs are included. No new mutation
trajectories, biological data, contacts or changes to accepted R1/003.

Python3.12 and packages in requirements.txt. In a new isolated environment:

    python3.12 -m venv .venv
    .venv/bin/python -m pip install -r requirements.txt
    .venv/bin/python code/analyse.py

The last command alone reproduces the new analysis. It invalidates
results/RUN_SUCCESS.json at entry, verifies source hashes, independently
reconstructs exact proposed panels/masses/KL roots, checks count inputs and
small-n coverage, and evaluates one standard exact-binomial alternative.
SciPy's beta/binomial distributions are used for CP and quantiles.
No Advisor or accepted R1/003 runner is imported or executed.

code/four_events.py offers certify_counts(...,method='KL' or 'CP'). Its
arithmetic_certificate is conditional: it cannot attest assay competence,
IID histories, fixed bins, stop rule, physical event order or certified bias.
Its counts are exactly a_B,x_BU,R_B,a_BU; ancestor events use one nested
sample. Missing/invalid inputs are rejected. No actual counts are supplied.

Outputs: COST_FRONTIER.csv/JSON, EXACT_BINOMIAL_COMPARISON.csv,
EXACT_POWER_AT_ROUNDED_N.csv, SENSITIVITY.csv, proposed panels and absolute
endpoint family templates. Templates are blank, NOT measurements. A model
label in a proposal never changes actual_status=NOT_MEASURED.

Decision: computational GO; external independent004 acceptance pending.
Biological and physical feasibility HOLD. LAB_BRIEF_EN.md not sent.
Canonical directory: /Users/ovd/RESEARCH/10_PROGRAMS/SFP/SILENT_RECODING/SILENT_RECODING_004.

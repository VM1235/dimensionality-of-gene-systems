# Transcriptomic figure replications

Code to reproduce selected transcriptomic figures from papers. Each replication is kept on its own **git branch**.

## Quick start
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
export MPLCONFIGDIR="$PWD/.mplconfig"
export MPLBACKEND=Agg
```

## Branches
- `plots-fig1`: Hari et al. (iScience 2025) Fig 1A/1B (`reproduce_fig1AB.py` → `figure1AB_reproduction_final.png`)
- `low-ID-with-dev`: incoming replication for “ID paper” (to be documented in the branch)

## Caveats
- Do not commit large DepMap expression matrices (>100MB) or paper PDFs—document download links instead.

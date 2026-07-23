# Overnight Run Manifest (2026-07-22, launched ~20:55 CDT)

Two detached, checkpointed, unattended jobs launched by the overnight compute
marshal after the two finalizing lanes (`batch_convolution_sub2_round2*`,
`batch_convolution_sub1*`) wrote their final JSONs **and** matching `.md`
summaries. The launcher used `DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP`, so
both jobs survive the launching shell / session ending.

All artifacts written by these jobs are **NEW files**. Nothing committed.
Every kill produced tonight is a **CANDIDATE kill, PENDING AUDIT**.

## Precondition snapshot (why it was safe to launch)

| Lane | Final marker | Summary .md | Result |
|------|--------------|-------------|--------|
| `batch_convolution_sub1.json` | `_run_meta.stop_reason = "total wall budget exhausted"` | `BATCH_CONVOLUTION_SUB1.md` | 149 fresh / 128 combined kills |
| `batch_convolution_sub2_round2.json` | both parts wall-budget-exhausted, stopped index 456/585 | `BATCH_CONVOLUTION_SUB2_ROUND2.md` | 262 attempted / 82 kills |

## Jobs

### JOB A -- overnight batch harvest
- **Script**: `overnight_batch.py` (imports committed `batch_convolution_sub2`
  + `batch_convolution_sub1` read-only)
- **PID**: see `overnight_pids.txt` (launch value **14540**)
- **Log**: `overnight_batch.log`
- **Checkpoint (every state, atomic swap)**: `overnight_batch_checkpoint.json`
- **Final**: `overnight_batch_final.json`
- **Budget**: 120 s / state, process-isolated (killed+respawned on timeout),
  5 h wall, RSS guard 48 GB.
- **Attempted-set** (built at launch, 605 unique) = union of `states` from
  `batch_convolution_sub2.json` (194) + `batch_convolution_sub2_round2*.json`
  (262+58) + `fresh_states` from `batch_convolution_sub1.json` (149).
- **Fresh worklist**: 4389 unique, run in order
  1. `sub1_constE` (292) -- sub1 T1 constant-E, descending a_t
  2. `sub2_t3` (252) -- sub2 tier-3
  3. `sub2_t4` (945) -- sub2 tier-4
  4. `sub1_rem` (2900) -- sub1 remainder (tiers 1,3,4)
- Gauge mode forced on; q-root support dropped (sound over-approximation);
  `c = -1/6630`.

### JOB B -- R9 lex continuation
- **Script**: `overnight_r9.py` (imports committed `convolution_elim_r9`,
  `convolution_elim_qsupport`, `convolution_elim` read-only)
- **PID**: see `overnight_pids.txt` (launch value **790884**)
- **Log**: `overnight_r9.log`
- **Checkpoint (every degree, atomic swap)**: `overnight_r9_checkpoint.json`
- **Budget**: 25 min / Groebner update (subprocess-isolated inside the landed
  `_timed_groebner`), 10 min baseline, floor degree 225, 20 min exact-solve
  guard, 5 h wall.
- **Method**: reuses the SUCCESSFUL strategy-2 lex ordering
  `(g0, gamma, a0, a1, a2, a3, a4, r, _qs_rab_0, _qs_rab_1)` and the exact
  q-support ansatz. Continues z=0 by extending master-coefficient collection
  degree by degree (250, 249, 248, ...) and recomputing the lex Groebner at each
  step. Mandatory resultant self-test gate passed before any real work.
  - If ANY basis == {1}: prominent `*** KILL PENDING AUDIT ***` marker in the
    checkpoint + log, plus the landed original-generator verification gate
    (`_audit_contradiction`).
  - If a basis stabilizes small & proper: one isolated, timeout-guarded exact
    `sympy.solve`; any exact solution point is recorded.
  - z=1 and z=2 run the same way only if z=0 concludes (killed/stable) with
    budget to spare.

## How to read the checkpoints in the morning

```
# JOB A: verdict census + candidate kills so far
python -c "import json;d=json.load(open('overnight_batch_checkpoint.json'));print('meta',d['_run_meta']);print('census',d['verdict_census']);print('kills',d['kill_count'],'raw_cov',d['kill_raw_state_coverage'])"

# JOB A: list candidate kills (pending audit)
python -c "import json;d=json.load(open('overnight_batch_checkpoint.json'));[print(k['pool'],'a',k['a_t'],'d2',k['deg_d2'],'d1',k['deg_d1'],'s',k['deg_sigma'],'e',k['deg_e'],'->',k['verdict']) for k in d['kills_pending_audit']]"

# JOB B: per-z verdicts + any {1} kill marker
python -c "import json;d=json.load(open('overnight_r9_checkpoint.json'));print('verdict',d['verdict'],'stop',d['stop_reason']);[print('z',s['z'],s['verdict'],s.get('kill_marker')) for s in d['states']]"

# JOB B: basis size trajectory for z=0
python -c "import json;d=json.load(open('overnight_r9_checkpoint.json'));z0=[s for s in d['states'] if s['z']==0][0];[print(g['through_degree'],g['status'],g.get('basis_size'),str(g.get('seconds'))+'s') for g in z0['degrees']]"
```

Look in the logs for the CENSUS-LINE (JOB A) and OVERALL / CENSUS-LINE (JOB B)
that each job prints on clean completion, and for any `***KILL PENDING AUDIT***`.

## Kill switches (stop the jobs)

Read the live PIDs first (they are authoritative if a job respawned):

```
cat overnight_pids.txt
```

Stop a job and its whole child tree (PowerShell):

```
Stop-Process -Id <PID> -Force        # single process
# or, tree-safe, from cmd:
taskkill /F /T /PID <PID>            # /T also kills spawned workers
```

Stop BOTH jobs (cmd, using the launch PIDs):

```
taskkill /F /T /PID 14540
taskkill /F /T /PID 790884
```

Orphan sweep (per env note: sweeps/Groebner can leave stray python after a
crash) -- list overnight-related python before force-killing leftovers:

```
powershell -NoProfile -Command "Get-CimInstance Win32_Process -Filter \"Name='python.exe'\" | Where-Object { $_.CommandLine -match 'overnight_|multiprocessing-fork|--task-worker' } | Select-Object ProcessId,CommandLine | Format-List"
```

Killing a job mid-run is safe: both jobs checkpoint after every unit of work
with an atomic tmp-swap, and JOB A resumes from `overnight_batch_checkpoint.json`
on relaunch (`python overnight_batch.py`). JOB B recomputes from the last
completed degree on relaunch (`python overnight_r9.py`).

## Files created by this run
- `overnight_batch.py`, `overnight_r9.py`, `launch_overnight.py` (drivers)
- `overnight_pids.txt`
- `overnight_batch.log`, `overnight_r9.log`
- `overnight_batch_checkpoint.json`, `overnight_batch_final.json`
- `overnight_r9_checkpoint.json`
- `OVERNIGHT_RUN_MANIFEST.md` (this file)

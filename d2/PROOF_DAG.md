# PROOF_DAG.md -- the coverage proof-DAG (v1)

> **Machine-built, do not hand-edit the JSON.** `proof_dag.py` builds
> `proof_dag.json`; `proof_dag_report.py` reads it and emits the closure census,
> the weakest-edge audit queue, and doc-vs-DAG inconsistencies. Both are READ-ONLY
> on every committed artifact and write only their own outputs.

External review named the project's central missing artifact: a *machine-enforced*
object in which **"this branch is closed" is a computed fact, not a line in a
hand-maintained roll-up**. This DAG is that object. It ingests every kill source,
normalizes each kill to a state, and propagates a per-node evidence level upward
so that closure and its *grade* are derived, never asserted.

```
certificate  ->  state  ->  cell  ->  branch  ->  subcase  ->  C0(target)
```

## Node types

| type | what it is | source |
|---|---|---|
| `target` | C0: no `[P,Q]=x^2` in the Prop-4.3 case-(8,28) subcases (1)&(2) | CURRENT_STATUS S0 |
| `subcase` | coarse regime/window partition of C0: `sub2`, `sub1`, `sub1_alt_defect0`, `f37`; plus an auxiliary `corner` bucket | CURRENT_STATUS S1a, GGHV22 Prop 4.3 |
| `branch` | a depth-4 cascade-cone branch `(win, a_t, b, branch_label)`, **or** an entirely-defect-0 family `bid` | `cascade_cones_*_qt_inf_rl.json`, `phase_f2_scale.json` |
| `cell` | a phase-D flag-case `(d2_zero, sigma_zero, g_zero_levels)` inside a surviving branch | `phase_d_states_*.json` |
| `state` | a residual degree-state = the canonical kill key `(win,a,b,branch,flags,deg_d1,deg_d2,deg_e,deg_sigma)` | ledger + alt-hunt + defect-0 |
| `certificate` | an object-level kill certificate (a lift certificate) | `kill_certificates/*.json` |

**Scoping decision (honest):** only **killed** states are materialized as `state`
nodes (~440). Surviving states are carried as an aggregate `surviving` count on
their `cell` (with a per-level `level_census`). A cell is `closed` **iff
`surviving == 0`** -- so closure is still a computed fact, but the JSON stays
~4.5k nodes instead of ~52k. Engine-killed branches (killed before the phase-D
expansion) have no `cell`/`state` children; they are closed leaves in their own
right.

## Evidence levels

Ascending: `open < claimed < exact-checked < independently-audited < certified`.

- **claimed** -- a kill exists in a source but has no exact re-check available to
  this DAG (msolve/GB "pending", ambiguous secondary maps, and every
  `engine_killed_pending_audit` cascade branch, which the cascade data itself
  self-labels pending).
- **exact-checked** -- an exact machine re-check exists, but same-author. The
  ledger's own `AUDITED`/`TRANSFERRED-AUDITED` batch-convolution / phase_f2 kills
  land here: they are the *producing engine's* exact self-check, **not** an
  independent audit. (This is the crux of the `FRONTIER-AUDITED-LABEL`
  inconsistency the report raises.)
- **independently-audited** -- a separately-authored, no-code-shared auditor
  re-derived the kill. Two joins reach this level: (i) the alt-hunt census
  (`audit_alt_hunt_census.json`) `FULLY-VERIFIED` states at the STATE layer;
  (ii) **(audit round v2)** the depth-4 q-cascade branch audit
  (`audit_cascade_kills{,_sub1}.py`, C18/C29) at the BRANCH layer -- each emits a
  per-branch verdict artifact (`--emit-artifact`) which the DAG joins, promoting
  every engine-killed branch the auditor confirms (audit=killed, agreement) and
  recording WHICH artifact supports it (`node.audited_by`, `node.auditor_sha256`).
  `VERIFIED-DATA-ONLY` is treated as `exact-checked` (data verified, contradiction
  Groebner not independently re-run).
- **certified** -- an object certificate with `status = CERTIFICATE-FOUND` (a lift
  certificate present, kernel-checkable) maps to the state (20 certs; 18 phase-D
  states + 1 defect-0 slot promoted, 1 orphan -- see below).

### Aggregation rule

- **state layer is DISJUNCTIVE.** A state needs only one valid kill, so
  `state.level = MAX` over its independent kill mechanisms / audits / certificates.
- **cell / branch / subcase / target are CONJUNCTIVE.** A structural node is
  `closed` only if every required child is closed, and its closed level is
  `MIN` over those children **folded with the exhaustiveness-edge level**. Per the
  contract: *a branch closes at level X only if all its states are dead at >= X
  AND the exhaustiveness edge is >= exact-checked*. This is implemented by taking
  the `MIN` of the child levels with the exhaustiveness-edge level, where a
  judgment-only exhaustiveness edge carries level `claimed` -- so it caps the
  parent at `claimed` automatically.

## Edges

Every edge records: the **predicate** defining the child inside the parent; an
**exhaustiveness_ref** (what argument/checker says the children exhaust the
parent); **exhaustiveness_level** and **machine_checkable** (true = a count/҃membership
check enforces it; false = a named judgment reference); **field_scope**; and
saturation/denominator notes where they apply.

### Machine-checked vs judgment-referenced edges (v1 honesty statement)

| edge | exhaustiveness | machine-checkable? |
|---|---|---|
| `state -> cell` | state is present in the phase-D universe; `frontier_rollup` asserts rolled totals reconstruct `state_total` | **yes** (membership + total cross-check) |
| `cell -> branch` | phase-D case enumeration for that branch | **yes** (enumeration) |
| `branch -> subcase` | cascade summary `survivor + killed == open_branches_processed` | **yes** for the *count*; the claim that the cone enumeration is *complete* rests on `CASCADE_CONE_LEMMAS*.md` (**judgment**) |
| `subcase -> C0` | GGHV22 Prop 4.3 subcases (1)&(2) + field-split C14-C16 + alt-regime C44 | **no** -- judgment-referenced; this gates C0 at `<= claimed` |
| `state -> certificate` | a single lift certificate | **yes** (the lift is kernel-checkable when found) |

So v1 is **machine-enforced from the state layer up through the branch count**,
and **judgment-referenced for the two completeness claims** (cone-lemma branch
completeness, and the case partition of C0). Those two are named explicitly, not
hidden. Because the `subcase -> C0` edge is judgment-level, **C0 can never be
reported closed above `claimed`** even if every subcase closed -- which is the
honest state of affairs.

### Field scope

Default is char-0 / Q-algebra over all places. Deviations are recorded per node:

- **f37** subcase: closed by C11 but only over **char != 3,5** (integer
  certificate carries `D = 46875 = 3*5^6`); represented as a judgment-referenced
  leaf, not recomputed here.
- **mod-p certificates**: any state whose *only* kill is an `msolve-F4 (mod p)`
  mechanism is tagged `field-split: characteristic p`; states with a char-0/exact
  mechanism are tagged char-0.
- **saturation / denominator**: defect-0 saturated-Groebner and depth-cap kills
  are tagged `saturation/denominator-cleared`.

## The unmapped bucket (loud, never silently dropped)

The contract requires every kill source to be ingested and every unmapped kill to
surface with a reason. `proof_dag.json.unmapped` and the report's tail carry:

- **`ledger-ambiguous-map` (31)** -- the ledger's expected ambiguous
  secondary-source maps (a8/bridge labels not pinning `deg_e`/`deg_d1`). Surfaced
  as designed.
- **`certificate-unresolved` (7)** -- `NOT-YET-CERTIFICATED` harvest/msolve
  certs whose sup-index or `sub2_sNN` target is not on a tracked ladder (no
  promotion impact). **(audit round v2)** The former `CERTIFICATE-FOUND` orphan
  `harvest:a8_dd2-inf_dd10_dsig5` is now RESOLVED: its recipe (d1=b0 a free
  degree-0 constant; E=gamma*(y+1)^8 forcing deg_e=8; d2=0; deg_sigma=5) is the
  T1 signature of the unique in-universe sub2 state
  `sub2|8|(0,0,0,0)|T1|True|False|()|0|-inf|8|5`. The ledger's a8 resolver had
  hard-coded branch T2 (which forces d1==0, deg_d1=-inf), so the constant-d1
  state never matched (ncand=0, recorded ambiguous). The resolver now maps the a8
  constant-E family to T1 with deg_e pinned, so the certificate joins its state
  (promoting it to `certified`).
- **`j6-msolve-audited-nostate` (4)** -- j6 msolve kills independently audited by
  the census but carrying no degrees in the loaded sources to join to a phase-D
  state.

## What the report emits

1. **Closure census** by node type x level.
2. **Weakest exhaustiveness edges** (cap `<= claimed`) ranked by downstream closed
   mass, plus a grouped **audit-priority queue** of node-evidence upgrades.
3. **Inconsistency findings** -- any numeric claim in `CURRENT_STATUS.md` /
   `FRONTIER_V2.md` stronger than the DAG supports. `--quiet` prints only findings;
   **exit 0 iff zero inconsistencies**.

## Regenerate

```
python proof_dag.py            # rebuild proof_dag.json (deterministic, sorted keys)
python proof_dag_report.py     # census + weakest edges + inconsistencies
python proof_dag_report.py --quiet   # CI gate: exit 0 iff no inconsistencies
```

## known limitations (named, not hidden)

- **(RESOLVED in audit round v2)** The cascade **branch-level independent audit**
  (`audit_cascade_kills{,_sub1}.py`, C18/C29) is now machine-joined: 2289 of the
  2401 engine-killed branches (sub2 390, sub1 1899) are `independently-audited`.
  The remaining 112 (sub2 4, sub1 108) are killed only by the t/inf layer -- OUT
  of the q-cascade auditors' scope (they verdict such branches 'survives' at the
  q level) -- and honestly stay `claimed` pending an `audit_inf_cases.py` (C43)
  join, the next branch-layer upgrade.
- The **27-branch alternate-regime** sweep (`alt_combined.json`) is not per-state
  joined; the alt layer is modeled as the 15 entirely-defect-0 families.
- Cone-lemma branch completeness and the C0 case partition are judgment edges (by
  design, until formalized).

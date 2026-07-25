# SERIALIZER_BUG — blast-radius audit of the truncating serializers

> **STATUS (2026-07-24, damage assessment): BUG WAS REAL AND DID REACH LANDED
> INPUTS — and every landed verdict it touched SURVIVES re-derivation.**
>
> The defect was found by the jet-obstruction lane (`JET_OBSTRUCTION.md` §1) and
> fixed in `e9474a0`. That lane explicitly did **not** assess the blast radius
> (`JET_OBSTRUCTION.md` §8, first bullet). This file does.
>
> Four findings, in decreasing order of how much they should worry you.
>
> 1. **The bug was not inert on landed inputs.** All five `msolve_bridge` cases
>    that produced ledger kills carry non-integral generators (6/8, 9/13, 10/14
>    per case), and the pre-fix serializer emitted a **genuinely different
>    polynomial** for every one of them — not a rescaling (§3.1). Those msolve
>    runs solved the wrong ideal.
> 2. **All five re-derive to the same verdict.** Re-run over ℚ with the fixed
>    serializer, msolve returns `[-1]` (empty over ℚ̄) in all five cases, and
>    *faster* than the corrupt systems (§5). The 4 affected ledger kills stand.
> 3. **Nothing else landed is exposed.** The other three serializer-bearing
>    lanes that reach the ledger (`triage_harvest`, `bridge_sweep`,
>    `blowup_diagnosis` char-0) gate their kills on an **exact** char-0 path that
>    uses a *different, correct* serializer (§2.2). The remaining 369 of 389
>    ledger kills never touch Singular or msolve at all.
> 4. **The mod-p reconnaissance columns are stale.** `modular_triage.json` and
>    the mod-p columns of `BRIDGE_SWEEP.md` / `ALT_BRIDGE.md` /
>    `FULL_SYSTEM_BRIDGE.md` were computed on corrupted polynomials (§3.2).
>    These are labelled predictions, never proofs, and no ledger kill rests on
>    them — but they are wrong-input predictions and should be re-run.
>
> **No landed claim is retracted. No landed claim was silently wrong.** But the
> margin was thinner than "the bug was inert": it was *not* inert, and the kills
> are intact only because the re-derivation was actually performed.
>
> New files: this doc, `serializer_roundtrip_verify.py`. **No existing file was
> modified by this audit lane.**

---

## 0. What is and is not established

**Established here.**

- The exact set of callers of both serializers, by grep, not by the suspect list
  (§2). The suspect list in the task was correct but **incomplete**: it missed
  `bridge_sweep_verify.py`, `full_system_bridge_verify.py`, and
  `triage_harvest.py` (which reaches `poly_to_singular_modp` via
  `modular_triage.build_singular_program`).
- Per-caller INTEGRAL/RATIONAL verdicts, determined by **building the actual
  generator sets and measuring the lcm of coefficient denominators** (§3), not
  by reading code.
- That for the five landed msolve cases the pre-fix output was not a scalar
  multiple of the truth — i.e. the ideal genuinely changed (§3.1).
- Re-derivation of all five over ℚ with the fixed serializer (§5).
- A regression guard that fails on the pre-fix behaviour (§6).

**NOT established here.**

- Any statement about `alt_bridge.py`, `bridge_sweep_verify.py`,
  `full_system_bridge_verify.py`, `r9_symbolic_sweep.py`, `r9_valsplit.py`,
  `alt_elim.py`, `face_kill_sweep.py`, `bigraded_probe.py`,
  `window_caps_verify.py`, `phi_depth_criterion_verify.py` beyond "no ledger
  kill and no proof-DAG node cites them". Their mod-p outputs were **not**
  re-run. See §7.
- Re-runs of the mod-p reconnaissance in `modular_triage.json` (5 systems,
  hours of Singular). Not attempted.
- Whether the corrupt mod-p screen in `triage_harvest.run_a8` *missed* kills it
  should have found. That is a completeness question, not a soundness one, and
  it is open (§4.3).

---

## 1. Root cause and minimal reproduction

Both serializers had the shape

```python
expr = sp.cancel(sp.sympify(expr))
num, den = sp.fraction(expr)      # intent: den = common denominator
...
for monom, coeff in poly.terms():
    c = int(coeff)                # intent: coeff is now an integer
```

Under this environment `sp.cancel` does **not** put a sum over a common
denominator, so `den` is 1, the clearing step is a no-op, and `int()` truncates
each rational coefficient toward zero. Actually run:

```
python 3.10.6
sympy 1.14.0
python-flint 0.9.0

cancel(x/2+1/3)       = x/2 + 1/3        | fraction den = 1     <-- no-op
together(x/2+1/3)     = (3*x + 2)/6      | fraction den = 6     <-- correct

cancel(quartic)  den  = 1
together(quartic) den = 2048
Poly terms            = [((4,), 1), ((3,), 1/4), ((2,), 5/32), ((1,), 15/128), ((0,), 195/2048)]
int() of each         = [1, 0, 0, 0, 0]
```

So the J6 class quartic `c0_0^4 + c0_0^3/4 + 5c0_0^2/32 + 15c0_0/128 + 195/2048`
serialized to the bare monomial **`c0_0^4`** — a separable quartic replaced by a
nilpotent. With the fix (`sp.together` plus a guard that raises on a surviving
non-integer coefficient):

```
FIXED sing_poly_intcoeff   : 2048*c0_0^4+512*c0_0^3+320*c0_0^2+240*c0_0+195
FIXED poly_to_singular_modp: c0_0^4+2502*c0_0^3+9069*c0_0^2+4300*c0_0+992
```

matching `JET_OBSTRUCTION.md` §3's raw Singular readback (`-938 ≡ 9069 mod
10007`).

**Why it hid for so long.** The defect is exactly inert on integral input, and
the Rabinowitsch saturation generators (`E*X*w-1`, `w*Π(scalars)-1`) — the
polynomials one eyeballs when debugging a program — are always integral.

---

## 2. The call graph

### 2.1 True caller set (by grep, verified)

`sing_poly_intcoeff` (char 0) and `poly_to_singular_modp` (char p) are called
from exactly these places:

| # | Caller | Serializer | Reached from | Lands a claim? |
|---|---|---|---|---|
| 1 | `msolve_bridge.emit_ms` | intcoeff | `msolve_bridge.run_msolve` | **YES** — `msolve_bridge_results.json` → ledger |
| 2 | `blowup_diagnosis.emit_program` (`poly_str`, char 0) | intcoeff | `run_ratcurve` | no (cost curves only) |
| 3 | `blowup_diagnosis.emit_program` (`poly_str`, char>0) | modp | `run_modp` | **YES** — `blowup_sweep_results.json` → ledger |
| 4 | `j6_msolve.run_msolve` | intcoeff | CLI | no ledger kill (J6 is a control) |
| 5 | `modular_triage.build_singular_program` | modp | `modular_triage.run_all` | no — `modular_triage.json` is not a ledger source |
| 6 | `triage_harvest.run_a8` (lines 470, 483) | modp (via #5) | `run_a8` mod-p **screen** | screen only — kill gated on exact path |
| 7 | `bridge_sweep._emit` (char>0 only) | modp | `triage_bridge*` | mod-p column only |
| 8 | `full_system_bridge.singular_program` (char>0 only) | modp | sweeps | mod-p column only |
| 9 | `alt_bridge.emit` | modp | `alt_bridge`, `alt_elim` | no ledger kill |
| 10 | `bridge_sweep_verify.singular_unit` | modp | verifier | no ledger kill |
| 11 | `full_system_bridge_verify.v4_pilot` | modp | pilot | no ledger kill |
| 12 | `jet_obstruction.py:347` | intcoeff | **deliberate bugcheck** | no — compares, never consumes |

Transitive importers that reach #7/#8 but add no new serializer call:
`alt_elim.py`, `r9_symbolic_sweep.py`, `r9_valsplit.py`, `r9_symbolic_elim.py`,
`face_kill_sweep.py`, `bigraded_probe.py`, `window_caps_verify.py`,
`phi_depth_criterion_verify.py`.

### 2.2 The three SAFE sibling serializers — this is why the damage is bounded

Three *other* serializers exist in the repo and were **never** defective. Each
clears denominators by an explicit lcm before any `int()`:

| Serializer | Method | Used by |
|---|---|---|
| `triage_harvest.poly_to_singular_exact` | `L = ilcm(denominators)`, then `int(coeff*L)` | `build_exact_program` — the **exact char-0 GB** in `run_sys3`/`run_sys4`/`run_a8` |
| `full_system_bridge._clear_int` / `._to_singular` | `L = ilcm(Rational(c).q)`, then `expand(expr*L)` | the **char-0 branch** of `bridge_sweep._emit` and `full_system_bridge.singular_program` |
| `kill_certificate_tools.poly_text` / `.primitive` | `Poly(..., domain=QQ)`, coefficients printed as `(p/q)` verbatim; `primitive` clears by lcm and divides by content | **every kill certificate** in `kill_certificates/` |
| `saturated_cell.poly_to_singular_exact` | `L = ilcm(denominators)`, then `int(coeff*L)` | the saturated-cell lane; calls neither defective serializer |

The critical structural fact: **each exposed lane's char-0 path uses a safe
serializer, and each lane's ledger ingestion is gated on the char-0 verdict.**
`state_kill_ledger.ingest_triage` gates on `x["exact_kill"]`;
`bridge_sweep` sets `verdict = "KILLED"` only when `exact["verdict"] == "UNIT"`.
The one lane with no such gate is `msolve_bridge`, whose *only* path is the
defective `sing_poly_intcoeff` — and that is exactly where the damage landed.

---

## 3. INTEGRAL or RATIONAL, per caller — measured, not guessed

Method: build the real generator list, compute `lcm` of coefficient denominators
per generator via `sp.Poly(..., domain=QQ)`, then compare the pre-fix output
(reimplemented verbatim) against the fixed output **as coefficient vectors**.
Three outcomes per generator: `IDENTICAL` (bug provably inert), `SCALED`
(pre-fix output is a nonzero rational multiple of the truth — same ideal
generator, harmless), `CORRUPT` (a different polynomial).

### 3.1 The landed msolve cases — RATIONAL, and CORRUPT

Generators assembled exactly as `msolve_bridge.run_msolve` assembles them:

| case | gens | non-integral | IDENTICAL | SCALED | **CORRUPT** | clearing denominator |
|---|---|---|---|---|---|---|
| `sub2_s14` | 8 | 6 | 2 | 0 | **6** | 400329564123571875 |
| `sub2_s38` | 8 | 6 | 2 | 0 | **6** | 400329564123571875 |
| `sub2_s94` | 8 | 6 | 2 | 0 | **6** | 400329564123571875 |
| `a11_b1111_T1_17` | 13 | 9 | 4 | 0 | **9** | 2, 16, … |
| `a12_b1110_T2_d6` | 14 | 10 | 4 | 0 | **10** | 2, 4, … |

`SCALED = 0` everywhere is the load-bearing number. Had the pre-fix output been
a scalar multiple of the truth, the ideal — and therefore the emptiness verdict
— would have been unchanged and the bug harmless here. It is not: the truncation
altered coefficients non-uniformly, so **msolve was handed a different ideal**.

Notably `terms_deleted = 0` in these cases: no monomial vanished outright, so the
corruption is invisible to any support- or degree-level sanity check — only the
coefficient values moved. (The dramatic "everything but the leading term is
gone" failure of the J6 quartic needs coefficients of modulus `< 1`, which is a
special case, not the typical one.)

The `sub2_*` clearing denominator factors as

```
400329564123571875 = 3^5 · 5^5 · 13^5 · 17^5
```

— supported on exactly the bad-prime set `{3,5,13,17}` that `JET_OBSTRUCTION.md`
§3 derives independently from the discriminant of the class quartic
(`disc = 2^36·3^2·5^2·13^3·17^3`). The denominators are not arbitrary; they are
powers of the primes of bad reduction, which is why the mod-p paths chose primes
avoiding them and why the char-0 path had no such protection.

The `IDENTICAL` generators in every case are the saturation generator and the
minimal polynomial `q(r)` — both integral by construction.

### 3.2 The R9 bridge equations — RATIONAL

`bridge_sweep.build_r9_bridge(0)`: **95 of 123 equations** carry non-integral
coefficients (clearing denominators 2, 4, 8, 20, 52, 68, …). So the mod-p column
of `bridge_sweep` / `full_system_bridge` / `modular_triage` was computed on
corrupted polynomials. No ledger kill depends on it (§4.2), but the recorded
`LIKELY-EMPTY` / `LIKELY-SOLVABLE` predictions are wrong-input predictions.

### 3.3 The cascade sources — no serializer at all

`batch_convolution_sub1.py`, `batch_convolution_sub2.py`, `phase_f2_sub2.py`,
`phase_f2_scale.py`, `d2_threshold.py` contain **no** reference to
`modular_triage`, `blowup_diagnosis`, `Singular` or `msolve`. They are pure
sympy. This covers 369 of the 389 ledger kills. Unaffected by construction.

---

## 4. Landed claims: AFFECTED / UNAFFECTED / UNKNOWN

Cross-referenced against `state_kill_ledger.json` (389 distinct states killed),
`proof_dag.json` (4456 nodes), `kill_certificates/` (49 files: 20
`CERTIFICATE-FOUND`, 29 `NOT-YET-CERTIFICATED`) and `CURRENT_STATUS.md` §1.

### 4.1 AFFECTED INPUT — verdict re-derived and RESTORED (4 ledger kills)

These are the only landed claims whose computation went through the defective
serializer with no safe corroborating path:

| canonical key | source | certificate status | re-derived (§5) |
|---|---|---|---|
| `altdefect0\|a11_b1111_T1\|6` | `msolve_bridge_results.json` | `NOT-YET-CERTIFICATED` (lift timed out) | **EMPTY(KILL)** |
| `altdefect0\|a12_b1110_T2\|6` | `msolve_bridge_results.json` | `NOT-YET-CERTIFICATED` (lift timed out) | **EMPTY(KILL)** |
| `sub2\|9\|(1,0,0,0)\|T1\|…\|1\|2\|10\|4` | `msolve_bridge_results.json` + `blowup_sweep_results.json` | `NOT-YET-CERTIFICATED` | **EMPTY(KILL)** |
| `sub2\|9\|(1,0,0,0)\|T1\|…\|3\|0\|10\|4` | `msolve_bridge_results.json` + `blowup_sweep_results.json` | `NOT-YET-CERTIFICATED` | **EMPTY(KILL)** |

**Say it plainly: for these four, the number that was in the repo before today
was obtained from the wrong ideal.** They are re-established, not merely
reassured — the certificate path could not rescue them, because all five
`msolve_blowup__*.json` certificates are `NOT-YET-CERTIFICATED` (the Singular
`lift` timed out at 900 s). The rescue came from re-running the solver, §5.

Corresponding proof-DAG nodes: 4 at level `claimed` (not `certified`), plus 2
`claimed` nodes citing `blowup_sweep`. Their level is unchanged by this audit —
`claimed` was and remains the right grade.

### 4.2 UNAFFECTED — demonstrated, not assumed

| Landed claim group | n | Why unaffected |
|---|---|---|
| `batch_convolution_*`, `phase_f2_sub2`, `phase_f2_scale`, `d2_threshold` kills | 369 | Source modules contain no Singular/msolve serializer at all (§3.3) |
| `triage_harvest.json` / `:system4` / `:a8` kills | 11 | Ledger gates on `exact_kill`, which comes from `build_exact_program` → `poly_to_singular_exact` (safe lcm serializer, §2.2). The corrupt mod-p call is a *screen* only |
| `bridge_sweep.json` → `corner\|R9_z0` | 1 | `verdict = "KILLED"` requires `exact["verdict"] == "UNIT"`; the char-0 path uses `fsb._to_singular` → `_clear_int` (safe). Recorded exact verdict is `UNIT`, wall 3.36 s |
| `altdefect0\|a11_b3100_T2\|5`, `altdefect0\|a12_b1110_T2\|5` | 2 | Two independent sources (`D2_THRESHOLD.md` + `triage_harvest.json`), both safe; one carries `CERTIFICATE-FOUND` |
| All 20 `CERTIFICATE-FOUND` certificates | 20 | Built and checked by `kill_certificate_tools` (`primitive` + `poly_text`, exact, §2.2). `PROOF_PATTERNS.md` records `expansion sum c_i g_i == 1 verified : 20/20`. **This path never touches either defective serializer**, so these are safe irrespective of how the state was first found |
| `CURRENT_STATUS.md` §1 proven / §1b checker-enforced nodes (C6, C11, C12, C14–C16, C18, C21–C29, C33–C34, C41–C46) | all | Grep of `CURRENT_STATUS.md` finds **no** reference to `msolve_bridge`, `triage_harvest`, `bridge_sweep`, `blowup_sweep`, `modular_triage` or `blowup_diagnosis`. The suite checkers are a disjoint lane. C11's Lean-checked certificate comes from `f37_sat_verify.py` / `generators.json` |

### 4.3 UNKNOWN — honestly open

- **Completeness of `triage_harvest.run_a8`.** Its mod-p screen ran on corrupt
  polynomials and *selected* which states got the exact-Q confirmation. A
  corrupt screen cannot make a bad kill land (the exact confirmation is
  independent and sound), but it can have **missed** states that should have
  flipped to `UNIT`. `a8_flip_to_unit: 7 / 13 screened`. That census is unsound
  as a census. No landed claim depends on it; a re-screen could yield *more*
  kills, never fewer.
- **The mod-p reconnaissance corpus.** `modular_triage.json` (systems 1–4), the
  mod-p columns in `BRIDGE_SWEEP.md`, `ALT_BRIDGE.md`, `FULL_SYSTEM_BRIDGE.md`,
  `ALT_ELIM.md`. All computed on corrupted input per §3.2. All are labelled
  predictions/triage, never proofs. **Not re-run.**
- **`j6_msolve.py`'s recorded timings** in `J6_MSOLVE.md` §2 are timings of the
  truncated systems. `JET_OBSTRUCTION.md` §6 already established that the J6
  *verdicts* survive; the *costs* recorded there are not costs of the intended
  systems.

---

## 5. Re-run with the fixed serializers

All five landed `msolve_bridge` cases, re-run over ℚ with the corrected
`sing_poly_intcoeff`, from cached generator lists (so the reconstruction cost is
not repeated). `msolve_bridge_results.json` was **not** written to; output went
to scratch.

| case | pre-fix verdict / wall | **post-fix verdict / wall** | msolve output |
|---|---|---|---|
| `sub2_s38` | EMPTY(KILL) / 3.3 s | **EMPTY(KILL)** / 0.4 s | `[-1]:` |
| `sub2_s14` | EMPTY(KILL) / 6.2 s | **EMPTY(KILL)** / 0.9 s | `[-1]:` |
| `sub2_s94` | EMPTY(KILL) / 15.7 s | **EMPTY(KILL)** / 1.0 s | `[-1]:` |
| `a11_b1111_T1_17` | EMPTY(KILL) / 65.2 s | **EMPTY(KILL)** / 16.1 s | `[-1]:` |
| `a12_b1110_T2_d6` | EMPTY(KILL) / 231.2 s | **EMPTY(KILL)** / 53.1 s | `[-1]:` |

**5/5 verdicts unchanged.** `[-1]` is msolve's "no solution in ℚ̄".

Two observations worth recording.

- The corrected systems are **4–8× faster**, uniformly. The truncated systems
  were *harder*, not easier — truncation destroyed structure (the leading-term
  collapse turns separable relations into nilpotent ones) and made F4 work
  harder for the same answer. This runs opposite to `JET_OBSTRUCTION.md` §6,
  where the corrected depth-2 J6 system was *slower* (1.9 s → 233.1 s). Cost is
  not a reliable tell in either direction; that is why the coefficient-level
  diff of §3.1, not the timings, is the evidence.
- That all five survive is a genuine result, not a foregone conclusion. Emptiness
  of a truncated system implies nothing about the original: the ideals are
  incomparable. Had any come back `HAS_SOL`, a landed kill would have been
  retracted today.

---

## 6. Regression guard

`serializer_roundtrip_verify.py` — `--quiet`, exit 0 on success.

It feeds a 10-polynomial rational corpus (including the two J6 polynomials from
`JET_OBSTRUCTION.md` §1, negative rationals — `int()` truncates toward zero, so
sign matters — sub-unit coefficients that vanish entirely pre-fix, and two
integral controls where the bug is provably inert) through **both** serializers
at char 0, mod 10007 and mod 100019, and **re-parses the emitted string**,
comparing coefficient vectors against expectations computed independently from
`sp.Poly` dictionaries. The serializers' own arithmetic is never trusted.

Three parts, all required for exit 0:

- **A** — live serializers round-trip on the corpus. `10 × 3 = 30` checks.
- **B** — a *verbatim reimplementation of the pre-fix code* is pushed through the
  **same checker**, which must report it as FAILING. This is what proves the
  guard has teeth: if part B ever passes, the checker has stopped detecting the
  original defect and is worthless. Currently catches 7 char-0 and 7 mod-p
  corruptions (the 3 integral corpus entries are skipped — the bug is inert on
  them by definition).
- **C** — with `sp.together` monkeypatched to the identity (precisely how
  `sp.cancel` behaved under sympy 1.14.0 + python-flint 0.9.0, and the reason
  the clearing step became a no-op), the live functions must **RAISE**, never
  silently truncate. This checks the fail-loud guard added in `e9474a0` is
  actually wired in, not just present in a comment.

Actual run on this machine:

```
== A. live serializers ==            30/30 ok
== B. pre-fix re-implementation must be DETECTED ==
  caught 7 char-0 and 7 mod-p pre-fix corruptions
== C. with sp.together neutralised the live code must RAISE ==
  ok   blowup_diagnosis.sing_poly_intcoeff  raised ValueError: non-integer coefficient 1/4 …
  ok   modular_triage.poly_to_singular_modp raised ValueError: non-integer coefficient 1/4 …
SERIALIZER ROUND-TRIP: PASS (10 polynomials x 3 encodings, pre-fix behaviour
detected, fail-loud guard live)
EXIT=0
```

Part C is the part that would have caught this bug on day one: the pre-fix code
had no guard at all, so a silent truncation was indistinguishable from success.

---

## 7. Honesty

- **The comfortable answer was not the true one.** The task anticipated
  "nothing landed is affected; every rational-input path was already
  denominator-cleared upstream". That is true for 385 of 389 ledger kills and
  for every proven/checker-enforced node — but **false for four ledger kills**,
  which were computed from genuinely corrupted ideals and are intact only
  because they were re-run today. The upstream-clearing hypothesis held for
  `triage_harvest`, `bridge_sweep` and `full_system_bridge` (via *safe sibling
  serializers*, §2.2 — not via upstream clearing) and failed for
  `msolve_bridge`, which had no safe path.
- **The certificates did not save the affected kills.** The task suggested
  expansion-verified certificates would cover the exposed set. They cover a
  different set: all five `msolve_blowup__*` certificates are
  `NOT-YET-CERTIFICATED` (lift timed out at 900 s). The 20 `CERTIFICATE-FOUND`
  files are real and safe, and they do rescue 10 of the exposed proof-DAG nodes
  — but not the four in §4.1.
- **`SCALED = 0` is a measurement, not a theorem.** I tested scalar-multiple
  equivalence per generator. I did **not** test the weaker and more relevant
  property "the two ideals are equal" — that would need a Gröbner comparison per
  case. `SCALED = 0` proves the *generators* differ; it does not by itself prove
  the *ideals* differ. The re-run of §5 is what makes this moot for the landed
  claims, and it is why I re-ran rather than arguing.
- **I did not re-run the mod-p reconnaissance** (§4.3, second bullet). Those
  files remain stale-on-corrupt-input and are not marked as such in their own
  docs. That is a known, unrepaired gap, not an oversight.
- **The `run_a8` completeness gap (§4.3, first bullet) is unresolved.** It can
  only cost us kills we never claimed, so it is not a soundness risk, but the
  `7/13` flip census in `TRIAGE_HARVEST.md` should not be cited as a census.
- **Scope of the "no ledger kill" claims in §2.1.** For callers 9–11 I verified
  that no ledger source and no proof-DAG node names them. I did not audit
  whether some prose claim in `ALT_BRIDGE.md` or similar rests on them.
- **A concurrent lane was committing during this audit.** `HEAD` moved from
  `124a733` to `79b4fa7` mid-session, and that lane's `git add` swept
  `serializer_roundtrip_verify.py` into commit `79b4fa7` (whose message is about
  the discriminant identity and does not mention it). The committed copy is
  byte-identical to what this audit wrote and still exits 0. That same lane also
  added an independent fail-loud guard to `saturated_cell.poly_to_singular_exact`
  — I did not write that change and have not reviewed it. `SERIALIZER_BUG.md` is
  left uncommitted.
- **One correction to `JET_OBSTRUCTION.md` §8.** Its caller list omits
  `triage_harvest.py`, which reaches `poly_to_singular_modp` transitively
  through `modular_triage.build_singular_program`. That omission mattered:
  `triage_harvest` is the single largest exposed contributor to the ledger (11
  kills), and it took reading the `exact_kill` gate to establish they are safe.

---

## 8. Reproduction

```
python serializer_roundtrip_verify.py            # guard, exit 0 iff healthy
python serializer_roundtrip_verify.py --quiet    # same, silent
```

The blast-radius measurements of §3 and the re-runs of §5 were produced by
throwaway scripts in the session scratchpad (generator caching, coefficient-vector
diffing, msolve re-invocation). They are not committed: everything they establish
is either recorded in the tables above or re-derivable from
`serializer_roundtrip_verify.py` plus `msolve_bridge.run_msolve`, which now uses
the fixed serializer.

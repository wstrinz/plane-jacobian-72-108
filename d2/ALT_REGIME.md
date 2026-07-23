# Alternate regime for f31/subcase (1)

Date: 2026-07-22  
Verifier: `alt_regime_verify.py`

## Flipped reduction

> **Post-L2 update (audit note N2):** the '6 strata die in both branches' figure below is the pre-L2 count; after `ALT_REGIME_L2.md` it is **10**.


Write `Phi~=t^30u`, `u=cq`, `e=t^a E`, and `E(-1)!=0`. For `a=11,...,15`, put

```text
v=30-3a<0,              w=|v|=3a-30.
```

The t-order of the f-term of

```text
F=sum_{f=0}^7 Phi~^f e^(21-3f) h_f(d~)
```

is `30f+a(21-3f)=21a+fv`. It is minimized at `f=7`, and

```text
21a+7v=21a+7(30-3a)=210.
```

Hence the alternate polynomial reduction is

```text
F=t^210 G',
G'=sum_{f=0}^7 t^((7-f)w) u^f E^(21-3f) h_f(d~).
```

The verifier checks this exactly on a seeded random subcase-(1) window with `a=12`, and checks every exponent for all five values of `a`.

## Descending cascade

Put `T=t^w`. A normalized descending recursion is

```text
h_7                         = T r_6,
E^3  h_6 + u r_6            = T r_5,
E^6  h_5 + u r_5            = T r_4,
             ...
E^18 h_1 + u r_1            = T r_0,
E^21 h_0 + u r_0            = 0.
```

It telescopes exactly to `G'=0`. The first forced condition is `T|u^7h_7`. Since `(T,u)=1` and `h_7=8192d1^2`,

```text
2v_t(d1)>=w.                                      (T1 anchor)
```

The terminal condition is now the bottom equation `E^21h_0+ur_0=0`, rather than a top equation reached by upward t-divisibility.

For T2 (`d1=0`), `r_6=0` and

```text
T r_5 = E^3 h_6 = -3072 E^3 sigma^2,
T r_4 = E^6 h_5 + u r_5,
h_5|_{d1=0} = -9216 d2 sigma^2 + 2048 e^2.
```

At `t`, the two terms of `h_5` are divisible by `T` after the anchor: `2v_t(sigma)>=w` and `2a>=w`. Thus `T|r_5`, giving the stronger conclusion

```text
v_t(sigma)>=w.                                    (T2 first level)
```

## Survival/replacement table

| Standard tool | Alternate-regime verdict |
|---|---|
| Graded identity `f31=sum Phi^f e^(21-3f)h_f` | Survives verbatim; it is window-independent. |
| `F=t^(21a)G` | Survives only as a Laurent equality. Its polynomial replacement is `F=t^210G'`. |
| Upward t-cascade from `h_0` | Fails; replaced by the descending `r_6,...,r_0` cascade. |
| Polynomial `g_l` and top terminal identity | Survive after defining `g_1=Th_0`, `g_{l+1}=T(E^3g_l+u^lh_l)`. Then `E^3g_7+u^7h_7=G'`. |
| q-place valuation transitions | Survive verbatim: `t`, hence `T`, is a unit at every q-root. |
| T1 q-terminal | Survives: `3b_i+v_i(g_7)=7+2v_i(d1)`. |
| T2 q-terminal | Survives: `3b_i+v_i(g_6)=6+2v_i(sigma)`. |
| Terminal caps | `deg g7<=46`, `deg g6<=48` survive, sharpening in a stratum to `46-3sum b_i`, `48-3sum b_i`. |
| Old lower-level global `g_l` caps / t coupling | Do not transfer; negative `v` reverses the degree recursion and t-adic edge. |

Thus the local q equations are regime-independent, but the old upward t-divisibility proof is not. The bottom-up polynomial definition supplies the q-local auxiliary equations without importing the invalid t-side.

## Terminal plus first-level local lemmas

Let `P` be `t` or a simple q-root, `m=v_P(e)`, and let the flipped spacing be `s=w=3a-30` at `t`, or `s=3b-1` at a q-root with `b>0`.

For T1 put `x=v_P(d1)`. The anchor gives `2x>=s`, while

```text
r_6=8192d1^2/P^s,
h_6=-3072sigma^2+14336d1^2d2+8192d1e.
```

If `s` is odd and `x<s`, `v_P(r_6)=2x-s` is odd. The `sigma^2` order is even and the other two h6 terms have strictly larger order, so the first congruence cannot cancel. Hence `x>=s` for odd `s`; for even `s`, this level gives `x>=s/2`. Therefore

```text
T1 at t: a=11,12,13,14,15 -> v_t(d1)>=3,3,9,6,15.
T1 at q: b=1,2,3,4       -> v_P(d1)>=1,5,4,11.
```

For T2 put `z=v_P(sigma)`. If `2m>=s`, both h5 terms are divisible by `P^s`, so the next congruence forces `z>=s`. If `2m<s`, h5 has exact order `2m`, which must equal `2z-s`. Odd `s` is then impossible; for even `s`, `z=(s+2m)/2`. Thus

```text
T2 at t: a=11,12,13,14,15 -> v_t(sigma)>=3,6,9,12,15.
T2 at q: b=1 -> v_P(sigma)>=2;
          b=2 -> impossible;
          b=3 -> v_P(sigma)=7 is necessary, with residue cancellation;
          b=4 -> impossible.
```

Orders at distinct places add toward polynomial degree. Combining these with `deg d1<=9`, `deg sigma<=12`, the surviving q-terminal identities, and `deg E<=15-a<=4` gives the exhaustive audit below.

## Per-stratum verdicts

`K` means proved killed at terminal/first level; `O` means open. An open branch passes these valuation and degree tests but still has residue congruences and lower flipped levels to solve.

| `a` | sorted `(b1,b2,b3,b4)` | T1 | T2 |
|---:|:---:|:---:|:---:|
| 11 | `(0,0,0,0)` | O | O |
| 11 | `(1,0,0,0)` | O | O |
| 11 | `(1,1,0,0)` | O | O |
| 11 | `(1,1,1,0)` | O | O |
| 11 | `(1,1,1,1)` | O | O |
| 11 | `(2,0,0,0)` | O | K: q parity |
| 11 | `(2,1,0,0)` | O | K: q parity |
| 11 | `(2,1,1,0)` | K: `deg d1>=10` | K: q parity |
| 11 | `(2,2,0,0)` | K: `deg d1>=13` | K: q parity |
| 11 | `(3,0,0,0)` | O | O (`deg sigma>=10`) |
| 11 | `(3,1,0,0)` | O | O (`deg sigma>=12`) |
| 11 | `(4,0,0,0)` | K: q order 11 | K: q parity |
| 12 | `(0,0,0,0)` | O | O |
| 12 | `(1,0,0,0)` | O | O |
| 12 | `(1,1,0,0)` | O | O |
| 12 | `(1,1,1,0)` | O | O (`deg sigma>=12`) |
| 12 | `(2,0,0,0)` | O | K: q parity |
| 12 | `(2,1,0,0)` | O | K: q parity |
| 12 | `(3,0,0,0)` | O | K: `deg sigma>=13` |
| 13 | `(0,0,0,0)` | O (`deg d1>=9`) | O |
| 13 | `(1,0,0,0)` | K: `deg d1>=10` | O |
| 13 | `(1,1,0,0)` | K: `deg d1>=11` | K: `deg sigma>=13` |
| 13 | `(2,0,0,0)` | K: `deg d1>=14` | K: q parity |
| 14 | `(0,0,0,0)` | O | O (`deg sigma>=12`) |
| 14 | `(1,0,0,0)` | O | K: `deg sigma>=14` |
| 15 | `(0,0,0,0)` | K: `deg d1>=15` | K: `deg sigma>=15` |

Totals:

- T1: 7 killed, 19 open.
- T2: 12 killed, 14 open.
- Overall: **19 of 52 branches killed**.
- Six strata die in both branches; **33 branches remain open in 20 strata**.

For every open branch, the next exact obligations are the residue congruences

```text
E^3h_6+ur_6=Tr_5                         (T1),
E^6h_5+ur_5=Tr_4                         (T2),
```

followed by the lower descending levels through `E^21h_0+ur_0=0`. The old standard-regime lower-level degree profiles cannot be used here.

## Proposed STATE.md entry

```text
- SUBCASE-(1) ALTERNATE REGIME (2026-07-22): `ALT_REGIME.md` and
  `alt_regime_verify.py` derive and exactly check the v<0 reduction for
  a=11..15. Since 21a+7(30-3a)=210, F=t^210 G' with
  G'=sum t^((7-f)|v|)u^f E^(21-3f)h_f, and the t-cascade descends from
  h7=8192d1^2. The old q-place transitions and terminal identities survive
  via polynomial bottom-up auxiliaries, but the old upward t-cascade and its
  lower-level degree profiles do not. Exact terminal+first-level parity
  lemmas, coupled to deg d1<=9 and deg sigma<=12, kill 19/52 branches
  (T1 7, T2 12); 6/26 strata die in both branches. Frontier: 33 branches in
  20 strata, each with explicit first residue congruences and the lower
  descending cascade still open. Verifier passes all checks.
```

#!/usr/bin/env bash
# Repo test runner (public release layout).
#   D3: fast, exact — always run.
#   D2: exact proof checkers — always run.
#   D2 harness positive control (~a few minutes) — run unless SKIP_SLOW=1.
set -uo pipefail
here="${RUN_TESTS_HERE:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
if [ -z "${RUN_TESTS_PRIVATE_COPY:-}" ]; then
  _tmp="$(mktemp)" && cp "${BASH_SOURCE[0]}" "$_tmp" &&
  RUN_TESTS_PRIVATE_COPY=1 RUN_TESTS_HERE="$here" exec bash "$_tmp" "$@"
fi
rc=0

echo "### CLEAN-CLONE guard — every file the suite reads must be git-tracked"
( cd "$here" && python3 tools/clean_clone_check.py ) || rc=1

echo "### D3 — exact counterexample verification"
( cd "$here/d3" && python3 verify.py ) || rc=1

echo "### D2 — exact field-split proof checks"
(
  cd "$here/d2" &&
  python3 t5_split_place_verify.py &&
  python3 test_split_place_proofs.py &&
  python3 test_split_place_ledger.py &&
  python3 test_cascade_signature.py &&
  python3 test_cascade_engine.py &&
  python3 test_cascade_inf.py &&
  python3 cascade_inf_ties_verify.py &&
  python3 test_cone_lemmas.py &&
  python3 audit_cascade_kills.py --quiet &&
  python3 audit_cascade_kills_sub1.py --quiet &&
  python3 audit_tplace_cases.py --quiet &&
  python3 audit_inf_cases.py --quiet &&
  python3 audit_alt_regime.py --quiet &&
  python3 audit_convolution_kills.py --quiet &&
  python3 audit_convolution_kills_r2.py --quiet &&
  python3 audit_reconstruction_kills.py --quiet &&
  python3 residue_lemmas_depth_verify.py &&
  python3 t5_90t2_verify.py &&
  python3 t5_90t1_verify.py &&
  python3 t5_90t1_constant_verify.py &&
  python3 convolution_descent.py &&
  python3 convolution_elim.py --gates-only &&
  python3 t5_90t1_local_verify.py &&
  python3 f37_free_family_verify.py &&
  python3 f37_sat_verify.py &&
  python3 sub1_cascade_verify.py &&
  python3 alt_regime_verify.py &&
  python3 t6_premises_verify.py &&
  python3 t5_t2_column_verify.py &&
  python3 t5_t2_infinity_verify.py &&
  python3 alt_regime_l2_verify.py &&
  python3 residue_lemmas_verify.py &&
  python3 alt_regime_audit_verify.py &&
  python3 alt_regime_inf_verify.py &&
  python3 alt_inf_sweep_verify.py &&
  python3 alt_combined_verify.py &&
  python3 alt_residue_congruences_verify.py &&
  python3 phase_f2_pilot_verify.py &&
  python3 envelope_bounds_verify.py &&
  python3 phi_75_125_verify.py &&
  python3 phi_corner4_verify.py &&
  python3 phi_f14_verify.py --quiet &&
  python3 galois_library_verify.py --quiet &&
  python3 phi_f7_verify.py --quiet &&
  python3 prior_art_postdiction_verify.py --quiet &&
  python3 composite_charts_verify.py --quiet &&
  python3 zeta_tail_verify.py --quiet &&
  python3 mu_rungs_f10_verify.py --quiet &&
  python3 c_series_75_125_verify.py --quiet &&
  python3 alok_crosscheck.py --quiet &&
  python3 g_system_75_125_verify.py --quiet &&
  python3 f2_family_verify.py --quiet &&
  python3 f2_tower_verify.py --quiet &&
  python3 family_grammar_verify.py --quiet &&
  python3 chain_survey_verify.py --quiet &&
  python3 polygon_reduction_verify.py --quiet &&
  python3 minimal_core.py --quiet &&
  python3 support_certificates.py --quiet &&
  python3 contact_lemma.py --quiet &&
  python3 toric_syzygy.py --quiet &&
  python3 toric_general.py --quiet &&
  python3 corner_atlas.py --quiet &&
  python3 corner_atlas_toric_34.py --quiet &&
  python3 gamma_from_corner_verify.py --quiet &&
  python3 monomial_window_law.py --quiet &&
  python3 ml_restriction_check.py --quiet &&
  python3 case_compiler_verify.py --quiet &&
  python3 window_caps_verify.py --quiet &&
  python3 window_functions_75_125_verify.py --quiet &&
  python3 proof_dag_report.py --quiet &&
  python3 phi_depth_criterion_verify.py --quiet &&
  python3 marked_polynomial.py --quiet &&
  python3 serializer_roundtrip_verify.py --quiet
) || rc=1

if [ "${SKIP_SLOW:-0}" = "1" ]; then
  echo "### D2 — harness control SKIPPED (SKIP_SLOW=1)"
else
  echo "### D2 — harness positive control (f31_sub2, must reach <=1e-5)"
  ( cd "$here/d2" && python3 jetlift.py control f31_sub2 ) || rc=1
fi

echo
[ "$rc" = "0" ] && echo "ALL TESTS PASSED" || echo "TEST FAILURES (rc=$rc)"
exit $rc

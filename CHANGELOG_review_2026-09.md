# Review fixes — branch `review/2026-09-line-by-line-fixes`

Base: `e6644e2` (main, 2026-01-28). Line-by-line code/logic review of all 25 files (177 findings: 2 critical, 44 major, 107 minor, 24 style); every finding rated critical/major/minor/style was addressed except where noted in the per-group notes (`fix_*_notes.md`, kept outside the repo).

## Behaviour-changing fixes (read before re-running old projects)

1. **FASTQ pairing (`01_preprocess_reads.pl`)** — the `<base>_1/_2.fastq(.gz)` rule is now tested *before* the bcl2fastq `_R1/_R2` rule, and a read number assigned twice to the same sample aborts with an explicit message. Sample names containing `_R1`/`_R2` (e.g. `NEG_R1`, `NEG_R2`) were previously mis-paired across samples.
2. **fastp failures propagate** — `01_preprocess_reads.pl` now exits non-zero if any sample's fastp run fails, and never writes the completion marker in that case.
3. **Stage completion markers** — `01_fastp_out/.preprocess_done` and `02_decoded/.decode_done` are written only on success and removed at stage start; `run_delegance_pipeline.py` uses them for resume decisions (legacy `All done.` log lines are still accepted for old runs).
4. **`run_delegance_pipeline.py`** — `--skip-fastp` pre-checks no longer require `BB_information_fixed.tsv` before preprocess has run; `preprocess_done` is computed before it is used; decode/hit cache hashes are computed after preprocess; `--no-stop-on-error` returns exit 1 when any stage failed; `--dry-run` no longer moves legacy files; the beginner QC report is written as `beginner_qc_report.html` (no longer overwrites 03's `report.html`).
5. **`03_call_hits.py`**
   - `--auto_tune 0` is honoured (pre- and post-GLM tuning skipped; CLI values used); `--auto_lock` applies to both phases; parameters actually used are written to `OUTDIR/hit_params_used.json` (`--auto_save_json`), and a minimal `hit_params.json` is written for stand-alone runs.
   - The pseudocount is applied to the observed counts as well as the DEL2 offset, so zero-count conditions give the bounded floor `log2(epsk/(DEL2+epsk))` on every backend (previously −20…−45 depending on numpy/torch), which also stabilises NEG centering.
   - `--dispersion shared` no longer breaks out of the condition loop after the first condition.
   - ReadScaler statistics (`avg_R1`, `mean_log2FC_*`, `q_*`, `HitScore_RS`) are computed on size-factor-normalised counts instead of raw counts.
   - q-value gates (`max_q_bead`, `max_q_bead_r2`, `max_q_paired_boost`; set by the `neg_heavy`/`very_hard` presets) are disabled with a warning when the condition has fewer than 4 replicates (no statistical power; previously every compound failed); effect-size gates and `neg_gate_mode` are unchanged. Add the parameter name to `--auto_lock` to keep the gate.
   - NB2 log-pmf switches to the Poisson limit for alpha < 1e-4 (float32-safe); one-sample t-test with sd = 0 and mean ≠ 0 now returns p = 0; `report.html` image paths are relative to the report; `tag_seq` whitespace stripping fixed; `read_csv` keeps literal `NA` strings; BB annotation lookup is O(1).
6. **`02_decode_reads.pl`** — `--auto-select-window` / `--auto-center-window` default to 0 and are now effective (`-W 1` selects a single 3/4-cycle window, `-C 1` re-centres it); reads are upper-cased and CR-stripped; 1-bp anchor variants that collide with another anchor are skipped with a warning instead of aborting; `mkdir` via File::Path.
7. **`04_build_interactive_report.py`** — `--img_fallback` defaults to `none` (SMILES are no longer sent to cactus.nci.nih.gov unless `http` is requested explicitly; a warning is printed); the HTML is saved with `save()` instead of `show()` (no browser launch).


### Orchestration group — interface/contract changes

- postprocess_after_hits.sh: 초보자 QC 리포트 출력 파일명 변경 — <preset_dir>/report.html → beginner_qc_report.html, report.tsv → beginner_qc_tophits.tsv (03_call_hits.py의 report.html은 더 이상 덮어쓰지 않음).
- postprocess_after_hits.sh: DEL2_COL 환경변수 기본값 DEL234 → DEL2.
- run_delegance_pipeline.py: 스테이지 완료 판정 마커 파일 신설 — <RUN_ROOT>/01_fastp_out/.preprocess_done, <RUN_ROOT>/02_decoded/.decode_done (내용: ISO 8601 시각). 마커가 없고 로그에 'All done.'이 있는 구 실행은 완료로 취급(하위호환). 01/02 Perl이 성공 시 같은 파일을 쓰도록 리드 검토자 패치 필요; orchestrator도 성공 시 기록함.
- run_delegance_pipeline.py: 03_normalized/<dir>/hit_params.json 스키마 version 1→2 (args가 hit 관련 17개 키만 포함, decode_hash → decoded_matrix 지문). 기존 v1 캐시는 해시 불일치로 1회 재계산됨.
- run_delegance_pipeline.py: 01_fastp_out/preprocess_params.json에서 threads 키 제거 → 기존 캐시는 1회 재실행됨.
- run_delegance_pipeline.py: --no-stop-on-error에서 어느 단계라도 실패하면 종료 코드 1(이전 0).
- run_delegance_pipeline.py: index.html이 <preset>/beginner_qc_report.html 링크를 (존재할 때) 추가하고, 결과가 없는 하위 디렉터리와 존재하지 않는 파일 링크는 생략.
- run_full_validation.py: undecoded 재검사 출력 파일명 undecoded_sample_recheck_<N>.tsv → undecoded_sample_recheck_<sample>_<N>.tsv; JSON에 undecoded_missing/decoded_missing 플래그, validation 카운터 cycles_out_of_range 추가; 새 CLI 옵션 --mismatch, --max-cp-cands, --max-anchor-cands(기본값은 디코더 기본과 동일: hp_op_cp, 0, 0 — 종전 하드코딩 6/5/5보다 완화).
- run_hits_then_postprocess.sh: hit 실행에 --auto-opt 0 --glm-mode full --device cpu --dtype float64 고정 전달(폴더명 glm_full_dev_cpu_fp64와 일치시킴).
- run_autopilot.sh: R1_COLS/R2_COLS/NEG_COLS 환경변수 덮어쓰기 지원(R2_COLS="" → R1-only).

### Reports group — interface/contract changes

- 04 `--img_fallback` 기본값 auto → none (CLI 기본 동작 변경; RDKit 미설치 시 이미지가 생략됨. 기존 동작을 원하면 `--img_fallback http` 명시).
- 04 출력 HTML 생성이 show()→save()로 바뀌어 브라우저를 열지 않음(파일 내용/이름 불변).
- 04 LibID 소스 우선순위: primary_lib_id_used → LIB_ID_x → LIB_ID → LIB_ID_y → lib_id* (기존은 LIB_ID가 LIB_ID_x보다 우선; 03 hybrid에는 둘 중 하나만 있으므로 실제 영향 없음).
- 04 Top-N 정렬이 안정 정렬 + 2차 키(HitScore_RS desc, id asc)로 바뀌어 동점 경계 행 선택이 결정적이 됨(기존 quicksort 결과와 동점 구간 순서가 달라질 수 있음).
- 05 `--top-n 0`(또는 음수): 기존 1개 → 전체 행; 출력 prefix가 `top1_*` 대신 `topall_*`.
- 05 점수(HitScore_GLM/RS) 결측 행은 순위에서 제외됨(기존: -inf로 하위 배치되어 top_n 부족 시 포함).
- 05 `top<N>_diverse.tsv` 선택 규칙 변경: 필터 통과 전체에서 클러스터별 1개(대표 우선) 후 rec_n 컷. 클러스터링 미수행(RDKit 없음/--cluster 0) 시 diverse=recommended(기존: 빈 파일 또는 동일).
- 05 순위 동점 처리: HitScore_RS desc → ID_x asc 2차 키(mergesort). rank 값이 동점 구간에서 기존과 달라질 수 있음.
- 05 내부 컬럼 정규화: 입력에 ID/LIB_ID/BB1..BB4/CP(무접미사)만 있어도 동작(ID_x 등 별칭 생성). 출력 컬럼명은 불변(rep_ID_x, rep_BB1_x ...). 단, top<N>_hits.tsv에는 원본 컬럼과 *_x 별칭이 함께 기록될 수 있음(03이 접미사를 없앤 경우에 한함).
- 05 top<N>_hits.tsv에 *_CPM 및 raw 카운트 컬럼이 추가로 포함됨(usecols 확장).

### Compare group — interface/contract changes

- tier_report_*_all_candidates.tsv / *_specific.tsv / *_diverse.tsv / *_other.tsv 에 컬럼 2개 추가: source_run(행을 제공한 run 라벨), source_role(active|inactive|both). 기존 컬럼명은 변경 없음.
- tier_report: active_score 값 의미가 'active run 의 점수'로 교정됨(README 정의와 일치); rank 는 active_score 기준. 기존 결과와 수치가 달라질 수 있음(inactive/both 전용 후보).
- tier_report: cluster_rep/cluster_medoid 가 (group, cluster) 단위로 재선정되어 *_diverse.tsv 행 수가 기존보다 늘어날 수 있음. cluster_id 는 동일.
- --labels 생략 + 유도 라벨 중복 시 라벨(및 출력 파일명 compare_top<N>_<label>_*.tsv, union/overlap 컬럼 접두어)에 _1,_2 접미사가 붙음.
- --out-prefix 가 디렉터리 성분 없는 bare 이름이면 이제 --out-dir 아래에 생성됨(기존: cwd).
- 신규 CLI 플래그 --bbavg-missing {skip,zero} (기본 skip = 기존 동작). 기존 플래그는 모두 유지(미사용 3개는 help 문구만 변경).
- compare_top<N>_<label>_bb_frequency.tsv 는 --include-summary 0 일 때 더 이상 빈 파일로 생성되지 않음.
- union/overlap 테이블의 <label>_score 에 -inf 대신 NaN 기록; <label>_rank/_rank_pct 계산에서 NaN score compound 제외(기존엔 하위 순위로 포함).
- 입력: ID/LIB_ID/BB1..BB4/CP 가 *_x 없이 plain 이름으로 들어와도 읽힘(읽은 직후 *_x 로 정규화, [WARN] 출력).

### Tiered group — interface/contract changes

- 새 CLI 인자 --del2(기본 DEL2), --exclude-samples(기본 DEL234): 기본값에서 DEL2_CPM이 enrichment 집계에서 추가 제외되어 *_enrich 및 enrichment 기반 group_rank가 달라질 수 있음
- 새 출력 컬럼 active_present/inactive_present/both_present(0/1)가 모든 TSV/XLSX에 추가
- group_rank_score, <score_col>, active_score/inactive_score/both_score의 -inf가 NaN(빈칸)으로 기록; 동점 score의 *_rank/*_rank_pct가 method='min'으로 동일값
- coalesced 샘플 컬럼 값의 우선순위가 항상 active→inactive→both(이전: 행 출처 run 우선)
- --both-run 지정 시 <prefix>_all_candidates_<both_label>.xlsx 추가 생성
- --labels 생략 + 중복 라벨 시 _2/_3 접미사로 prefixed 컬럼명·group 표시명이 바뀔 수 있음
- NEG 필터: percentile 0인 컬럼은 NEG==0 행 유지(이전 전부 제거); 필터 후 후보 0이면 exit 1
- 파일명·기존 컬럼명·기존 CLI 플래그명 변경 없음

### Utilities group — interface/contract changes

- make_top_hit_shift_report.py 출력 TSV/Excel에 컬럼 `in_cur_file`, `in_prev_file`(bool) 추가; 상대 파일에 없는 compound도 행으로 유지되어 행 수가 늘어날 수 있음(해당 행의 *_prev/*_cur 값은 NA).
- export_full_runs.py `compound_key` 정의 변경: `CP_x|BB1_x|...|BB4_x`(LIB 유지, NaN→'nan') → `BB1|BB2|BB3|BB4`(LIB 접미사 제거, CP 제외, 빈 BB→''). 06/07 원본은 빈 BB를 'NA'로 표기하므로 4-cycle이 아닌 라이브러리에서 키를 맞추려면 `KEY_NA` 상수와 06/07 구현을 한쪽으로 통일해야 함.
- verify_random_reads.py: count 검증 키가 `id` → `(lib_id, id)`로 변경; stdout 예시 튜플 형식 변경. raw_counts_matrix.tsv에 `lib_id` 컬럼이 필수(02_decode_reads.pl은 항상 기록).
- export_final_excel.py: 시트 이름 규칙 변경 — 접미사(`_All_Core`, `_Top{N}`, `_Consensus`, `_Params`)를 보존하고 run 부분을 절단; 충돌 시 `~N` 접미사. 31자 미만 run 이름은 기존과 동일.
- export_beginner_qc_report.py: `--del2_col` 미지정 시 annot 옆 `hit_params.json`이 있으면 그 `normalized_columns.del2`를 사용(orchestrator 실행 결과에서 DEL2_raw 컬럼이 raw count로 바뀔 수 있음). 휴리스틱은 `*_norm/*_sum/*_CPM`을 더 이상 선택하지 않으며, DEL 이름 컬럼이 전혀 없으면 DEL2_raw 컬럼·LOW_DEL2 플래그가 생략됨.
- subsample_fastq_pairs.py: 인식 파일명 범위 확대(`-R1`, `_R1.<tok>`, `_R1_L001` 등); 출력 파일명 규칙(`{base}{suffix}{read}{ext}`)은 그대로. head 모드 로그 라벨 `total_pairs=` → `pairs_read=`.
- anonymize_with_map.py: `--dry-run` 옵션 추가; `.git/` 등 제외; 실패가 있으면 exit code 1.
- analyze_lib_suffix_missing.py: pandas 의존 제거(csv 모듈로 동일 TSV 출력).
- environment.yml: openpyxl/networkx/perl 추가, pandas>=2.0·bokeh>=3.4 pin.

## Not changed on purpose

- The `_x/_y` column suffixes in `05_hybrid_annot.tsv` are kept; downstream scripts (04–07, exporters) now accept both `ID`/`ID_x` etc., so the suffixes can be removed from 03 in a later release without breaking them.
- Preset threshold values are unchanged; only the q-gate disabling at n<4 replicates was added.
- The logistic slope heuristic constant (2.772) is kept for backward compatibility (comment added).
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
- export_full_runs.py `compound_key` 정의 변경: `CP_x|BB1_x|...|BB4_x`(LIB 유지, NaN→'nan') → `BB1|BB2|BB3|BB4`(LIB 접미사 제거, CP 제외, 빈 BB→'NA'; 06/07 `_make_compound_key`와 동일).
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

## 2차 검토(2026-09-06) 수정

1차 수정 트리를 다시 라인-바이-라인 검토하여 발견한 회귀(R)와 주요 잔여 결함(M)의 수정.

### 회귀 (1차 수정으로 도입)

- R1 `03_call_hits.py` — one-sample t-test에서 sd = 0, mean ≠ 0일 때 p = 0으로 바꾼 변경을 되돌려 p = 1(정보 없음)로 복귀; n < 2 도 p = 1.
- R2 `03_call_hits.py` — 관측 카운트에 pseudocount를 더하던 방식을 제거하고 `fit_one(delta_min=)`으로 log2FC 하한을 클램프(0-카운트 조건의 floor는 유지, 비-0 조건의 추정치는 왜곡하지 않음).
- R3 `03_call_hits.py` — `bh_fdr`가 NaN p-value를 유지(NaN을 1로 치환하지 않음; 순위·분모에서 제외).
- R4 `03_call_hits.py` — NEG 센터링 null subset에서 NEG 원카운트가 0인 행을 제외(`_neg_nonzero_mask`); n < 4 q-gate 규칙에서 `< 0.2` 조건 제거, 부재 조건에 대한 경고 억제; `hit_params_used.json` 키 확장.
- R5 `01_preprocess_reads.pl` — bcl2fastq `_R1/_R2` 규칙의 base를 greedy로(마지막 `_R[12]` 토큰에서 분리), `--fastq-regex` 중복 배정 가드, 짝 없는 파일·0바이트 파일·페어 0개는 실패로 처리, fixed BB 파일의 `type`을 대문자로 기록, 헤더 판정은 첫 non-blank 줄 기준.
- R6 `subsample_fastq_pairs.py` — 01과 동일한 2단계 페어링 규칙(`match_read_name`/`READ_RE`)으로 통일; `out_r1 == out_r2` 가드; 중간 빈 줄 skip.
- R7 유틸리티 5종(`export_beginner_qc_report.py`, `export_final_excel.py`, `export_full_runs.py`, `make_display_hybrid_tsv.py`, `make_top_hit_shift_report.py`) — 1차에서 `_LIB[^_]+$`로 바뀐 BB 접미사 정규식을 04–07과 같은 `_LIB[\w.-]+$`로 복귀(`_`를 포함한 lib_id에서 접미사가 제거되지 않던 회귀). 전체 ID용 fallback은 LIB_ID 컬럼에서 확인된 lib_id의 정확한 `_LIB<lib_id>` 토큰만 제거하고, lib_id를 전혀 알 수 없을 때만 `_LIB[^_]+(?=_|$)` 정규식을 사용(export_beginner/make_display/make_top_hit_shift 공통).

### 주요 잔여 결함

- M1 `export_full_runs.py` — KEY_NA 주석·docstring의 빈 BB→'' 서술을 코드 동작('NA', 06/07과 동일)에 맞게 정정; `ID/LIB_ID/CP/BB1..BB4`가 `_x` 접미사 없이 들어와도 읽은 직후 `_x`로 정규화(`[WARN]`), BB 컬럼이 전혀 없으면 compound_key가 전부 `NA|NA|NA|NA`로 붕괴하므로 exit 1; `*_all.tsv`/`merged_all.tsv`에 `na_rep="NA"`; `--active-label`과 `--inactive-label`이 같으면 exit 1.
- M2 `export_beginner_qc_report.py` — `--del2_col`이 컬럼에 없으면 무경고로 휴리스틱에 떨어지던 것을 `[WARN]` 후 `hit_params.json` → 휴리스틱 순서로 폴백(postprocess_after_hits.sh가 항상 `--del2_col`을 전달하므로 실질적으로 hit_params.json 경로가 살아남); `_qc_flags`가 쓰는 LFC 컬럼을 한 번 `to_numeric(errors="coerce")`로 변환(문자 값에서 ValueError 방지); NaN LibID(`NA`) 행은 fallback ID 사용; HTML 플래그 목록에 NEG_R1_HIGH 추가, 임계값 노트를 항목별로 표기.
- M3 `export_final_excel.py` — `hit_params.json`을 UTF-8·예외 처리로 읽음(손상 JSON은 `[WARN]` 후 무시); 03 단독 실행의 최소 `hit_params.json`(`args` 없음)에서는 빈 `_Params` 시트를 만들지 않고 `[INFO]`만 출력; LibID 결측 행의 `ID`가 `nan_...`이 되지 않도록 BB-only ID로 폴백; 헤더 폴백 `[WARN]` 문구가 파일 부재/`normalized_columns` 부재를 구분.
- M4 `make_top_hit_shift_report.py` — 표시 필드 Series를 df.index로 생성(필터된 프레임에서도 정렬 유지); 결측 LibID를 `nan` 문자열 대신 `NA`로 출력.
- M5 `anonymize_with_map.py` — 매핑 파일의 `.json`/`.tsv` 쌍도 rename 보호; `--dry-run` 출력이 부모 디렉터리 rename 전 경로임을 명시.
- M6 `README.md` — Beginner QC 산출물명을 `beginner_qc_report.html`/`beginner_qc_tophits.tsv`로 정정(`report.html`은 03의 all-in-one 리포트); `top<N>_hits.tsv`가 실제로 추가하는 컬럼(`rank` + cluster_*)만 기재하고 `compound_key`/`rank_pct`/`score_z` 서술 삭제; 05 diverse 선택 규칙, 06 `source_run`/`source_role`·`--bbavg-missing`·라벨 중복 접미사·`--out-prefix` 규칙·NaN 점수, 07 `--del2`/`--exclude-samples`·`*_present`·NEG 필터·`_all_candidates_<both>.xlsx`, 04 `--img_fallback` 기본값, 01 FASTQ 페어링 규칙(`_R1` 포함 샘플명, `--fastq-regex`), run_full_validation `validation_cycles_out_of_range`·`undecoded_sample_recheck_<sample>_<N>.tsv` 및 `--mismatch/--max-cp-cands/--max-anchor-cands`가 디코드 설정과 일치해야 한다는 주의를 반영.
- M7 `CHANGELOG` — Utilities 그룹의 export_full_runs 항목(빈 BB→'')을 실제 동작('NA')으로 정정.
- M8 정적 분석 확인 — pyflakes의 `export_final_excel.py`/`make_top_hit_shift_report.py` openpyxl·xlsxwriter 'imported but unused'는 엔진 존재 확인용 import(`# noqa: F401`)로 오탐; ruff 목록에 유틸리티 항목 없음.

### 2차 수정의 인터페이스/계약 변경 (Utilities)

- make_top_hit_shift_report.py: `LibID` 컬럼의 결측이 `nan` → `NA`. 표시 `ID` fallback은 테이블에 존재하는 lib_id의 `_LIB<lib_id>` 토큰만 제거(알 수 없는 토큰은 유지).
- export_full_runs.py: `*_all.tsv`/`merged_all.tsv`의 NaN이 빈 칸 → `NA`; `_x` 없는 키 컬럼은 `_x`로 rename되어 기록됨; 라벨 동일·BB 컬럼 부재 시 exit 1.
- export_final_excel.py: `args` 없는 hit_params.json이면 `_Params` 시트가 생성되지 않음(이전: 빈 시트); LibID 결측 행의 `ID`/`ID_display`가 `BB1_BB2_BB3_BB4` 형식.
- export_beginner_qc_report.py: 존재하지 않는 `--del2_col`은 hit_params.json 값으로 대체됨(이전: 휴리스틱) — postprocess_after_hits.sh 경유 결과의 `DEL2_raw` 컬럼이 바뀔 수 있음. HTML 임계값 노트 문구 변경.

### 2차 수정의 인터페이스/계약 변경 (그 외 그룹)


**02_decode_reads.pl**

- 샘플명: fastp JSON 이 전혀 없는 --skip-fastp 경로에서 `<SAMPLE>_merged.fq(.gz)`/`<SAMPLE>.fpmerged.fq(.gz)` 의 샘플명(=raw/scaled_counts_matrix 컬럼명, decoded_reads_<sample>.tsv 파일명)이 `<SAMPLE>_merged` → `<SAMPLE>` 로 바뀜(README 의 '<SAMPLE> 이 컬럼명' 규약과 일치). JSON 이 있는 orchestrator 경로는 변화 없음.
- decoding_summary.tsv / undecoded_reads_*.tsv: revcomp 방향 리드의 실패 사유가 항상 no_cp 였던 것이 실제 결손 요소(no_op/no_hp/codon_fail)로 분류됨(cp_found/op_found/hp_found 플래그도 그 방향 기준). decoded 카운트·행렬은 불변.
- 종료 코드: 잘린/손상 .gz(gzip -cd 비정상 종료) → die(exit≠0), .decode_done 미기록(이전: exit 0 + 마커). 인자 전무 실행 → usage + exit 2(이전 0); -h 는 0 유지. `-p` 만 주면 이제 usage 대신 기본 경로로 실행 시도.
- BB 로딩: 편집거리 ≤2 앵커 쌍이 있어도 die 하지 않음(모호 변이체 제거·WARN); 같은 lib/type 의 perfect 완전 중복만 die(메시지 형식 변경: 'Duplicate <type> anchor sequence ... (line N, tag ...)'). CP 행 0 → die, 소문자/공백 type 허용, 무음 폐기 행이 WARN 으로 노출. 편집거리 2 케이스에서 이전 패치 트리는 첫 앵커의 변이체로 디코딩했으나 이제 해당 리드는 undecoded(no_hp 등)로 집계됨.
- length_window_stats.tsv: -W 1 로 mode 가 3/4/both 가 되면 active_centers 가 빈 값(이전: Top-K 피크 나열). -U 1 에서 두 번째 이후 샘플의 tol 이 CLI 기본값에서 다시 시작(이전: 앞 샘플 값 상속). 동률 피크 순서가 결정적(길이 오름차순)으로 고정.
- qc_checks.tsv: codon_len_not9/cp_len_out_of_27_29 컬럼명은 유지하되 기준이 -c/-Q canonical 길이(±1)로 변경(기본값 9/28 에서는 불변).
- 로그: 'Trie build complete' 줄에 anchors=/variants_dropped= 필드 추가; eval 내부 예외는 [DIE] 로 기록되지 않음; 'No fastp JSON' 메시지는 WARN→INFO 1회.

**03_call_hits.py (1–1140)**

- streaming_group_counts (--streaming_agg 1): lib_id/id 가 'NA','None','null' 등 pandas 기본 NA 토큰인 행이 더 이상 탈락하지 않고 비스트리밍 경로와 동일하게 집계됨(정상 파이프라인에서는 02 가 이런 lib_id 를 금지하므로 실질 변화 없음). 반환 프레임 컬럼 순서가 [LIB_ID, ID, DEL2, count_cols...] 로 고정됨(기존 dict 삽입 순서와 동일).
- load_bbinfo_auto: SMILES 셀이 비어 있거나 NA/nan/None/null 이면 '' 로 정규화 → 05_hybrid_annot.tsv 의 bb*_smiles 는 'nan' 대신 'NA', BB_SMILES_CONCAT 에서 해당 BB 생략(기존 'CCO.nan.C' → 'CCO.C'). 하류 05/06/07 의 SMILES 파싱 실패·클러스터 키 오염이 사라짐.
- BB 메타 헤더 판정이 토큰 완전일치로 바뀜: 첫 행의 어느 필드도 알려진 헤더 단어와 정확히 같지 않으면(예: 'Sequence_ID' 같은 변형 헤더) 헤더 없음으로 간주될 수 있음. 01_preprocess 가 쓰는 표준 헤더(type/seq/bb_id_fixed/cycle/tag_id/lib_id/smiles)는 그대로 인식.
- compute_synthon_scores: 3/4-cycle 혼합 실행에서 3-cycle 화합물의 SynthonScore 가 BB4='NA' 그룹 z 만큼 달라짐(단일 사이클 수 라이브러리에서는 수치 동일).
- --device cuda 가 CUDA 미가용 환경에서 오류 대신 경고 후 cpu 로 실행됨. 이때 --dtype auto 는 이미 float32 로 결정된 뒤이므로(main, B 범위) CPU 에서 float32 로 적합됨 — 정밀도가 필요하면 --dtype float64 명시.
- --r3_cols 기본값 None → []; hit_params_used.json 에는 r3_cols 가 포함되지 않아 오케스트레이터 캐시 해시 영향 없음.

**03_call_hits.py (main)**

- 03_glm_results.csv: glm_mode=top/skip 에서 비적합 행의 alpha_for_penalty 가 0.0 대신 빈칸(NaN) 으로 기록됨. 점수 계산(05_hybrid_annot.tsv 의 alpha_for_penalty/Penalty/HitScore_GLM)은 glm_mode=top 에서 적합 행 중앙값으로 대치된 값을 사용(skip/full 은 기존과 동일).
- --pseudocount_k 0 또는 음수는 이제 즉시 [ERROR] 종료(이전: DEL2=0 행에서 예외로 중단).
- auto_tune=0 + neg_strict=1 에서 --auto_lock w_neg / fix_neg_poisson 이 이제 존중됨(이전: neg_strict 가 lock 을 무시).
- --auto_syn_target 이 실제로 rho_syn 자동튠에 반영됨(기본 0.20 은 기존 동작과 동일).
- 로그 문구: '[AUTO] auto_tune=0 → using CLI values:' 접두 통일, '[NEG] neg_strict=1 → ...' 는 실제 적용 항목만 나열, hard_filter/preset 충돌 [WARN] 신설, '[GLM] alpha_for_penalty imputed ...' 신설.

**run_delegance_pipeline.py + .sh**

- run_delegance_pipeline.py: 02_decoded/decode_params.json 스키마 version 1→2 — `preprocess_hash` 제거, `merged_fingerprint`(01_fastp_out의 merged FASTQ path/size/mtime 목록) 추가. 기존 v1 캐시는 해시 불일치로 decode가 1회 재실행됨(--only hit에는 영향 없음).
- run_delegance_pipeline.py: 01_fastp_out/preprocess_params.json — --skip-fastp 실행에서 `fastq_fingerprint`가 항상 [] 이고 `merged_fingerprint`가 추가됨(skip-fastp 캐시 1회 재실행; fastp를 실제 실행한 캐시는 그대로 유효). skip-fastp이면 --fastq-dir 없이 실행 가능하며 01에 `-f`를 넘기지 않음(01은 $Bin/00_original_files 기본값을 갖지만 읽지 않음).
- run_delegance_pipeline.py: raw FASTQ 디렉터리가 없어도 출력+완료 마커가 있으면(--force-preprocess 아님) preprocess를 경고 후 건너뜀(이전: exit 2). 단계 시작 시 params.json 캐시가 삭제되므로 실패한 재실행 뒤에는 항상 재실행됨.
- run_delegance_pipeline.py: 레거시(마커 없음) 로그는 마지막 'Starting preprocess.'(01) / 'Merged dir = ' 또는 'Effective config:'(02) 이후에 'All done.'이 있어야 완료로 인정. 실패한 재실행이 있는 구 run은 이제 미완료로 판정되어 재실행됨(의도된 변경).
- run_delegance_pipeline.py: `--neg`에 3개 이상 컬럼 → exit 2(이전: 3번째 이후 무음 탈락). `--only hit`은 --skip-fastp 여부와 무관하게 merged FASTQ를 요구하지 않음. fastp 부재+merged 부재 신규 run의 오류 메시지 변경. `--dry-run`은 디렉터리/index.html을 만들지 않음. auto-opt는 --only all/hit에서만 동작하며 glm top+CUDA 조합은 device=cpu로 조정됨(디렉터리 토큰 dev_cpu_fp64).
- run_delegance_pipeline.py: index.html — 프리셋 토큰 없는 디렉터리 제목이 'Preset: glm (…)' → 'Preset: (none) (<dir>)'; 중복 `<p>Interactive: …</p>` 링크 제거.
- run_autopilot.sh: FASTQ_DIR이 비어 있으면(기본 후보 디렉터리 없음) 오류 대신 --fastq-dir을 생략하고 진행. 깨진 심볼릭 링크는 교체됨.
- postprocess_after_hits.sh: DEL2_COL 미설정 시 --del2_col을 넘기지 않음(exporter가 hit_params.json normalized_columns.del2 → 이름 휴리스틱 순으로 결정). 이전 기본값 'DEL2'가 필요하면 DEL2_COL=DEL2를 명시. pgrep 대기 패턴이 심볼릭 링크를 해소하지 않은 abspath 기준으로 바뀜.
- run_hits_then_postprocess.sh: 값 누락 옵션은 usage + exit 2.

**run_full_validation.py / verify_random_reads.py**

- run_full_validation.py JSON: undecoded_sample_validation 및 undecoded_sample_totals에 `len_conf_missing` 키 추가(샘플당 1; totals는 length_window_stats가 없는 샘플 수). 해당 샘플에서는 `reason_len_mismatch`가 더 이상 집계되지 않음(이전: len_out_of_range 행 전부가 mismatch로 계산). TSV 컬럼 `undecoded_reason_len_mismatch`는 그 경우 0.
- run_full_validation.py validation 카운터: 결손 cycle이 있는 손상 행의 gap 라벨이 실제 cycle 번호(`gap_C1_C3` 등)로 바뀔 수 있음. 정상 행의 라벨은 불변.
- run_full_validation.py: BB 표의 Codon cycle 값이 1..4 밖이면 무시(디코더와 동일). 그런 표에서 lib_expected_cycles가 4→3으로 바뀌어 cycles_mismatch 집계가 달라질 수 있음(디코더와 일치하는 방향).
- run_full_validation.py: 존재하지 않는 run 디렉터리는 exit 1([ERROR])로 종료(이전: FileNotFoundError traceback).
- run_full_validation.py: --max-cp-cands/--max-anchor-cands > 0일 때 동점 후보 절단 순서가 디코더와 동일하게 고정(기본값 0에서는 영향 없음).
- verify_random_reads.py: decoded 파일/raw_counts_matrix.tsv/--sample 컬럼 부재 시 [ERROR] 메시지와 exit 1(이전: FileNotFoundError / pandas ValueError traceback).

**04/05**

- 04: main_source CDS no longer has IMG1..IMG4 columns; new CDS img_source (columns smiles, img) shared by the detail callbacks. Rendered text/images unchanged.
- 04: --img_http_cap applies per unique SMILES instead of per unique BB.
- 04: --only_passed with no pass_filters column exits 1 (previously silently used all rows).
- 04: default --debug 0 no longer emits browser console logs; --debug 1 shows the last log line in the debug panel.
- 05: new flag --force-cluster; clustering is skipped with a WARN when top-N rows > 20,000 (no top<N>_clusters.tsv; cluster_* columns NA/0).
- 05: --preset that does not exist exits 1 (previously silent fallback to the newest preset).
- 05: zero rows with a numeric score exits 1 (previously empty TSV/HTML with exit 0).
- 05: bare --out-prefix NAME is written as <hybrid_dir>/NAME_* (previously cwd).
- 05: cluster_medoid (=0) column is emitted even when clustering is skipped/unavailable (fixed schema).
- 05: *_CPM/raw count columns in top<N>_hits.tsv are appended after the other columns (same set, different order); if the second-pass alignment check fails they are omitted with a WARN.
- 05: 'N/A' is treated as a missing BB in bb_frequency aggregation.
- README wording changes are proposed in fix2_Reports_README.txt (not applied).

**06_compare_top_hits.py**

- --bbavg-missing zero: denominator = positions where at least one of the two compounds has a fingerprint (was: number of BB columns); 3-cycle libraries can now reach similarity 1.0 (was capped at 0.75). Default skip unchanged.
- Derived labels (no --labels): 05_hybrid_annot.tsv paths are labelled <run>_<preset> instead of 05_hybrid_annot.tsv(_N); duplicate suffixes may skip numbers (_2,_3) to stay unique. Directory inputs unchanged.
- *_diverse.tsv (tier_report and per-run): when clustering is not performed (--cluster 0, --include-summary 0, RDKit/bb*_smiles/valid fingerprints absent) the file is written with header only plus an [INFO] line; per-run diverse file is now always created.
- tier_report inactive_score/both_score: -inf -> NaN (same as active_score).
- RDKit/smiles-absent fallback outputs gain a cluster_medoid (=0) column.
- --roles: a role assigned to more than one run is now an error (was: second run silently ignored); missing active+inactive pair prints a warning.
- --score-col not present in a run header -> explicit [ERROR] exit (was: pandas usecols exception).
- HTML: Download CSV omits index/structure_html/*_img columns; specificity tables no longer show duplicate group/group_code columns.
- bokeh missing -> [WARN] and TSV-only run (was: NameError crash).

**07_tiered_report.py**

- --cluster 0: tier_report_*_diverse.tsv가 그룹 전체 목록 → 빈 표(헤더만), HTML Diverse 표 생략, summary n_diverse='n/a (clustering disabled)'. RDKit 부재 시 동작(빈 diverse)과 final_hits 선택 결과는 불변.
- --labels 생략 + 05_hybrid_annot.tsv 파일 경로 입력: 기본 라벨이 '05_hybrid_annot.tsv(_2)' → run 루트 디렉터리명(또는 부모 디렉터리명) → prefixed 컬럼명/그룹 표시명 변경. 디렉터리 입력은 불변.
- HTML Final hits 탭(group_code가 여러 값인 표): 검색/페이지 이동/CSV 다운로드 정렬이 group_code 순서 → 선택 rank 컬럼 순(이전: rank 단독으로 그룹 인터리브). TSV/XLSX 불변.
- 동점 score의 compound_key 중복(다중 LIB) 행 채택이 ID 오름차순으로 고정 → 해당 행 샘플 카운트/NEG 필터 결과가 이전 실행과 달라질 수 있음(이전엔 비결정적).
- bokeh 3.x에서 DataTable 셀 CSS(.struct-popup hover, .group-badge 색)가 실제 적용됨(시각 변화만).
- 내부 시그니처: _make_table_controls(group_sort_col=), _tie_break_arrays(exclude_prefixes=), _cluster_candidates(sample_prefixes=) 선택 인자 추가(기본값으로 기존 호출 호환).

### 2차 회귀 검증 요약 (gpu-4080, 2026-09-06)

- 기존 회귀(정적·GLM 단위·합성 03 8구성·FASTQ 엔드투엔드·리포트/익스포터 12종): FAIL 0.
- 3버전 비교(동일 합성 6,000-태그 행렬, numpy): alpha_for_penalty 중앙값 pristine 0.178 / 1차 0.043 / 2차 0.150; 태그별 alpha 비율(2차/pristine) 중앙값 1.000 (1차 0.325) — delta 클램프가 관측치 pseudocount의 과산포 축소 편향을 제거. HitScore_GLM Spearman vs pristine: 1차 0.9505 → 2차 0.9937. 0-카운트 LFC floor −1.58 유지. 진짜 hit 30개 중 GLM 29·Consensus 27(1차 30/27, pristine 29/28).
- 신규 케이스: 편집거리 1 앵커 쌍(die 없이 변이체 12개 제거 후 정상 디코딩), 소문자 type(94.35 % 디코딩), 잘린 .gz(exit≠0, 마커 미기록, 원인 메시지 decompressor), 실패 재실행 후 resume(decode 재실행), subsample 자동 페어링(bcl2fastq·`_1/_2` 모두 정상). bcl2fastq 이름 `R1C1_S1_L001_R1_001`은 페어링·디코딩 정상이며 샘플명이 `R1C1_S1_L001`로 남는 것은 문서화된 동작(README 참조).

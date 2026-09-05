#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
04_build_interactive_report.py

DELeGANce Master (master_annotated.tsv) → Interactive HTML (Bokeh, Top-N by HitScore)

수정사항:
- 클릭 이벤트 핸들링 개선: 단순화된 selection 기반 방식 사용
- JavaScript 코드 정리 및 디버깅 개선
- 이벤트 핸들러 충돌 해결
"""

import os
import re
import math
import time
import argparse
import base64
from io import BytesIO
from typing import Dict, List, Tuple, Set, Optional

import pandas as pd
import numpy as np

# -------- RDKit (optional) --------
_HAS_RDKIT = True
_RDKIT_IMPORT_ERR = ""
try:
    from rdkit import Chem
    from rdkit.Chem import Draw
except Exception as _e:
    _HAS_RDKIT = False
    _RDKIT_IMPORT_ERR = str(_e)

# -------- Bokeh --------
from bokeh.plotting import figure, output_file
from bokeh.io import save  # write HTML only; never open a browser (batch-safe)
from bokeh.models import (
    ColumnDataSource, CDSView, BooleanFilter, HoverTool, TapTool,
    Select, Div, DataTable, TableColumn, NumberFormatter, CustomJS,
    ColorBar, LinearColorMapper, CustomAction, Spacer, Button
)
from bokeh.palettes import Blues9
from bokeh.layouts import gridplot
from bokeh.layouts import row as bk_row, column as bk_column
from bokeh.events import Tap, DocumentReady

# Hard cap for performance
TOP_CAP = 10000


# ------------------------- CLI -------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="DELeGANce master_annotated.tsv → Interactive HTML (Top-N by HitScore)"
    )
    p.add_argument("--master_tsv", dest="master_tsv",
                   default="hit_results/master_annotated.tsv",
                   help="master_annotated.tsv 경로(압축 .gz 가능)")
    # alias
    p.add_argument("--master", dest="master_tsv", help=argparse.SUPPRESS)

    p.add_argument("--top_hitscore", type=int, default=10000,
                   help="HitScore 기준 상위 N만 시각화 (성능상 최대 10000 강제 적용)")
    p.add_argument("--only_passed", action="store_true",
                   help="pass_filters==True 인 레코드만 사용")
    p.add_argument("--bbinfo",
                   default="DELeGANce_out/BB_information_fixed.tsv",
                   help="(선택) master에 SMILES/LibID 비어있는 경우 보완용 BB 정보(.gz 가능)")
    p.add_argument("--out", default="bokeh_master_topN.html", help="출력 HTML 파일명")
    p.add_argument("--top_table", type=int, default=500, help="오른쪽 Top 테이블 행 수(기본 500)")
    p.add_argument("--imgsize", type=int, default=150, help="RDKit 분자 이미지 한 변 픽셀(기본 150)")
    p.add_argument("--debug", type=int, choices=[0,1], default=0,
                   help="1이면 디버그 패널/콘솔 로깅 활성화")
    # HTTP fallback for images when RDKit is unavailable
    p.add_argument("--img_fallback", choices=["auto","none","http"], default="none",
                   help="RDKit 미설치 시 이미지 대체 방식(기본 none→이미지 생략). "
                        "http/auto 선택 시 BB SMILES가 외부 서버(cactus.nci.nih.gov)로 전송됨 — 명시적 opt-in 필요")
    p.add_argument("--img_http_cap", type=int, default=400,
                   help="HTTP 대체 이미지 최대 개수(유니크 BB 기준, 기본 400)")
    p.add_argument("--min_dot", type=float, default=4, help="포인트 최소 크기")
    p.add_argument("--max_dot", type=float, default=16, help="포인트 최대 크기")
    p.add_argument("--plot_height", type=int, default=220, help="BB x BB 산점도 각 패널 높이(px)")
    # 사용자 지정 메트릭 목록(공백/콤마 구분). 지정 안하면 자동 탐지
    p.add_argument("--metrics", type=str, default="",
                   help="색/크기 기준으로 선택할 수치 메트릭들(콤마/공백 구분). 예: HitScore,mean_log2FC_BEAD,log2Boost_R2vsR1")
    return p.parse_args()


# ------------------------- 경로 유틸 -------------------------
def _is_file(p: str) -> bool:
    return bool(p) and os.path.isfile(p)

def _with_gz(p: str) -> Optional[str]:
    if not p:
        return None
    gz = p + ".gz"
    return gz if os.path.isfile(gz) else None

def resolve_path(primary: str, candidates: List[str]) -> str:
    """Return the first on-disk path among primary/candidate options, accepting .gz variants."""
    pool: List[str] = []
    if primary:
        pool.append(primary)
    pool.extend([c for c in candidates if c])
    for p in pool:
        if _is_file(p):
            return p
        gz = _with_gz(p)
        if gz:
            return gz
    # Nothing found: report the path the user actually asked for (not the last fallback candidate)
    return primary if primary else (pool[-1] if pool else "")


# ------------------------- 공용 유틸 -------------------------
# Missing-value markers: 03_call_hits writes NA as the literal "NA" (na_rep="NA"); pandas may
# also read them as NaN which becomes "nan" after astype(str).
_MISSING_MARKERS = {"", "NA", "N/A", "nan", "None"}

def _is_missing(s: pd.Series) -> pd.Series:
    """Boolean mask of missing cells (NaN or any textual NA marker)."""
    return s.isna() | s.astype(str).str.strip().isin(_MISSING_MARKERS)

def _clean_str(v: object) -> str:
    """Return '' for NaN/None/NA-marker values, else the stripped string."""
    if v is None:
        return ""
    try:
        if pd.isna(v):
            return ""
    except (TypeError, ValueError):
        pass
    s = str(v).strip()
    return "" if s in _MISSING_MARKERS else s

def pick_col(df: pd.DataFrame, *candidates: str) -> Optional[str]:
    """Return the first candidate column present in df (e.g. pick_col(df, 'ID_x', 'ID')).

    Backward/forward compatible accessor for hybrid tables whose merge suffixes (_x/_y)
    may be removed in a future 03_call_hits release.
    """
    for c in candidates:
        if c and c in df.columns:
            return c
    return None


def sanitize_field_name(s: str) -> str:
    s2 = re.sub(r"\W+", "_", s)
    if re.match(r"^\d", s2):
        s2 = "_" + s2
    return s2

def _strip_lib_suffix(bb: str) -> str:
    if bb is None:
        return "NA"
    s = str(bb)
    if s in ("", "NA", "nan", "None"):
        return "NA"
    return re.sub(r"_LIB[\w\.-]+$", "", s)

def _strip_lib_anywhere(s: str) -> str:
    if s is None:
        return ""
    return re.sub(r"_LIB[\w\.-]+", "", str(s))

def _strip_series(s: pd.Series) -> pd.Series:
    return s.astype(str).str.replace("\u00A0", " ", regex=False).str.strip()

def smiles_to_base64(smiles: str, size=(150, 150)) -> str:
    if not isinstance(smiles, str) or not smiles.strip():
        return "No SMILES"
    if not _HAS_RDKIT:
        return f"RDKit not available ({_RDKIT_IMPORT_ERR})"
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return "Invalid SMILES"
        img = Draw.MolToImage(mol, size=size)
        buf = BytesIO()
        img.save(buf, format="PNG")
        b64 = base64.b64encode(buf.getvalue()).decode("utf-8")
        return f"data:image/png;base64,{b64}"
    except Exception as e:
        return f"Image Gen Error: {e}"

def parse_id_to_fields(id_str: str) -> Tuple[int, List[str]]:
    """
    기대형식: "<cycles>_<BB1>_<BB2>_<BB3>_<BB4>"
    - 구분자 _, |, :, ;, , , 공백 혼용 허용
    - 3-cycle이면 BB4='NA'
    """
    s = str(id_str).strip()
    if not s:
        return 3, ["NA", "NA", "NA", "NA"]
    m = re.match(r"^\s*(\d+)\s*[\|_,:;/\s]+(.+?)\s*$", s)
    if not m:
        toks = re.split(r"[\|_,:;/\s]+", s)
        try:
            cyc = int(toks[0]); parts = toks[1:]
        except Exception:
            cyc = 3; parts = toks
    else:
        cyc = int(m.group(1))
        parts = re.split(r"[\|_,:;/\s]+", m.group(2))
    parts = [p for p in parts if p != ""]

    # Stitch namespaced BB IDs like XBA0038_LIBDEL004 that get split by '_' tokenisation.
    merged: List[str] = []
    i = 0
    while i < len(parts) and len(merged) < 4:
        t = parts[i]
        if (
            i + 1 < len(parts)
            and parts[i + 1].startswith("LIB")
            and (t not in ("NA", "") and not t.startswith("LIB"))
        ):
            merged.append(f"{t}_{parts[i + 1]}")
            i += 2
        else:
            merged.append(t)
            i += 1
    while len(merged) < 4:
        merged.append("NA")
    return (3 if cyc not in (3, 4) else cyc), merged[:4]


def _normalize_bb_value(value: object) -> str:
    s = "" if value is None else str(value).strip()
    if s in ("", "NA", "N/A", "nan", "None"):
        return "NA"
    return s


def sort_bb_categories(values: List[str]) -> List[str]:
    uniq = []
    seen = set()
    for v in values:
        s = _normalize_bb_value(v)
        if s not in seen:
            uniq.append(s)
            seen.add(s)

    def _key(s: str) -> Tuple[int, int, str, int, str]:
        if s == "NA":
            return (1, 1, "", 1_000_000_000, "")
        m = re.match(r"^([A-Za-z]+)(\d+)(.*)$", s)
        if m:
            return (0, 0, m.group(1), int(m.group(2)), m.group(3))
        return (0, 1, s, 1_000_000_000, "")

    return sorted(uniq, key=_key)


# ------------------------- bbinfo 로더 (SMILES & lib_id 동시 추출) -------------------------
def _pick_by_names(df: pd.DataFrame, name_candidates: List[str]) -> Optional[str]:
    lower_map = {c.lower(): c for c in df.columns}
    for cand in name_candidates:
        c = lower_map.get(cand.lower())
        if c is not None:
            return c
    return None

def _guess_smiles_col(df: pd.DataFrame) -> Optional[str]:
    by_name = _pick_by_names(df, ["smiles", "SMILES", "bb_smiles"])
    if by_name:
        return by_name
    pat_allowed = re.compile(r"^[A-Za-z0-9@+\-\[\]\(\)=#\\/\.]+$")
    def smiles_like_score(series: pd.Series) -> int:
        cnt = 0
        for v in series.astype(str).values:
            if not v or v.strip() == "":
                continue
            if not pat_allowed.match(v):
                continue
            if any(ch in v for ch in ("=", "#", "[", "]", "(", ")")):
                cnt += 1
        return cnt
    best = (0, None)
    for col in df.columns:
        s = _strip_series(df[col])
        score = smiles_like_score(s)
        if score > best[0]:
            best = (score, col)
    if best[0] >= max(5, int(0.05 * len(df))):
        return best[1]
    return None

def _guess_bb_col_by_intersection(df: pd.DataFrame, needed_bbs: Set[str]) -> Optional[str]:
    best_match_count = -1
    best_unique = -1
    best_col = None
    exclude = set()
    for c in df.columns:
        cl = c.lower()
        if "smile" in cl or "seq" in cl or "sequence" in cl or "desc" in cl:
            exclude.add(c)
    for col in df.columns:
        if col in exclude:
            continue
        s = _strip_series(df[col])
        matches = s.isin(needed_bbs)
        match_count = int(matches.sum())
        unique_match = int(s[matches].nunique())
        if (match_count > best_match_count) or (match_count == best_match_count and unique_match > best_unique):
            best_match_count = match_count
            best_unique = unique_match
            best_col = col
    if best_match_count <= 0 or best_col is None:
        return None
    return best_col

def load_bbinfo(path_bbinfo: str, needed_bbs: Set[str]) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    BB 파일 로더:
    - 00_7DEL_BB_information_20241126.txt 형식이면 type=='Codon'만 사용
    - bb_id 열과 smiles, lib_id 열을 자동 탐지
    - 반환: (bb_to_smiles, bb_to_lib)
    """
    p = path_bbinfo
    if not _is_file(p):
        gz = _with_gz(p)
        if gz:
            p = gz
    if not _is_file(p):
        raise FileNotFoundError(
            f"[bbinfo] 파일을 찾을 수 없습니다: {path_bbinfo}\n"
            f"  해결: --bbinfo /절대/또는/상대/경로.tsv[.gz] 로 지정해 주세요.\n"
            f"  (예) --bbinfo DELeGANce_out/BB_information_fixed.tsv  또는\n"
            f"       --bbinfo 00_7DEL_BB_information_20241126.txt"
        )
    df = pd.read_csv(p, sep="\t", dtype=str, compression="infer")
    for c in df.columns:
        df[c] = _strip_series(df[c])
    # Codon만
    if "type" in df.columns:
        df_type = df["type"].astype(str).str.lower()
        if any(df_type.str.contains("codon", regex=False)):
            df = df[df_type.str.contains("codon", regex=False)].copy()

    col_id = _pick_by_names(df, ["bb_id_fixed", "bb_id", "id"])
    col_smiles = _pick_by_names(df, ["smiles", "SMILES", "bb_smiles"])
    col_lib = _pick_by_names(df, ["lib_id", "library", "lib"])

    # If both bb_id_fixed and bb_id exist, prefer the one that best matches needed_bbs
    try:
        guessed = _guess_bb_col_by_intersection(df, needed_bbs)
        if guessed is not None and col_id is not None and guessed != col_id:
            # Only switch if guessed provides strictly more matches
            def _match_count(c):
                try:
                    return int(_strip_series(df[c]).isin(needed_bbs).sum())
                except Exception:
                    return -1
            if _match_count(guessed) > _match_count(col_id):
                col_id = guessed
    except Exception:
        pass

    if col_id is None:
        col_id = _guess_bb_col_by_intersection(df, needed_bbs)
    if col_smiles is None:
        col_smiles = _guess_smiles_col(df)

    if col_id is None:
        raise ValueError("[bbinfo] BB 식별 컬럼을 찾지 못했습니다.")
    # col_smiles 없을 수도 있음
    # col_lib 없으면 모두 빈 문자열

    sub = df[df[col_id].isin(needed_bbs)].copy()
    bb_to_smiles: Dict[str, str] = {}
    bb_to_lib: Dict[str, str] = {}

    for _, r in sub.iterrows():
        bb = str(r[col_id])
        # NOTE: _strip_series() above turned NaN into the string "nan"; pd.isna() would never fire.
        smi = "" if col_smiles is None else _clean_str(r[col_smiles])
        lib = "" if col_lib is None else _clean_str(r[col_lib])
        bb_to_smiles.setdefault(bb, smi)
        bb_to_lib.setdefault(bb, lib)

        # If we used a fixed id like 'XBA0040.1', also map its base 'XBA0040'
        if "." in bb:
            base = bb.split(".", 1)[0]
            if base in needed_bbs:
                if base not in bb_to_smiles or not bb_to_smiles[base]:
                    bb_to_smiles[base] = smi
                if base not in bb_to_lib or not bb_to_lib[base]:
                    bb_to_lib[base] = lib

    # 누락 보정
    for bb in needed_bbs:
        bb_to_smiles.setdefault(bb, "")
        bb_to_lib.setdefault(bb, "")

    return bb_to_smiles, bb_to_lib


# ------------------------- Master 로딩 & 준비 -------------------------
def load_master(path_master: str) -> pd.DataFrame:
    if not _is_file(path_master):
        gz = _with_gz(path_master)
        if gz:
            path_master = gz
    if not _is_file(path_master):
        raise FileNotFoundError(
            f"[master] 파일을 찾을 수 없습니다: {path_master}\n"
            f"  해결: --master_tsv /절대/또는/상대/경로.tsv[.gz] 로 지정해 주세요."
        )
    # 문자열로 유지해야 하는 열 (ID/LIB/BB/CP 계열은 숫자형으로 보여도 float 변환 금지)
    keep_str = {
        "id", "ID", "ID_x", "ID_y", "id_x", "id_y",
        "primary_lib_id_used", "LIB_ID", "LIB_ID_x", "LIB_ID_y", "lib_id", "lib_id_x", "lib_id_y",
        "BB1", "BB2", "BB3", "BB4", "BB1_x", "BB2_x", "BB3_x", "BB4_x", "CP", "CP_x",
        "bb1_id","bb2_id","bb3_id","bb4_id",
        "bb1_smiles","bb2_smiles","bb3_smiles","bb4_smiles",
        "fail_reasons", "pass_filters",
    }
    # Memory: only force str on the columns we must keep textual; let pandas infer the rest
    # (hybrid tables have ~100 columns; dtype=str on all of them costs GBs at 10^6-10^7 rows).
    header = pd.read_csv(path_master, sep="\t", nrows=0, compression="infer").columns.tolist()
    dtype_map = {c: str for c in header if c in keep_str}
    df = pd.read_csv(path_master, sep="\t", dtype=dtype_map, compression="infer", low_memory=False)
    # Preserve the previous all-string semantics for boolean-like columns (e.g. GLM_hit "True"/"False")
    for c in df.columns:
        if c in keep_str:
            continue
        if pd.api.types.is_bool_dtype(df[c]):
            df[c] = df[c].astype(str)
        elif df[c].dtype == object:
            df[c] = df[c].where(df[c].isna(), df[c].astype(str))

    # Normalize common column variants produced by 03_call_hits (legacy names still supported)
    if "id" not in df.columns:
        k = pick_col(df, "ID_x", "ID", "ID_y", "id_x", "id_y")
        if k is not None:
            df["id"] = df[k].astype(str)
    # HitScore: prefer HitScore, else fallback to HitScore_GLM or HitScore_RS
    if "HitScore" not in df.columns:
        k = pick_col(df, "HitScore_GLM", "HitScore_RS")
        if k is not None:
            df["HitScore"] = df[k]
    if "id" not in df.columns or "HitScore" not in df.columns:
        raise ValueError("[master] 'id'와 'HitScore' 컬럼이 필요합니다. (또는 'ID'/'ID_x'와 'HitScore_GLM'/'HitScore_RS')")

    # 선택적 수치 변환(해당 열의 90% 이상이 숫자로 변환될 수 있는 경우에만) — 문자열(object) 열만 대상
    for c in df.columns:
        if c in keep_str or c == "id":
            continue
        s = df[c]
        if s.dtype != object:
            continue  # already numeric by inference
        try:
            conv = pd.to_numeric(s, errors="coerce")
        except Exception:
            continue
        valid = conv.notna().sum()
        if valid >= 0.9 * len(conv):
            df[c] = conv
        # else: 그대로 문자열 유지
    return df


def build_records_from_master(df_master: pd.DataFrame,
                              top_hitscore: int,
                              only_passed: bool) -> pd.DataFrame:
    df = df_master.copy()

    # pass_filters 옵션
    if only_passed and "pass_filters" in df.columns:
        if df["pass_filters"].dtype != bool:
            df["pass_filters"] = df["pass_filters"].astype(str).str.lower().isin(["true", "1", "t", "yes", "y"])
        df = df[df["pass_filters"]]

    # HitScore 정렬 & 상위 N (deterministic: stable sort + secondary keys, matching 03's top-K ordering)
    df["HitScore"] = pd.to_numeric(df["HitScore"], errors="coerce")
    df = df[np.isfinite(df["HitScore"])].copy()
    sort_keys = ["HitScore"]
    sort_asc = [False]
    if "HitScore_RS" in df.columns and pd.api.types.is_numeric_dtype(df["HitScore_RS"]):
        sort_keys.append("HitScore_RS"); sort_asc.append(False)
    sort_keys.append("id"); sort_asc.append(True)
    df = df.sort_values(sort_keys, ascending=sort_asc, kind="mergesort")
    if top_hitscore is not None and top_hitscore > 0:
        df = df.head(top_hitscore).copy()

    # cycles & BB1..BB4 — accept bb1_id.. (03 annotation), BB1.. or BB1_x.. (hybrid merge suffix)
    bb_src = {}
    for i in (1, 2, 3, 4):
        bb_src[i] = pick_col(df, f"bb{i}_id", f"BB{i}", f"BB{i}_x")
    has_bb_cols = all(bb_src[i] is not None for i in (1, 2, 3, 4))
    if has_bb_cols:
        for i in (1, 2, 3, 4):
            s = df[bb_src[i]]
            s = s.where(~_is_missing(s), "NA").astype(str)
            df[f"BB{i}"] = s
        cyc = []
        for s, b4 in zip(df["id"].astype(str).tolist(), df["BB4"].astype(str).tolist()):
            m = re.match(r"^\s*(\d+)[_\|:,;/\s]", s)
            if m:
                c = int(m.group(1))
                c = c if c in (3,4) else (3 if b4 in ("", "NA") else 4)
            else:
                c = 3 if b4 in ("", "NA") else 4
            cyc.append(c)
        df["cycles"] = cyc
    else:
        # id에서 파싱
        BB1,BB2,BB3,BB4,cycles = [],[],[],[],[]
        for s in df["id"].astype(str):
            c, parts = parse_id_to_fields(s)
            cycles.append(c)
            b1,b2,b3,b4 = parts
            BB1.append(b1); BB2.append(b2); BB3.append(b3); BB4.append(b4)
        df["BB1"],df["BB2"],df["BB3"],df["BB4"],df["cycles"] = BB1,BB2,BB3,BB4,cycles

    # LibID (1) primary_lib_id_used 우선, (2) LIB_ID_x / LIB_ID (또는 _y, lib_id*) 보조
    # NOTE: Series.replace(dict) is NOT chained ("nan"->"" would never reach np.nan); mask explicitly.
    def _lib_series(s: pd.Series) -> pd.Series:
        s2 = s.astype(str).str.strip()
        return s2.where(~_is_missing(s), np.nan)
    lib_src = pick_col(df, "primary_lib_id_used", "LIB_ID_x", "LIB_ID", "LIB_ID_y", "lib_id", "lib_id_x", "lib_id_y")
    if lib_src is not None:
        df["LibID"] = _lib_series(df[lib_src])
    else:
        df["LibID"] = np.nan

    # 반환
    return df


# ------------------------- 메인 -------------------------
def main():
    args = parse_args()
    # Enforce a hard cap for speed
    if args.top_hitscore is None or args.top_hitscore <= 0 or args.top_hitscore > TOP_CAP:
        print(f"  - [info] top_hitscore={args.top_hitscore} → capped to {TOP_CAP} for performance")
        top_n = TOP_CAP
    else:
        top_n = min(args.top_hitscore, TOP_CAP)
    print("=" * 90)
    print("DELeGANce master_annotated.tsv → Interactive HTML (Top-N by HitScore)")
    print("=" * 90)

    # 0) 경로 자동 탐색
    resolved_master = resolve_path(
        args.master_tsv,
        candidates=[
            "hit_results/master_annotated.tsv",
            "./master_annotated.tsv",
        ],
    )
    resolved_bbinfo = resolve_path(
        args.bbinfo,
        candidates=[
            "DELeGANce_out/BB_information_fixed.tsv",
            "BB_information_fixed.tsv",
            "00_7DEL_BB_information_20241126.txt",
        ],
    )

    # 1) master 로딩
    print(f"[{time.strftime('%H:%M:%S')}] Load master: {resolved_master}")
    master = load_master(resolved_master)
    print(f"  - master 전체 행 수: {len(master)} (열 {len(master.columns)})")

    # 2) 상위 N HitScore 선택(+옵션 pass_filters)
    print(f"  - 옵션: only_passed={args.only_passed}, top_hitscore={top_n} (capped)")
    df = build_records_from_master(master, top_n, args.only_passed)
    print(f"  - 사용 행 수: {len(df)}")
    if df.empty:
        raise SystemExit("[master] 선택된 행이 0개입니다 (only_passed / top_hitscore / HitScore 결측 여부를 확인하세요).")

    # 3) LibID 보완 및 SMILES 보완 (bbinfo)
    # Missing SMILES arrive as "NA" (03 na_rep) or NaN, never as "" — use _is_missing on the selected rows.
    need_bb_map = False
    if df["LibID"].isna().any():
        need_bb_map = True
    for k in ("bb1_smiles","bb2_smiles","bb3_smiles","bb4_smiles"):
        if k in df.columns and _is_missing(df[k]).any():
            need_bb_map = True

    # 필요 BB 집합 수집
    needed_bbs: Set[str] = set()
    for bcol in ("BB1","BB2","BB3","BB4"):
        needed_bbs |= {b for b in df[bcol].astype(str).tolist() if b and b != "NA"}

    bb_to_smiles: Dict[str, str] = {}
    bb_to_lib: Dict[str, str] = {}
    if need_bb_map and len(needed_bbs) > 0:
        try:
            print(f"[{time.strftime('%H:%M:%S')}] Load BB info (SMILES & LibID): {resolved_bbinfo}")
            bb_to_smiles, bb_to_lib = load_bbinfo(resolved_bbinfo, needed_bbs)
        except Exception as e:
            print(f"  - [경고] bbinfo 로딩 실패: {e}")

    # LibID 보완: 첫 유효 BB의 lib_id로 채움
    if df["LibID"].isna().any() and bb_to_lib:
        def fill_libid(row):
            for k in ("BB1","BB2","BB3","BB4"):
                bb = str(row[k])
                if bb and bb != "NA":
                    li = bb_to_lib.get(bb, "")
                    if li:
                        return li
            return np.nan
        df["LibID"] = df["LibID"].fillna(df.apply(fill_libid, axis=1))

    # 최종 NA 채우기
    df["LibID"] = df["LibID"].fillna("NA")

    # 4) SMILES/이미지 구성
    # master에 SMILES 있으면 우선 사용, 없으면 보완
    def pick_smiles(i):
        # df keeps master's original index (subset), so align on df.index; NA markers -> ""
        if i in df.columns:
            s = df[i]
        elif i in master.columns:
            s = master[i].reindex(df.index)
        else:
            return pd.Series([""] * len(df), index=df.index)
        return s.where(~_is_missing(s), "").astype(str)
    df["SMILES1"] = pick_smiles("bb1_smiles")
    df["SMILES2"] = pick_smiles("bb2_smiles")
    df["SMILES3"] = pick_smiles("bb3_smiles")
    df["SMILES4"] = pick_smiles("bb4_smiles")

    if bb_to_smiles:
        for i, col in enumerate(("SMILES1","SMILES2","SMILES3","SMILES4"), start=1):
            mask_empty = _is_missing(df[col])
            if mask_empty.any():
                bbcol = f"BB{i}"
                df.loc[mask_empty, col] = df.loc[mask_empty, bbcol].map(lambda bb: bb_to_smiles.get(bb, ""))

    # RDKit 이미지
    # Fallback helper: HTTP fetch from NCI Cactus (no external deps)
    def http_smiles_to_base64(smiles: str, size: int) -> str:
        try:
            import urllib.request as urlreq, urllib.parse as urlparse
            if not isinstance(smiles, str) or not smiles.strip():
                return "No SMILES"
            enc = urlparse.quote(smiles, safe='')
            # NCI Cactus image endpoint
            url = f"https://cactus.nci.nih.gov/chemical/structure/{enc}/image?format=png&width={size}&height={size}"
            req = urlreq.Request(url, headers={"User-Agent":"DELeGANce/1.0"})
            with urlreq.urlopen(req, timeout=5) as resp:
                data = resp.read()
            b64 = base64.b64encode(data).decode("utf-8")
            return f"data:image/png;base64,{b64}"
        except Exception as e:
            return f"HTTP Image Error: {e}"

    if not _HAS_RDKIT and (args.img_fallback in ("auto","http")):
        print("[WARN] SMILES will be sent to cactus.nci.nih.gov (--img_fallback=%s). "
              "Proprietary building-block structures leave this machine; use --img_fallback none to disable."
              % args.img_fallback)
        print(f"  - [알림] RDKit 미설치: HTTP 대체 이미지 사용 시도({args.img_http_cap}개 제한). ({_RDKIT_IMPORT_ERR})")
    elif not _HAS_RDKIT:
        print(f"  - [알림] RDKit 미설치: 이미지 없이 진행합니다. ({_RDKIT_IMPORT_ERR})")
    unique_bb = set(df["BB1"]) | set(df["BB2"]) | set(df["BB3"]) | set(df["BB4"])
    unique_bb.discard("NA")
    bb_to_img: Dict[str, str] = {}
    if _HAS_RDKIT:
        print(f"[{time.strftime('%H:%M:%S')}] Generate RDKit Base64 images (unique BB set: {len(unique_bb)})...")
        ok, ng = 0, 0
        for bb in sorted(unique_bb):
            # 해당 BB의 첫 번째 SMILES를 찾아 그립니다
            smi = ""
            for i in (1,2,3,4):
                s = df.loc[df[f"BB{i}"] == bb, f"SMILES{i}"]
                if len(s) and isinstance(s.iloc[0], str) and s.iloc[0].strip():
                    smi = s.iloc[0]; break
            img = smiles_to_base64(smi, (args.imgsize, args.imgsize))
            bb_to_img[bb] = img
            if isinstance(img, str) and img.startswith("data:image"):
                ok += 1
            else:
                ng += 1
        print(f"  - RDKit 이미지: 성공 {ok}, 실패/결측 {ng}")
    elif args.img_fallback in ("auto","http"):
        cap = max(0, int(args.img_http_cap))
        sel = sorted(unique_bb)[:cap]
        print(f"[{time.strftime('%H:%M:%S')}] HTTP 이미지 페치 (유니크 BB {len(unique_bb)} → {len(sel)} 제한)...")
        ok, ng = 0, 0
        consecutive_http_err = 0
        for bb in sel:
            smi = ""
            for i in (1,2,3,4):
                s = df.loc[df[f"BB{i}"] == bb, f"SMILES{i}"]
                if len(s) and isinstance(s.iloc[0], str) and s.iloc[0].strip():
                    smi = s.iloc[0]; break
            img = http_smiles_to_base64(smi, args.imgsize)
            bb_to_img[bb] = img
            if isinstance(img, str) and img.startswith("data:image"):
                ok += 1
                consecutive_http_err = 0
            else:
                ng += 1
                if isinstance(img, str) and img.startswith("HTTP Image Error"):
                    consecutive_http_err += 1
                    # Offline/blocked network: stop early instead of paying cap x 5s timeouts
                    if consecutive_http_err >= 3:
                        print("  - [경고] HTTP 이미지 요청이 연속 실패하여 나머지 요청을 중단합니다 (네트워크 차단?).")
                        break
        print(f"  - HTTP 이미지: 성공 {ok}, 실패/결측 {ng} (캡 {cap})")

    def _img(bb): return bb_to_img.get(bb, "N/A") if (bb and bb != "NA") else "N/A"
    df["IMG1"] = df["BB1"].map(_img)
    df["IMG2"] = df["BB2"].map(_img)
    df["IMG3"] = df["BB3"].map(_img)
    df["IMG4"] = df["BB4"].map(_img)

    # 4.5) Normalize BB/ID for display (remove LIB suffix)
    for k in ("BB1", "BB2", "BB3", "BB4"):
        df[k] = df[k].astype(str).map(_strip_lib_suffix)
    def _make_display_id_row(r: pd.Series) -> str:
        lib = str(r.get("LibID", "")).strip()
        if lib and lib not in ("NA", "nan", "None"):
            return f"{lib}_{r.get('BB1','NA')}_{r.get('BB2','NA')}_{r.get('BB3','NA')}_{r.get('BB4','NA')}"
        return _strip_lib_anywhere(r.get("id", ""))
    df["id"] = df.apply(_make_display_id_row, axis=1)

    # 5) 시각화용 데이터 가공
    df["CombinedID"] = (
        df["LibID"].astype(str) + "_" +
        df["BB1"].astype(str) + "_" + df["BB2"].astype(str) + "_" + df["BB3"].astype(str) + "_" + df["BB4"].astype(str)
    )
    df["DisplayID"] = (
        df["LibID"].astype(str) + "_" +
        df["BB1"].astype(str).map(_strip_lib_suffix) + "_" +
        df["BB2"].astype(str).map(_strip_lib_suffix) + "_" +
        df["BB3"].astype(str).map(_strip_lib_suffix) + "_" +
        df["BB4"].astype(str).map(_strip_lib_suffix)
    )
    df["BB1_disp"] = df["BB1"].astype(str).map(_strip_lib_suffix)
    df["BB2_disp"] = df["BB2"].astype(str).map(_strip_lib_suffix)
    df["BB3_disp"] = df["BB3"].astype(str).map(_strip_lib_suffix)
    df["BB4_disp"] = df["BB4"].astype(str).map(_strip_lib_suffix)

    # 메트릭 후보(자동 탐지 or 사용자 지정)
    def _is_numeric_col(s: pd.Series) -> bool:
        return pd.api.types.is_numeric_dtype(s)

    if args.metrics.strip():
        metric_candidates = [m.strip() for m in re.split(r"[,\s]+", args.metrics.strip()) if m.strip()]
        metric_candidates = [m for m in metric_candidates if m in df.columns and _is_numeric_col(df[m])]
    else:
        preferred = [
            "HitScore",
            "mean_log2FC_BEAD","mean_log2FC_DEL2","mean_log2FC_BEAD_R2",
            "log2Boost_R2vsR1","mean_log2Boost_R2vsR1_paired",
            "var_penalty","avg_R1","avg_R2","sd_R1",
            "p_BoostPaired","q_BoostPaired","p_BEAD","q_BEAD","p_DEL2","q_DEL2","p_BEAD_R2","q_BEAD_R2"
        ]
        metric_candidates = [m for m in preferred if (m in df.columns and _is_numeric_col(df[m]))]
        numerics = [c for c in df.columns if _is_numeric_col(df[c])]
        for c in numerics:
            if c not in metric_candidates:
                metric_candidates.append(c)
            if len(metric_candidates) >= 40:
                break
    if not metric_candidates:
        raise ValueError("색/크기 기준으로 사용할 수치 메트릭을 찾지 못했습니다.")

    # 안전한 필드명으로 매핑
    metric_field_map: Dict[str, str] = {}
    for m in metric_candidates:
        safe = f"M__{sanitize_field_name(m)}"
        metric_field_map[m] = safe
        df[safe] = pd.to_numeric(df[m], errors="coerce")

    default_metric = "HitScore" if "HitScore" in metric_candidates else metric_candidates[0]
    default_safe = metric_field_map[default_metric]

    # 크기/색 스케일: 각 메트릭의 1%~99% 분위수
    metric_minmax: Dict[str, Dict[str, float]] = {}
    for m in metric_candidates:
        v = pd.to_numeric(df[metric_field_map[m]], errors="coerce").replace([np.inf,-np.inf], np.nan).dropna()
        if len(v) == 0:
            lo, hi = 0.0, 1.0
        else:
            lo = float(np.nanpercentile(v, 1))
            hi = float(np.nanpercentile(v, 99))
            if not np.isfinite(lo): lo = float(np.nanmin(v))
            if not np.isfinite(hi): hi = float(np.nanmax(v))
            if lo == hi:
                hi = lo + (1.0 if lo == 0 else abs(lo)*0.01 + 1e-6)
        metric_minmax[m] = {"low": lo, "high": hi}

    # ColumnDataSource
    df_cds = df.reset_index(drop=True).copy()
    df_cds["value"] = pd.to_numeric(df_cds[default_safe], errors="coerce").fillna(0.0)
    lo, hi = metric_minmax[default_metric]["low"], metric_minmax[default_metric]["high"]
    rng = hi - lo if (hi > lo) else 1.0
    t = (df_cds["value"] - lo) / rng
    t = t.clip(0, 1)
    df_cds["size"] = float(args.min_dot) + t * (float(args.max_dot) - float(args.min_dot))
    df_cds["_rank_default"] = df_cds["value"].rank(ascending=False, method="first")

    cds_cols = [
        "id","LibID","cycles","BB1","BB2","BB3","BB4",
        "BB1_disp","BB2_disp","BB3_disp","BB4_disp",
        "SMILES1","SMILES2","SMILES3","SMILES4",
        "IMG1","IMG2","IMG3","IMG4","CombinedID","DisplayID",
        *metric_field_map.values(), "value", "size", "_rank_default",
    ]
    source = ColumnDataSource(data=df_cds[cds_cols].to_dict("list"), name="main_source")

    # BooleanFilter & View (deprecation 해결: filter=...)
    filter_obj = BooleanFilter([True] * len(df_cds))
    view = CDSView(filter=filter_obj)

    # 색상 매핑
    palette = Blues9[::-1]
    color_mapper = LinearColorMapper(palette=palette, low=lo, high=hi)
    # Use a single shared color_mapper for all plots so low/high updates reflect everywhere

    # Hover
    TOOLTIPS = (
        "<div>"
        "<b>LibID:</b> @LibID | <b>cycles:</b> @cycles<br>"
        "<b>BB1:</b> @BB1_disp | <b>BB2:</b> @BB2_disp | <b>BB3:</b> @BB3_disp | <b>BB4:</b> @BB4_disp<br>"
        "<b>Metric value:</b> @value{0.000} (<i><b>metric</b> 드롭다운 선택 기준</i>)<br>"
        "<span style='font-size: 0.8em; color: #666;'>Click to see details (Shift/Ctrl/Cmd+Click to append)</span>"
        "</div>"
    )
    hover = HoverTool(tooltips=TOOLTIPS)

    # 6) 플롯 생성
    def make_scatter(x_col, y_col, x_cats, y_cats, title, height=None):
        ttool = TapTool()
        tools = [hover, ttool, "pan", "wheel_zoom", "box_zoom", "reset", "save", "lasso_select", "box_select"]
        p = figure(
            title=title, tools=tools, x_range=x_cats, y_range=y_cats,
            sizing_mode="stretch_width", height=int(height or args.plot_height),
            x_axis_label=f"{x_col}", y_axis_label=f"{y_col}",
            tooltips=None
        )
        p.scatter(x=x_col, y=y_col, source=source, view=view,
                  size="size", color={"field":"value","transform": color_mapper}, alpha=0.7,
                  nonselection_alpha=0.1, nonselection_color="lightgray",
                  selection_color="red", selection_alpha=0.85,
                  hover_color="orange", hover_alpha=0.9)
        p.xaxis.major_label_orientation = math.pi / 4
        p.yaxis.major_label_orientation = math.pi / 4
        p.xaxis.major_label_text_font_size = "7pt"
        p.yaxis.major_label_text_font_size = "7pt"
        p.xaxis.axis_label_text_font_size = "9pt"
        p.yaxis.axis_label_text_font_size = "9pt"
        p.title.text_font_size = "11pt"
        p.grid.grid_line_alpha = 0.3
        # Prefer Tap as the active tool so single-click selects a point
        try:
            p.toolbar.active_tap = ttool
        except Exception:
            pass
        return p

    bb1_cats = sort_bb_categories(df["BB1_disp"].unique())
    bb2_cats = sort_bb_categories(df["BB2_disp"].unique())
    bb3_cats = sort_bb_categories(df["BB3_disp"].unique())
    bb4_cats = sort_bb_categories(df["BB4_disp"].unique())

    p1 = make_scatter("BB1_disp", "BB2_disp", bb1_cats, bb2_cats, "BB1 vs BB2")
    p2 = make_scatter("BB1_disp", "BB3_disp", bb1_cats, bb3_cats, "BB1 vs BB3")
    p3 = make_scatter("BB2_disp", "BB3_disp", bb2_cats, bb3_cats, "BB2 vs BB3")
    p4 = make_scatter("BB1_disp", "BB4_disp", bb1_cats, bb4_cats, "BB1 vs BB4 (4-cycle)")

    # ColorBar
    color_bar = ColorBar(color_mapper=color_mapper, label_standoff=10, location=(0, 0),
                         title=f"{default_metric}", title_text_font_size="8pt",
                         major_label_text_font_size="7pt")
    p1.add_layout(color_bar, "left")

    # 7) 위젯/테이블/상세
    detail_init = "<p style='color:#888; font-style:italic; padding:10px;'>Click a point to see details here. <b>(Shift/Ctrl/Cmd+Click to append)</b></p>"
    detail_div = Div(
        text=detail_init,
        height=720,
        sizing_mode="stretch_width",
        styles={"overflow-y":"auto", "overflow-x":"hidden", "border":"1px solid #eee"}
    )
    clear_btn = Button(label="🧹 Clear details", width=130)
    home_btn = Button(label="↻ Home (initial)", width=160)

    libs = ["All", *sorted(pd.Series(df["LibID"]).astype(str).unique().tolist())]
    cycles_opts = ["All", "3", "4"]
    bb1_opts = ["All"] + [x for x in bb1_cats if x != "NA"]
    bb2_opts = ["All"] + [x for x in bb2_cats if x != "NA"]
    bb3_opts = ["All"] + [x for x in bb3_cats if x != "NA"]
    bb4_opts = ["All"] + [x for x in bb4_cats if x != "NA"]

    select_metric = Select(title="Metric (color/size):", value=default_metric, options=list(metric_field_map.keys()), width=260)
    select_sort   = Select(title="Sort:", value="Score↓",
                           options=["Score↓","Score↑","LibID","Cycle","BB1","BB2","BB3","BB4"], width=120)
    select_lib = Select(title="LibID:", value="All", options=libs, width=160)
    select_cycle = Select(title="Cycle:", value="All", options=cycles_opts, width=120)
    select_bb1 = Select(title="BB1:", value="All", options=bb1_opts, width=150)
    select_bb2 = Select(title="BB2:", value="All", options=bb2_opts, width=150)
    select_bb3 = Select(title="BB3:", value="All", options=bb3_opts, width=150)
    select_bb4 = Select(title="BB4:", value="All", options=bb4_opts, width=150)

    # TopN 테이블(기본 메트릭 기준)
    TOP_TABLE = int(args.top_table)
    top_df = df_cds.nsmallest(TOP_TABLE, "_rank_default")[["CombinedID", "DisplayID", "value"]].copy()
    top_df.rename(columns={"value": "metric_value"}, inplace=True)
    top_source = ColumnDataSource(top_df, name="top_source")
    initial_top_data = {k: list(v) for k, v in top_source.data.items()}

    columns = [
        TableColumn(field="DisplayID", title="ID (Lib_BB1_BB2_BB3_BB4)", width=280),
        TableColumn(field="metric_value", title=f"Value ({default_metric})",
                    formatter=NumberFormatter(format="0.0000"), width=130)
    ]
    top_table = DataTable(source=top_source, columns=columns, width=430, height=360,
                          selectable='checkbox', editable=False, sortable=True,
                          index_position=None, sizing_mode='stretch_both')

    # Optional debug panel (created once here; placed into the layout in step 8 when --debug 1)
    if int(args.debug) == 1:
        debug_panel = Div(text="<pre style='margin:0'>[debug] ready</pre>", height=120, sizing_mode="stretch_width",
                          styles={"overflow-y":"auto","border":"1px solid #eee","background":"#fafafa"})
    else:
        debug_panel = None

    # ---------------- JS 콜백 ----------------
    
    # 개선된 JavaScript 유틸리티 함수들
    js_utils = """
    // HTML escape utility
    function escapeHtml(text) {
        const map = {
            '&': '&amp;',
            '<': '&lt;',
            '>': '&gt;',
            '"': '&quot;',
            "'": '&#039;'
        };
        return String(text || '').replace(/[&<>"']/g, function(m) { return map[m]; });
    }
    
    // Sort indices based on current sort selection
    function sortIndices(indices, sortMode, data) {
        function getValue(i, field) {
            let f = field;
            if ((field === 'BB1' || field === 'BB2' || field === 'BB3' || field === 'BB4') && data[field + '_disp']) {
                f = field + '_disp';
            }
            return data[f] ? data[f][i] : '';
        }
        
        function compare(a, b) {
            switch(sortMode) {
                case 'Score↓':
                    return (+getValue(b, 'value') || 0) - (+getValue(a, 'value') || 0);
                case 'Score↑':
                    return (+getValue(a, 'value') || 0) - (+getValue(b, 'value') || 0);
                case 'LibID':
                    const libA = String(getValue(a, 'LibID') || '');
                    const libB = String(getValue(b, 'LibID') || '');
                    return libA.localeCompare(libB);
                case 'Cycle':
                    const cycA = String(getValue(a, 'cycles') || '');
                    const cycB = String(getValue(b, 'cycles') || '');
                    return cycA.localeCompare(cycB);
                default:
                    const valA = String(getValue(a, sortMode) || '');
                    const valB = String(getValue(b, sortMode) || '');
                    return valA.localeCompare(valB);
            }
        }
        
        return indices.slice().sort(compare);
    }
    
    // Generate HTML for a single compound detail block
    function generateCompoundDetailHtml(index, data, metricName) {
        const lib = escapeHtml(data['LibID'][index]);
        const cyc = escapeHtml(data['cycles'][index]);
        const val = data['value'][index];
        const formattedVal = (typeof val === 'number') ? val.toFixed(4) : escapeHtml(val);
        
        const bb1 = escapeHtml(data['BB1_disp'] ? data['BB1_disp'][index] : data['BB1'][index]);
        const bb2 = escapeHtml(data['BB2_disp'] ? data['BB2_disp'][index] : data['BB2'][index]);
        const bb3 = escapeHtml(data['BB3_disp'] ? data['BB3_disp'][index] : data['BB3'][index]);
        const bb4 = escapeHtml(data['BB4_disp'] ? data['BB4_disp'][index] : data['BB4'][index]);
        
        const img1 = data['IMG1'][index];
        const img2 = data['IMG2'][index];
        const img3 = data['IMG3'][index];
        const img4 = data['IMG4'][index];
        
        const smi1 = escapeHtml(data['SMILES1'][index]);
        const smi2 = escapeHtml(data['SMILES2'][index]);
        const smi3 = escapeHtml(data['SMILES3'][index]);
        const smi4 = escapeHtml(data['SMILES4'][index]);
        
        function generateCellHtml(img, bb, label, smiles) {
            const imgStyle = 'width:150px;height:150px;object-fit:contain;border:1px solid #ddd;margin-bottom:5px;background-color:#fff;';
            const labelStyle = 'font-size:0.9em;font-weight:bold;margin-bottom:3px;';
            const errStyle = 'color:red;font-size:0.8em;width:150px;height:150px;display:flex;align-items:center;justify-content:center;border:1px dashed #f5c6cb;background-color:#f8d7da;margin-bottom:5px;padding:5px;box-sizing:border-box;word-wrap:break-word;';
            const smilesStyle = 'font-size:0.75em;color:#555;word-break:break-all;margin-top:3px;max-width:150px;';
            
            let html = '<div style="display:flex;flex-direction:column;align-items:center;text-align:center;flex:1;min-width:160px;">';
            html += '<span style="' + labelStyle + '">' + escapeHtml(label) + ' (' + escapeHtml(bb) + ')</span>';
            
            if (img && String(img).startsWith('data:image')) {
                html += '<img src="' + img + '" style="' + imgStyle + '">';
            } else {
                const msg = img || 'No SMILES/Image';
                html += '<div style="' + errStyle + '">' + escapeHtml(msg) + '</div>';
            }
            
            html += '<span style="' + smilesStyle + '">' + escapeHtml(smiles || '') + '</span>';
            html += '</div>';
            return html;
        }
        
        let html = '<div style="border-bottom:1px solid #eee;margin-bottom:10px;padding-bottom:10px;display:flex;flex-direction:column;gap:10px;">';
        html += '<div><b>LibID:</b> ' + lib + ' | <b>cycles:</b> ' + cyc + ' | <b>' + escapeHtml(metricName) + ':</b> ' + formattedVal + '</div>';
        html += '<hr style="margin:5px 0;border:0;border-top:1px solid #eee;">';
        html += '<div style="display:flex;justify-content:space-around;align-items:flex-start;gap:10px;flex-wrap:wrap;">';
        html += generateCellHtml(img1, bb1, 'BB1', smi1);
        html += generateCellHtml(img2, bb2, 'BB2', smi2);
        html += generateCellHtml(img3, bb3, 'BB3', smi3);
        html += generateCellHtml(img4, bb4, 'BB4', smi4);
        html += '</div></div>';
        
        return html;
    }
    
    // Update detail panel with current selection
    function updateDetailPanel(selectedIndices, data, detailDiv, sortMode, metricName, initialText, accumulate = false) {
        try {
            if (!selectedIndices || selectedIndices.length === 0) {
                if (!accumulate) {
                    detailDiv.text = initialText;
                    if (window._DELE_detail_set) {
                        window._DELE_detail_set = {};
                    }
                }
                return;
            }
            
            if (!window._DELE_detail_set) {
                window._DELE_detail_set = {};
            }
            
            let html = '';
            if (accumulate && detailDiv.text && detailDiv.text.indexOf('Click a point') === -1) {
                html = detailDiv.text;
            }
            
            const sortedIndices = sortIndices(selectedIndices, sortMode, data);
            
            for (const index of sortedIndices) {
                if (!accumulate || !window._DELE_detail_set[index]) {
                    html += generateCompoundDetailHtml(index, data, metricName);
                    window._DELE_detail_set[index] = true;
                }
            }
            
            detailDiv.text = html || initialText;
            
        } catch (error) {
            console.error('Error updating detail panel:', error);
        }
    }
    """

    # 단순화된 selection change 핸들러
    selection_change_js = js_utils + """
    try {
        const selectedIndices = source.selected.indices || [];
        const data = source.data;
        const sortMode = select_sort.value;
        const metricName = select_metric.value;
        
        console.log('Selection changed:', selectedIndices.length, 'items selected');
        
        updateDetailPanel(selectedIndices, data, detail_div, sortMode, metricName, detail_init, false);
        
        // Update top table selection to match
        if (selectedIndices.length > 0) {
            const selectedCombinedIds = selectedIndices.map(i => String(data['CombinedID'][i]));
            const topData = top_source.data;
            const topIndices = [];
            
            for (let i = 0; i < (topData['CombinedID'] || []).length; i++) {
                const topCombinedId = String(topData['CombinedID'][i]);
                if (selectedCombinedIds.includes(topCombinedId)) {
                    topIndices.push(i);
                }
            }
            
            top_source.selected.indices = topIndices;
        }
        
    } catch (error) {
        console.error('Selection change error:', error);
    }
    """

    # Tap 이벤트 핸들러 (accumulate 모드 지원)
    tap_handler_js = js_utils + """
    try {
        const selectedIndices = source.selected.indices || [];
        const data = source.data;
        const sortMode = select_sort.value;
        const metricName = select_metric.value;
        
        // Check for modifier keys for accumulate mode
        const event = cb_obj;
        const accumulate = event && (
            event.modifiers && (event.modifiers.ctrl || event.modifiers.meta || event.modifiers.shift) ||
            event.ctrlKey || event.metaKey || event.shiftKey
        );
        
        console.log('Tap event:', selectedIndices.length, 'items, accumulate:', accumulate);
        
        updateDetailPanel(selectedIndices, data, detail_div, sortMode, metricName, detail_init, accumulate);
        
    } catch (error) {
        console.error('Tap event error:', error);
    }
    """

    debug_prefix = f"const DEBUG = {str(int(args.debug) == 1).lower()};\n"
    
    # Selection change 콜백 등록
    selection_cb = CustomJS(
        args=dict(
            source=source, 
            detail_div=detail_div, 
            select_sort=select_sort, 
            select_metric=select_metric,
            top_source=top_source,
            detail_init=detail_init
        ), 
        code=debug_prefix + selection_change_js
    )
    
    # Tap 이벤트 콜백
    tap_cb = CustomJS(
        args=dict(
            source=source, 
            detail_div=detail_div, 
            select_sort=select_sort, 
            select_metric=select_metric,
            detail_init=detail_init
        ), 
        code=debug_prefix + tap_handler_js
    )

    # 이벤트 핸들러 등록
    try:
        # NOTE: a single registration on selected.indices; js_on_change('selected', ...) was a duplicate
        # that only fires when the Selection object is replaced.
        source.selected.js_on_change('indices', selection_cb)
    except Exception as e:
        print(f"Warning: Could not register selection callbacks: {e}")

    # 각 플롯에 Tap 이벤트 등록
    for plot in [p1, p2, p3, p4]:
        try:
            plot.js_on_event(Tap, tap_cb)
        except Exception as e:
            print(f"Warning: Could not register tap callback: {e}")

    # (B) 필터 변경 → BooleanFilter 갱신 + TopN 갱신 + 상세 초기화
    filter_js = """
    const N = source.data["LibID"].length;
    const lib_v = select_lib.value;
    const cyc_v = select_cycle.value;
    const b1_v = select_bb1.value, b2_v = select_bb2.value, b3_v = select_bb3.value, b4_v = select_bb4.value;

    const L = source.data["LibID"], C = source.data["cycles"];
    const B1 = source.data["BB1_disp"] || source.data["BB1"];
    const B2 = source.data["BB2_disp"] || source.data["BB2"];
    const B3 = source.data["BB3_disp"] || source.data["BB3"];
    const B4 = source.data["BB4_disp"] || source.data["BB4"];

    const flags = Array(N).fill(true);
    for (let i=0;i<N;i++){
      if (lib_v !== "All" && L[i] !== lib_v) { flags[i] = false; continue; }
      if (cyc_v !== "All" && (""+C[i]) !== cyc_v) { flags[i] = false; continue; }
      if (b1_v !== "All" && B1[i] !== b1_v) { flags[i] = false; continue; }
      if (b2_v !== "All" && B2[i] !== b2_v) { flags[i] = false; continue; }
      if (b3_v !== "All" && B3[i] !== b3_v) { flags[i] = false; continue; }
      if (b4_v !== "All" && B4[i] !== b4_v) { flags[i] = false; continue; }
    }
    filter_obj.booleans = flags;

    // TopN 재계산(현재 metric value 기준, 필터 통과만)
    let idxs = [];
    for (let i=0;i<N;i++) if (flags[i]) idxs.push(i);
    idxs.sort((a,b)=> (source.data["value"][b] - source.data["value"][a]) );
    const K = Math.min(TOP_TABLE, idxs.length);
    const new_data = { CombinedID:[], DisplayID:[], metric_value:[] };
    for (let k=0;k<K;k++){
      const i = idxs[k];
      new_data["CombinedID"].push( source.data["CombinedID"][i] );
      new_data["DisplayID"].push( source.data["DisplayID"][i] );
      new_data["metric_value"].push( source.data["value"][i] );
    }
    top_source.data = new_data; top_source.change.emit();

    // 선택/상세 초기화
    source.selected.indices = [];
    detail_div.text = detail_init;
    if (window._DELE_detail_set) window._DELE_detail_set = {};
    """
    filter_cb = CustomJS(args=dict(
        source=source, filter_obj=filter_obj, detail_div=detail_div,
        select_lib=select_lib, select_cycle=select_cycle,
        select_bb1=select_bb1, select_bb2=select_bb2, select_bb3=select_bb3, select_bb4=select_bb4,
        top_source=top_source, TOP_TABLE=TOP_TABLE, detail_init=detail_init
    ), code=filter_js)
    for w in (select_lib, select_cycle, select_bb1, select_bb2, select_bb3, select_bb4):
        w.js_on_change("value", filter_cb)

    # (C) 메트릭 변경: value/색/크기 갱신 + TopN 재계산
    metric_js = """
    const label = select_metric.value;
    const safe  = metric_field_map[label];
    const data  = source.data;

    const low  = metric_minmax[label]["low"];
    const high = metric_minmax[label]["high"];
    const min_dot = SIZE_P.min_dot, max_dot = SIZE_P.max_dot;

    const N = data[ safe ].length;
    for (let i=0;i<N;i++){
      const v = +data[ safe ][i];
      data["value"][i] = (isFinite(v) ? v : 0.0);
      let t = (data["value"][i] - low) / (high - low);
      if (!isFinite(t)) t = 0.5;
      t = Math.max(0, Math.min(1, t));
      data["size"][i] = min_dot + t * (max_dot - min_dot);
    }
    source.change.emit();
    color_mapper.low  = low;
    color_mapper.high = high;

    // 필터 통과만 TopN 갱신
    const flags = filter_obj.booleans;
    let idxs = [];
    for (let i=0;i<N;i++) if (flags[i]) idxs.push(i);
    idxs.sort((a,b)=> (data["value"][b] - data["value"][a]) );
    const K = Math.min(TOP_TABLE, idxs.length);
    const new_data = { CombinedID:[], DisplayID:[], metric_value:[] };
    for (let k=0;k<K;k++){
      const i = idxs[k];
      new_data["CombinedID"].push( data["CombinedID"][i] );
      new_data["DisplayID"].push( data["DisplayID"][i] );
      new_data["metric_value"].push( data["value"][i] );
    }
    top_source.data = new_data; top_source.change.emit();

    // 선택/상세 초기화
    source.selected.indices = [];
    detail_div.text = detail_init;
    if (window._DELE_detail_set) window._DELE_detail_set = {};

    // ColorBar 타이틀 업데이트
    try { p1_left_colorbar.title = label; } catch(e) {}
    
    // 테이블 컬럼 타이틀 업데이트
    try { 
        if (top_table.columns && top_table.columns.length > 1) {
            top_table.columns[1].title = "Value (" + label + ")";
        }
    } catch(e) {}
    """
    metric_cb = CustomJS(args=dict(
        source=source, filter_obj=filter_obj, detail_div=detail_div,
        select_metric=select_metric, metric_field_map=metric_field_map,
        metric_minmax=metric_minmax, color_mapper=color_mapper,
        top_source=top_source, TOP_TABLE=TOP_TABLE, top_table=top_table,
        SIZE_P=dict(min_dot=float(args.min_dot), max_dot=float(args.max_dot)),
        p1_left_colorbar=color_bar, detail_init=detail_init
    ), code=metric_js)
    select_metric.js_on_change("value", metric_cb)
    
    # Sort change: re-render details based on current selection
    sort_change_js = js_utils + """
    try {
        const selectedIndices = source.selected.indices || [];
        const data = source.data;
        const sortMode = select_sort.value;
        const metricName = select_metric.value;
        
        if (selectedIndices.length > 0) {
            updateDetailPanel(selectedIndices, data, detail_div, sortMode, metricName, detail_init, false);
        }
    } catch (error) {
        console.error('Sort change error:', error);
    }
    """
    
    select_sort.js_on_change("value", CustomJS(
        args=dict(
            source=source, 
            detail_div=detail_div, 
            select_metric=select_metric, 
            select_sort=select_sort,
            detail_init=detail_init
        ), 
        code=debug_prefix + sort_change_js
    ))

    # (D) TopN 테이블 선택 → 상세 패널 누적 + 플롯 하이라이트(누적)
    table_select_js = js_utils + """
    try {
        const tableIndices = top_source.selected.indices || [];
        const tableData = top_source.data;
        const plotData = source.data;
        
        // Map table selections to plot indices
        const selectedPlotIndices = [];
        for (const tableIdx of tableIndices) {
            const combinedId = tableData['CombinedID'][tableIdx];
            for (let i = 0; i < plotData['CombinedID'].length; i++) {
                if (plotData['CombinedID'][i] === combinedId) {
                    selectedPlotIndices.push(i);
                    break;
                }
            }
        }
        
        // Accumulate with existing selections
        const existingIndices = source.selected.indices || [];
        const allIndices = [...new Set([...existingIndices, ...selectedPlotIndices])];
        
        source.selected.indices = allIndices;
        
        const sortMode = select_sort.value;
        const metricName = select_metric.value;
        
        updateDetailPanel(allIndices, plotData, detail_div, sortMode, metricName, detail_init, true);
        
    } catch (error) {
        console.error('Table selection error:', error);
    }
    """
    
    table_select_cb = CustomJS(
        args=dict(
            top_source=top_source, 
            source=source, 
            detail_div=detail_div, 
            select_metric=select_metric, 
            select_sort=select_sort,
            detail_init=detail_init
        ), 
        code=debug_prefix + table_select_js
    )
    top_source.selected.js_on_change("indices", table_select_cb)

    # (E) Reset 액션
    reset_js = """
    select_lib.value = 'All'; select_cycle.value = 'All';
    select_bb1.value = 'All'; select_bb2.value = 'All'; select_bb3.value = 'All'; select_bb4.value = 'All';
    filter_obj.booleans = Array(source.data["LibID"].length).fill(true);
    filter_obj.change.emit();
    source.selected.indices = [];
    top_source.selected.indices = [];
    detail_div.text = detail_init;
    if (window._DELE_detail_set) window._DELE_detail_set = {};
    p1.reset.emit(); p2.reset.emit(); p3.reset.emit(); p4.reset.emit();
    """
    reset_cb = CustomJS(args=dict(
        source=source, filter_obj=filter_obj, detail_div=detail_div,
        select_lib=select_lib, select_cycle=select_cycle,
        select_bb1=select_bb1, select_bb2=select_bb2, select_bb3=select_bb3, select_bb4=select_bb4,
        p1=p1, p2=p2, p3=p3, p4=p4, top_source=top_source, detail_init=detail_init
    ), code=reset_js)
    reset_action = CustomAction(description="Reset Filters, Selection, Zoom & Details", callback=reset_cb)
    for p in (p1, p2, p3, p4):
        p.tools.append(reset_action)

    # (F) Clear / Home 버튼
    clear_btn.js_on_click(CustomJS(
        args=dict(detail_div=detail_div, detail_init=detail_init),
        code="""
        if (window._DELE_detail_set) window._DELE_detail_set = {}; 
        detail_div.text = detail_init;
        """
    ))
    
    home_js = """
    select_metric.value = DEFAULT_METRIC;
    select_sort.value = 'Score↓';
    (function() {
      const label = DEFAULT_METRIC;
      const safe  = metric_field_map[label];
      const data  = source.data;
      const low  = metric_minmax[label]["low"];
      const high = metric_minmax[label]["high"];
      const min_dot = SIZE_P.min_dot, max_dot = SIZE_P.max_dot;
      const N = data[ safe ].length;
      for (let i=0;i<N;i++){
        const v = +data[ safe ][i];
        data["value"][i] = (isFinite(v) ? v : 0.0);
        let t = (data["value"][i] - low) / (high - low);
        if (!isFinite(t)) t = 0.5;
        t = Math.max(0, Math.min(1, t));
        data["size"][i] = min_dot + t * (max_dot - min_dot);
      }
      source.change.emit();
      color_mapper.low  = low;
      color_mapper.high = high;
      try { p1_left_colorbar.title = label; } catch(e) {}
    })();

    select_lib.value = 'All'; select_cycle.value = 'All';
    select_bb1.value = 'All'; select_bb2.value = 'All'; select_bb3.value = 'All'; select_bb4.value = 'All';
    const allTrue = Array(source.data["LibID"].length).fill(true);
    filter_obj.booleans = allTrue;
    filter_obj.change.emit();

    const cloned = {};
    for (const key in INITIAL_TOP_DATA) {
      if (!Object.prototype.hasOwnProperty.call(INITIAL_TOP_DATA, key)) continue;
      const value = INITIAL_TOP_DATA[key];
      cloned[key] = Array.isArray(value) ? value.slice() : value;
    }
    top_source.data = cloned;
    top_source.change.emit();

    source.selected.indices = [];
    top_source.selected.indices = [];
    if (window._DELE_detail_set) window._DELE_detail_set = {};
    detail_div.text = detail_init;
    p1.reset.emit(); p2.reset.emit(); p3.reset.emit(); p4.reset.emit();
    """
    home_btn.js_on_click(CustomJS(args=dict(
        source=source, select_metric=select_metric, color_mapper=color_mapper,
        metric_field_map=metric_field_map, metric_minmax=metric_minmax,
        select_lib=select_lib, select_cycle=select_cycle,
        select_bb1=select_bb1, select_bb2=select_bb2, select_bb3=select_bb3, select_bb4=select_bb4,
        filter_obj=filter_obj, top_source=top_source, detail_div=detail_div,
        INITIAL_TOP_DATA=initial_top_data, DEFAULT_METRIC=default_metric,
        SIZE_P=dict(min_dot=float(args.min_dot), max_dot=float(args.max_dot)),
        p1=p1, p2=p2, p3=p3, p4=p4, p1_left_colorbar=color_bar, detail_init=detail_init, select_sort=select_sort
    ), code=home_js))

    # 8) 레이아웃
    widgets_row = bk_row(
        Spacer(width=10), select_metric, select_sort, select_lib, select_cycle,
        select_bb1, select_bb2, select_bb3, select_bb4, Spacer(width=10),
        sizing_mode="stretch_width"
    )
    plots_layout = gridplot([[p1, p2], [p3, p4]], sizing_mode="stretch_width", merge_tools=False, toolbar_location="right")
    
    control_bar = bk_row(clear_btn, home_btn, sizing_mode="stretch_width")
    if debug_panel is not None:
        detail_col = bk_column(control_bar, debug_panel, detail_div, sizing_mode="stretch_both")
    else:
        detail_col = bk_column(control_bar, detail_div, sizing_mode="stretch_both")
    bottom_row = bk_row(detail_col, top_table, sizing_mode="stretch_width")
    final_layout = bk_column(widgets_row, plots_layout, bottom_row, sizing_mode="stretch_both")

    # 8.1) Responsive sizing for main page
    resize_main_js = """
      function _set(){
        try {
          const H = window.innerHeight || document.documentElement.clientHeight || 900;
          // Two rows of plots (top): scale by viewport
          const plotH = Math.max(220, Math.floor(H * 0.22));
          p1.height = plotH; p2.height = plotH; p3.height = plotH; p4.height = plotH;
          // Compute available height from the top of the detail panel down to the bottom of viewport
          let topY = 0;
          try { topY = detail_div.el.getBoundingClientRect().top; } catch(e) { topY = Math.floor(H*0.5); }
          const avail = Math.max(260, Math.floor(H - topY - 8));
          detail_div.height = avail;
          top_table.height = avail;
        } catch(e) {}
      }
      _set();
      try { window.addEventListener('resize', _set); } catch(e) {}
    """
    p1.js_on_event(DocumentReady, CustomJS(args=dict(p1=p1,p2=p2,p3=p3,p4=p4, detail_div=detail_div, top_table=top_table), code=resize_main_js))

    # 9) 출력
    print(f"[{time.strftime('%H:%M:%S')}] Write HTML: {args.out}")
    output_file(args.out, title="DELeGANce | Master Top-N Explorer (HitScore ranking)")
    save(final_layout)  # batch-safe: write only, never spawn a browser (show() could block on headless hosts)

    print(f"[{time.strftime('%H:%M:%S')}] Done.")
    print("=" * 90)


if __name__ == "__main__":
    main()

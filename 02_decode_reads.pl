#!/usr/bin/env perl
# 02_decode_reads.pl (refactored to unified style)
# - 기능/논리/입출력 규격: 변경 없음
# - 형태 통일점:
#   * Embedded LibraryTrie (공통) + die_on_dup=1 로 중복 삽입 방지 동작 유지
#   * ensure_dir: Windows/Unix 안전
#   * 로깅(ts/log_msg) + 경고/예외 트랩 공통
#   * 섹션 레이아웃/주석 규약 통일
# - 설계 고정:
#   * Codon은 perfect만, HP/OP/CP 1-bp mismatch 허용(치환/삽입/삭제; --mismatch로 제어)
#   * 길이 윈도우 정책/fastp 기반 Top-K 피크 감지/멀티 윈도우 로직 그대로
#   * 헤더/출력 컬럼 수(디코딩=35, 언디코딩=7) 검증 그대로

use strict;
use warnings;

# ★★★ 중요: 단축키 대소문자 구분(충돌 방지) ★★★
BEGIN { require Getopt::Long; Getopt::Long::Configure(qw(no_ignore_case)); }

use File::Spec;
use File::Basename qw(basename);
use FindBin qw($Bin);
use IO::Handle;
use Getopt::Long qw(GetOptions);
use JSON::PP;
use Scalar::Util qw(looks_like_number);

############################
# LibraryTrie.pm (embedded, unified)
############################
package LibraryTrie;
use strict; use warnings;

sub new {
    my ($class, %opt) = @_;
    my $die_on_dup = $opt{die_on_dup} // 0;
    return bless { children => {}, is_end => 0, values => {}, die_on_dup => $die_on_dup }, $class;
}
sub insert {
    my ($self, $sequence, $info, $lib_id) = @_;
    die "Error: 'info' must be defined.\n"   unless defined $info;
    die "Error: 'lib_id' must be defined.\n" unless defined $lib_id;
    my $node=$self;
    foreach my $ch (split //, $sequence) {
        $node->{children}{$ch} //= { children => {}, is_end => 0, values => {} };
        $node = $node->{children}{$ch};
    }
    $node->{is_end} = 1;
    if ($self->{die_on_dup}) {
        die "Error: dup sequence for lib_id in trie: $sequence / $lib_id\n"
            if exists $node->{values}{$lib_id} && @{$node->{values}{$lib_id}};
    }
    $node->{values}{$lib_id} //= [];
    push @{ $node->{values}{$lib_id} }, $info;
}
# 치환 기반 1-bp mismatch(디코딩에서는 자체 indel 생성 로직 사용; 필요 시 보조로 사용 가능)
sub insert_one_mismatch {
    my ($self, $sequence, $info, $lib_id) = @_;
    my @chars = split //, $sequence;
    my @bases = ('A','C','G','T');
    for (my $i = 0; $i < @chars; $i++) {
        my $orig = $chars[$i];
        foreach my $b (@bases) {
            next if $b eq $orig;
            my @mut = @chars; $mut[$i] = $b;
            my $mutseq = join("", @mut);
            $self->insert($mutseq, $info, $lib_id);
        }
    }
}
sub search_substrings {
    my ($self,$s)=@_;
    my @m;
    for (my $i=0; $i<length($s); $i++){
        my $node=$self; my $j=$i;
        while ($j<length($s)){
            my $c=substr($s,$j,1);
            last unless exists $node->{children}{$c};
            $node=$node->{children}{$c};
            if ($node->{is_end}){
                foreach my $lib_id (keys %{ $node->{values} }){
                    my $aref = $node->{values}{$lib_id} // [];
                    for my $info (@$aref){
                        push @m, { sequence=>substr($s,$i,$j-$i+1), lib_id=>$lib_id, info=>$info, position=>$i };
                    }
                }
            }
            $j++;
        }
    }
    return \@m;
}
sub search_substrings_by_lib_id {
    my ($self,$s,$target_lib)=@_;
    die "target_lib_id must be defined\n" unless defined $target_lib;
    my @m;
    for (my $i=0; $i<length($s); $i++){
        my $node=$self; my $j=$i;
        while ($j<length($s)){
            my $c=substr($s,$j,1);
            last unless exists $node->{children}{$c};
            $node=$node->{children}{$c};
            if ($node->{is_end} && exists $node->{values}{$target_lib}){
                my $aref = $node->{values}{$target_lib} // [];
                for my $info (@$aref){
                    push @m, { sequence=>substr($s,$i,$j-$i+1), lib_id=>$target_lib, info=>$info, position=>$i };
                }
            }
            $j++;
        }
    }
    return \@m;
}

############################
# main
############################
package main;
use strict; use warnings;

############################
# Usage / CLI parsing
############################
my %OPT = (
    target_protein     => $ENV{TARGET_PROTEIN}     // 'DELeGANce_out',
    fastp_outdir_name  => $ENV{FASTP_OUTDIR_NAME}  // '01_fastp_out',
    merged_dir         => undef,
    fixed_bb_file      => undef,
    out_dir            => undef,
    length_filter_mode => $ENV{LENGTH_FILTER_MODE} // 'reject_outside',
    length_tol         => defined $ENV{LENGTH_TOL} ? $ENV{LENGTH_TOL}+0 : 5,
    scaling_denom      => $ENV{SCALING_DENOM}      // 'sample_decoded',
    max_failed_dump    => defined $ENV{MAX_FAILED_DUMP} ? $ENV{MAX_FAILED_DUMP}+0 : 100000,
    adj_tol            => defined $ENV{ADJ_TOL} ? $ENV{ADJ_TOL}+0 : 3,
    max_cp_cands       => defined $ENV{MAX_CP_CANDS} ? $ENV{MAX_CP_CANDS}+0 : 0,
    max_anchor_cands   => defined $ENV{MAX_ANCHOR_CANDS} ? $ENV{MAX_ANCHOR_CANDS}+0 : 0,
    mismatch           => $ENV{MISMATCH_MODE} // 'hp_op_cp',

    # tag canonical lens (override)
    hp_len             => undef,   # default 8
    op_len             => undef,   # default 20
    codon_len          => undef,   # default 9
    cp_len             => undef,   # default 28

    # fastq matching
    fastq_regex        => undef,

    # fastp auto-detect
    auto_insert_detect => 1,
    auto_select_window => 0,   # -W 1: pick a single 3/4-cycle window (or 'both') from fastp support instead of Top-K 'multi'
    auto_center_window => 0,   # -C 1: (with -W 1) re-centre the selected window on the nearest fastp peak
    auto_tol_from_fastp=> 0,

    # Top-K peaks for mixed DEL
    peak_top_k         => 5,       # --peak-top-k / -K
    peak_min_sep       => 2,       # --peak-min-sep / -E (nt)
    peak_min_frac      => 0.01,    # --peak-min-frac / -F (>= 1% of hist total)

    help               => 0,
);

my $NO_CLI_ARGS = (@ARGV == 0) ? 1 : 0;   # bare invocation (no options at all) -> usage
GetOptions(
    'target-protein|p=s'     => \$OPT{target_protein},
    'fastp-outdir-name|n=s'  => \$OPT{fastp_outdir_name},
    'merged-dir|m=s'         => \$OPT{merged_dir},
    'fixed-bb-file|b=s'      => \$OPT{fixed_bb_file},
    'out-dir|o=s'            => \$OPT{out_dir},
    'length-filter-mode|M=s' => \$OPT{length_filter_mode},
    'length-tol|T=i'         => \$OPT{length_tol},
    'scaling-denom|S=s'      => \$OPT{scaling_denom},
    'max-failed-dump|x=i'    => \$OPT{max_failed_dump},
    'adj-tol|J=i'            => \$OPT{adj_tol},
    'max-cp-cands|P=i'       => \$OPT{max_cp_cands},
    'max-anchor-cands|Y=i'   => \$OPT{max_anchor_cands},
    'mismatch=s'             => \$OPT{mismatch},
    'hp-len|H=i'             => \$OPT{hp_len},
    'op-len|O=i'             => \$OPT{op_len},
    'codon-len|c=i'          => \$OPT{codon_len},
    'cp-len|Q=i'             => \$OPT{cp_len},
    'fastq-regex|R=s'        => \$OPT{fastq_regex},

    'auto-insert-detect|A=i'  => \$OPT{auto_insert_detect},
    'auto-select-window|W=i'  => \$OPT{auto_select_window},
    'auto-center-window|C=i'  => \$OPT{auto_center_window},
    'auto-tol-from-fastp|U=i' => \$OPT{auto_tol_from_fastp},

    'peak-top-k|K=i'          => \$OPT{peak_top_k},
    'peak-min-sep|E=i'        => \$OPT{peak_min_sep},
    'peak-min-frac|F=f'       => \$OPT{peak_min_frac},

    'help|h'                  => \$OPT{help},
) or die "Error parsing options. Use -h for usage.\n";

# 무옵션 실행 시 설명서 출력 (exit 2: 인자 부족). -p/-o 등 어떤 옵션이든 주어지면 -m/-b 는 기본 경로로 진행한다.
if ($OPT{help}) { print_usage_and_exit(0); }
if ($NO_CLI_ARGS) { print_usage_and_exit(2); }

# Validate
my %VALID_DEN = map { $_=>1 } qw(sample_decoded sample_length_passed sample_total_raw library_decoded);
die "Invalid --scaling-denom: $OPT{scaling_denom}\nAllowed: sample_decoded | sample_length_passed | sample_total_raw | library_decoded\n"
  unless $VALID_DEN{$OPT{scaling_denom}};
die "Invalid --mismatch: $OPT{mismatch} (use none|hp_op_cp)\n"
  unless $OPT{mismatch} =~ /^(?:none|hp_op_cp)$/;
my %VALID_LFM = map { $_=>1 } qw(reject_outside reject_inside);
die "Invalid --length-filter-mode: $OPT{length_filter_mode}\nAllowed: reject_outside | reject_inside\n"
  unless $VALID_LFM{$OPT{length_filter_mode}};
die "--length-tol must be >=0\n" if defined $OPT{length_tol} && $OPT{length_tol} < 0;
die "--max-failed-dump must be >=0\n" if defined $OPT{max_failed_dump} && $OPT{max_failed_dump} < 0;
die "--adj-tol must be >=0\n" if defined $OPT{adj_tol} && $OPT{adj_tol} < 0;
die "--max-cp-cands must be >=0\n" if defined $OPT{max_cp_cands} && $OPT{max_cp_cands} < 0;
die "--max-anchor-cands must be >=0\n" if defined $OPT{max_anchor_cands} && $OPT{max_anchor_cands} < 0;
die "--peak-top-k must be >=1\n" if $OPT{peak_top_k} < 1;
die "--peak-min-sep must be >=0\n" if $OPT{peak_min_sep} < 0;
die "--peak-min-frac must be between 0 and 1\n" if $OPT{peak_min_frac} < 0 || $OPT{peak_min_frac} > 1;

# --------- Paths / I/O ----------
my $TARGET_PROTEIN     = $OPT{target_protein};
my $FASTP_OUTDIR_NAME  = $OPT{fastp_outdir_name};

my $MERGED_DIR         = defined $OPT{merged_dir} ? $OPT{merged_dir}
                        : File::Spec->catdir($Bin, $TARGET_PROTEIN, $FASTP_OUTDIR_NAME);
my $FIXED_BB_FILE      = defined $OPT{fixed_bb_file} ? $OPT{fixed_bb_file}
                        : File::Spec->catfile($Bin, $TARGET_PROTEIN, "BB_information_fixed.tsv");
my $OUT_DIR            = defined $OPT{out_dir} ? $OPT{out_dir}
                        : File::Spec->catdir($Bin, $TARGET_PROTEIN, "02_decoded");
my $LOG_FILE           = File::Spec->catfile($OUT_DIR, "02_decode_reads.log");

my $RAW_MATRIX_FILE    = File::Spec->catfile($OUT_DIR, "raw_counts_matrix.tsv");
my $SCALED_MATRIX_FILE = File::Spec->catfile($OUT_DIR, "scaled_counts_matrix.tsv");
my $LIB_TOTALS_FILE    = File::Spec->catfile($OUT_DIR, "library_totals.tsv");
my $SUMMARY_FILE       = File::Spec->catfile($OUT_DIR, "decoding_summary.tsv");

# 추가 통계/검증 파일
my $SAMPLE_STATS_FILE      = File::Spec->catfile($OUT_DIR, "sample_stats.tsv");
my $LIB_SAMPLE_STATS_FILE  = File::Spec->catfile($OUT_DIR, "lib_sample_stats.tsv");
my $ANCHOR_QUAL_FILE       = File::Spec->catfile($OUT_DIR, "anchor_quality_counts.tsv");
my $GAP_STATS_FILE         = File::Spec->catfile($OUT_DIR, "gap_stats.tsv");
my $CYCLE_DIST_FILE        = File::Spec->catfile($OUT_DIR, "cycle_distribution.tsv");
my $BB_FREQ_BY_CYCLE_FILE  = File::Spec->catfile($OUT_DIR, "bb_frequency_by_cycle.tsv");
my $LEN_STATS_FILE         = File::Spec->catfile($OUT_DIR, "length_window_stats.tsv");
my $QC_CHECKS_FILE         = File::Spec->catfile($OUT_DIR, "qc_checks.tsv");

# 실패 리드 덤프 제한
my $MAX_FAILED_DUMP   = $OPT{max_failed_dump} // 100000;

# 길이 필터 정책/오차 (기본)
my $LENGTH_FILTER_MODE = $OPT{length_filter_mode};  # 'reject_outside' or 'reject_inside'
my $LENGTH_TOL_BASE    = $OPT{length_tol}+0;   # CLI value; never mutated
my $LENGTH_TOL         = $LENGTH_TOL_BASE;      # per-sample working value (may be changed by -U 1)

# 설계상 Canonical 길이 (옵션 오버라이드 지원)
my $HP_CAN_LEN    = defined $OPT{hp_len}    ? $OPT{hp_len}+0    : 8;
my $OP_CAN_LEN    = defined $OPT{op_len}    ? $OPT{op_len}+0    : 20;
my $CODON_CAN_LEN = defined $OPT{codon_len} ? $OPT{codon_len}+0 : 9;
my $CP_CAN_LEN    = defined $OPT{cp_len}    ? $OPT{cp_len}+0    : 28;
my $ADJ_TOL       = defined $OPT{adj_tol}  ? $OPT{adj_tol}+0   : 3;
my $MAX_CP_CANDS     = $OPT{max_cp_cands}    // 0;  # 0 = no limit
my $MAX_ANCHOR_CANDS = $OPT{max_anchor_cands}// 0;  # 0 = no limit

die "--hp-len must be >0\n"    if $HP_CAN_LEN    <= 0;
die "--op-len must be >0\n"    if $OP_CAN_LEN    <= 0;
die "--codon-len must be >0\n" if $CODON_CAN_LEN <= 0;
die "--cp-len must be >0\n"    if $CP_CAN_LEN    <= 0;

# 기본(설계) 기대 길이
my $EXP_LEN_3_BASE = $HP_CAN_LEN + $OP_CAN_LEN + 3*$CODON_CAN_LEN + $CP_CAN_LEN;  # e.g., 83
my $EXP_LEN_4_BASE = $HP_CAN_LEN + $OP_CAN_LEN + 4*$CODON_CAN_LEN + $CP_CAN_LEN;  # e.g., 92

# 실제 사용 길이 (샘플별로 변경될 수 있음)
my $EXP_LEN_3 = $EXP_LEN_3_BASE;
my $EXP_LEN_4 = $EXP_LEN_4_BASE;

# fastp 자동 감지 정책
my $SCALING_DENOM          = $OPT{scaling_denom};
my $AUTO_INSERT_DETECT     = $OPT{auto_insert_detect} ? 1 : 0;
my $AUTO_SELECT_WINDOW     = $OPT{auto_select_window} ? 1 : 0;
my $AUTO_CENTER_WINDOW     = $OPT{auto_center_window} ? 1 : 0;
my $AUTO_TOL_FROM_FASTP    = $OPT{auto_tol_from_fastp} ? 1 : 0;

# Top-K peak 추출 파라미터
my $PEAK_TOP_K    = $OPT{peak_top_k}+0;
my $PEAK_MIN_SEP  = $OPT{peak_min_sep}+0;
my $PEAK_MIN_FRAC = $OPT{peak_min_frac}+0.0;

# --------- Logging ----------
our $LOG_FH;
sub ensure_dir {
    my ($dir) = @_;
    return if -d $dir;
    # File::Path avoids interpolating the path into a shell command (quotes/special chars were unsafe)
    require File::Path;
    File::Path::make_path($dir, { error => \my $err });
    if ($err && @$err) { die "Failed to create directory: $dir\n"; }
}
ensure_dir($OUT_DIR);
open($LOG_FH, ">>", $LOG_FILE) or die "Cannot open log $LOG_FILE: $!\n";
$LOG_FH->autoflush(1); STDOUT->autoflush(1);
# Stage completion marker (log is append-mode; an old "All done." must not count as completion)
my $DONE_MARKER = File::Spec->catfile($OUT_DIR, ".decode_done");
unlink($DONE_MARKER) if -e $DONE_MARKER;
sub ts { scalar localtime() }
sub log_msg { my ($lvl,$msg)=@_; my $ln=ts()." [$lvl] $msg\n"; print $ln; print $LOG_FH $ln; }
$SIG{__WARN__} = sub { my $m=join("",@_); chomp($m); log_msg("WARN",$m); };
$SIG{__DIE__}  = sub { return if $^S; my $m=join("",@_); chomp($m); log_msg("DIE",$m); CORE::die("$m\n"); };   # $^S: inside eval -> caller handles it, do not log as DIE

sub dump_effective_config {
    my %conf = (
        target_protein     => $TARGET_PROTEIN,
        fastp_outdir_name  => $FASTP_OUTDIR_NAME,
        merged_dir         => $MERGED_DIR,
        fixed_bb_file      => $FIXED_BB_FILE,
        out_dir            => $OUT_DIR,
        length_filter_mode => $LENGTH_FILTER_MODE,
        length_tol         => $LENGTH_TOL,
        scaling_denom      => $SCALING_DENOM,
        max_failed_dump    => $MAX_FAILED_DUMP,
        hp_len             => $HP_CAN_LEN,
        op_len             => $OP_CAN_LEN,
        codon_len          => $CODON_CAN_LEN,
        cp_len             => $CP_CAN_LEN,
        exp_len_3_base     => $EXP_LEN_3_BASE,
        exp_len_4_base     => $EXP_LEN_4_BASE,
        fastq_regex        => ($OPT{fastq_regex}//'(auto)'),
        auto_insert_detect => $AUTO_INSERT_DETECT,
        auto_select_window => $AUTO_SELECT_WINDOW,
        auto_center_window => $AUTO_CENTER_WINDOW,
        auto_tol_from_fastp=> $AUTO_TOL_FROM_FASTP,
        peak_top_k         => $PEAK_TOP_K,
        peak_min_sep       => $PEAK_MIN_SEP,
        peak_min_frac      => $PEAK_MIN_FRAC,
        adj_tol            => $ADJ_TOL,
        max_cp_cands       => $MAX_CP_CANDS,
        max_anchor_cands   => $MAX_ANCHOR_CANDS,
        mismatch_mode      => $OPT{mismatch},
    );
    my @pairs = map { "$_=$conf{$_}" } sort keys %conf;
    log_msg("INFO","Effective config: ".join("; ", @pairs));
}

# --------- Tries & Maps ----------
# 중복방지(die_on_dup=1) 동작은 기존 디코딩 스크립트와 동일
my $HP_TRIE    = LibraryTrie->new(die_on_dup => 1);
my $OP_TRIE    = LibraryTrie->new(die_on_dup => 1);
my $CP_TRIE    = LibraryTrie->new(die_on_dup => 1);
my $CODON_TRIE = LibraryTrie->new(die_on_dup => 1);

my $INFO_SEP = "\x1f";
# info 포맷: "type<SEP>tag_id<SEP>bb_id_fixed<SEP>smiles<SEP>perf|miss[<SEP>cycle]"
my %LIB_EXPECTED_CYC;  # {lib_id} = 3 or 4
my %LIB_IDS;           # 실제 라이브러리 set
my %CP_VARIANT_OWNER;  # {seq} = lib_id (first seen) to detect CP collisions across libs

# --------- 1bp 변이 유틸 ----------
sub one_bp_deletion  { my ($s)=@_; my @o; for my $i(0..length($s)-1){ push @o, substr($s,0,$i).substr($s,$i+1) } @o }
sub one_bp_insertion {
    my ($s)=@_; my @o; my @nuc=qw(A G C T);
    my $L=length($s); my $first=substr($s,0,1); my $last=substr($s,-1,1);
    # Insertions at the two ends (i==0, i==L) are deliberately excluded: those strings contain the perfect
    # anchor as a substring and are already found by the substring search. Inserting $first at i==1 (or
    # $last at i==L-1) yields the same string as an end insertion, so those are skipped for the same reason.
    for my $i(1..$L-1){ for my $n(@nuc){
        next if ($i==1    && $n eq $first);
        next if ($i==$L-1 && $n eq $last);
        push @o, substr($s,0,$i).$n.substr($s,$i);
    }} @o
}
sub one_bp_substitution { my ($s)=@_; my @o; my @nuc=qw(A G C T);
    for my $i(0..length($s)-1){ my $orig=substr($s,$i,1);
        for my $n(@nuc){ next if $n eq $orig; my $t=$s; substr($t,$i,1,$n); push @o,$t; }
    } @o
}
sub remove_dups { my %s; my @u; for my $x(@_){ next if $s{$x}++; push @u,$x } @u }
sub rc { my ($s)=@_; $s=~tr/ACGT/TGCA/; scalar reverse $s }
sub mismatch_enabled_for_type {
    my ($type)=@_;
    return 0 if $type eq 'CODON';               # hard rule: CODON mismatch disabled
    return 0 if $OPT{mismatch} eq 'none';
    return 1 if $OPT{mismatch} eq 'hp_op_cp' && ($type eq 'HP' || $type eq 'OP' || $type eq 'CP');
    return 0;
}

# --------- BB 파일 로드 & Trie 구축 ----------
my %VARIANT_SKIPPED;   # type -> count of 1-bp variants dropped (collision with a perfect anchor or with another anchor's variant)
my %VARIANT_AMBIG;     # {type}{lib}{seq} = 1 : variant sequence shared by >=2 anchors of the same lib -> removed from the trie for good
my %ROWS_SKIPPED;      # reason -> count of BB rows silently ignored (reported once after loading)
my %ANCHOR_ROWS;       # type -> count of loaded perfect anchor rows

sub _trie_walk {
    # returns the node for $seq if every character exists on the path, else undef
    my ($trie, $seq) = @_;
    my $node = $trie;
    foreach my $ch (split //, $seq) {
        return undef unless exists $node->{children}{$ch};
        $node = $node->{children}{$ch};
    }
    return $node;
}
sub _info_qual { my ($info)=@_; my @f = split /\Q$INFO_SEP\E/, $info; return $f[4] // ''; }
sub _insert_perfect_checked {
    # Pass 1: perfect anchors only. The only fatal condition is a genuine duplicate of the same perfect sequence
    # within the same lib/type; report it with the offending tag ids instead of LibraryTrie's generic message.
    my ($trie, $type, $seq, $info, $lib, $tag_id, $ln) = @_;
    my $node = _trie_walk($trie, $seq);
    if ($node && $node->{is_end} && exists $node->{values}{$lib} && @{ $node->{values}{$lib} }) {
        my @prev_tags = map { (split /\Q$INFO_SEP\E/, $_)[1] // '?' } @{ $node->{values}{$lib} };
        die "Duplicate $type anchor sequence $seq in lib '$lib' (line $ln, tag $tag_id; already loaded as tag @prev_tags)\n";
    }
    $trie->insert($seq, $info, $lib);
    $ANCHOR_ROWS{$type}++;
}
sub _insert_variant_safe {
    # Pass 2: 1-bp variants. Because every perfect anchor is already in the trie, the outcome does not depend on
    # the row order of the BB file:
    #   (i)  variant == perfect sequence of another anchor (edit distance 1)  -> drop the variant, keep the perfect entry
    #   (ii) variant == variant of another anchor (edit distance 2)           -> ambiguous: drop BOTH (empty the lib's
    #        values at that node) and remember the sequence so a third anchor cannot re-insert it
    my ($trie, $type, $seq, $info, $lib, $tag_id) = @_;
    if ($VARIANT_AMBIG{$type}{$lib}{$seq}) { $VARIANT_SKIPPED{$type}++; return 0; }
    my $node = _trie_walk($trie, $seq);
    if ($node && $node->{is_end} && exists $node->{values}{$lib} && @{ $node->{values}{$lib} }) {
        my @vals = @{ $node->{values}{$lib} };
        if (grep { _info_qual($_) eq 'perf' } @vals) {
            $VARIANT_SKIPPED{$type}++;
            log_msg("WARN", "$type 1-bp variant '$seq' of tag $tag_id (lib '$lib') equals a perfect $type anchor; variant dropped (anchors within edit distance 1)")
                if $VARIANT_SKIPPED{$type} <= 20;
            return 0;
        }
        # only variants of other anchor(s) live here -> ambiguous sequence, remove them too
        my @other_tags = remove_dups(map { (split /\Q$INFO_SEP\E/, $_)[1] // '?' } @vals);
        delete $node->{values}{$lib};
        $node->{is_end} = 0 unless %{ $node->{values} };
        $VARIANT_AMBIG{$type}{$lib}{$seq} = 1;
        $VARIANT_SKIPPED{$type} += 1 + scalar(@vals);
        log_msg("WARN", "$type 1-bp variant '$seq' is shared by tags @other_tags and $tag_id (lib '$lib'); ambiguous variant removed for all of them (anchors within edit distance 2)")
            if $VARIANT_SKIPPED{$type} <= 20;
        return 0;
    }
    $trie->insert($seq, $info, $lib);
    return 1;
}
sub load_fixed_bb_and_build_tries {
    my ($file)=@_;
    log_msg("INFO","Loading fixed BB info: $file");
    open(my $IN,"<",$file) or die "Cannot open $file: $!\n";
    my $ln=0;
    my @anchor_rows;   # deferred: {trie,type,seq,info_p,info_m,lib,tag_id,ln} for HP/OP/CP (variants inserted in pass 2)
    my %tries = (HP => $HP_TRIE, OP => $OP_TRIE, CP => $CP_TRIE);
    while(my $line=<$IN>){
        $ln++; chomp($line);
        next if $line =~ /^\s*$/ || $line =~ /^\s*#/;
        my @d = split /\t/, $line, -1;
        if (@d < 7) {
            $ROWS_SKIPPED{short_row}++;
            log_msg("WARN", "BB line $ln has ".scalar(@d)." columns (<7); skipped") if $ROWS_SKIPPED{short_row} <= 20;
            next;
        }
        my ($type_raw,$seq_raw,$bb_id_fixed,$cycle,$tag_id,$lib_id,$smiles)=@d;
        my $type = uc($type_raw // ''); $type =~ s/^\s+|\s+$//g;   # 01_preprocess accepts type case-insensitively -> match that here
        if ($ln==1 && $type =~ /^(TYPE|BB[_-]?TYPE)$/){ next; } # header

        my $seq = uc(join("",(split /\s+/,$seq_raw)));
        unless ($seq =~ /^[ACGT]+$/) {
            $ROWS_SKIPPED{non_acgt_seq}++;
            log_msg("WARN", "BB line $ln ($type_raw $tag_id): sequence is not pure ACGT; skipped") if $ROWS_SKIPPED{non_acgt_seq} <= 20;
            next;
        }

        # Normalize lib_id: treat placeholders (HP/OP/CP/NA/NONE) as global ('')
        my $lib_norm = defined($lib_id) ? $lib_id : "";
        if ( $lib_norm =~ /^(?:HP|OP|CP|NA|NONE)$/i ) {
            $lib_norm = "";
        }

        if ($type =~ /^CODON/){
            die "Codon row missing lib_id at line $ln (type=Codon)\n" if $lib_norm eq "";
            my $info = join($INFO_SEP, $type_raw,$tag_id,$bb_id_fixed,$smiles,"perf",$cycle);
            $CODON_TRIE->insert($seq,$info,$lib_norm);
            $LIB_IDS{$lib_norm}=1 if $lib_norm ne "";
            if ($lib_norm ne "" && $cycle =~ /^[1-4]$/){
                $LIB_EXPECTED_CYC{$lib_norm}{$cycle}=1;
            } elsif ($lib_norm ne "") {
                $ROWS_SKIPPED{codon_cycle_not_1_4}++;
                log_msg("WARN", "BB line $ln: Codon cycle '$cycle' is not 1..4; this codon can never be selected") if $ROWS_SKIPPED{codon_cycle_not_1_4} <= 20;
            }
        } elsif ($type eq "CP" || $type eq "HP" || $type eq "OP"){
            die "CP row missing lib_id at line $ln (type=CP)\n" if $type eq "CP" && $lib_norm eq "";
            my $info_p = join($INFO_SEP, $type,$tag_id,$bb_id_fixed,$smiles,"perf");
            my $info_m = join($INFO_SEP, $type,$tag_id,$bb_id_fixed,$smiles,"miss");
            if ($type eq "CP") {
                # Enforce design invariant: CP sequences must not collide across libraries.
                if (exists $CP_VARIANT_OWNER{$seq} && $CP_VARIANT_OWNER{$seq} ne $lib_norm) {
                    die "CP collision across libs (seq=$seq): $CP_VARIANT_OWNER{$seq} vs $lib_norm\n";
                }
                $CP_VARIANT_OWNER{$seq} = $lib_norm;
                $LIB_IDS{$lib_norm}=1 if $lib_norm ne "";
            }
            # pass 1: perfect sequence now; variants deferred to pass 2
            _insert_perfect_checked($tries{$type}, $type, $seq, $info_p, $lib_norm, $tag_id, $ln);
            push @anchor_rows, { type=>$type, seq=>$seq, info_m=>$info_m, lib=>$lib_norm, tag_id=>$tag_id, ln=>$ln }
                if mismatch_enabled_for_type($type);
        } else {
            $ROWS_SKIPPED{unknown_type}++;
            log_msg("WARN", "BB line $ln: unknown type '$type_raw' (expected CODON/HP/OP/CP); skipped") if $ROWS_SKIPPED{unknown_type} <= 20;
        }
    }
    close($IN);

    # pass 2: 1-bp variants (sub/del/ins) of every anchor, in file order; collisions are resolved order-independently
    for my $r (@anchor_rows) {
        my $type = $r->{type};
        my @mis = remove_dups(one_bp_substitution($r->{seq}), one_bp_deletion($r->{seq}), one_bp_insertion($r->{seq}));
        if ($type eq 'CP') {
            # Enforce design invariant: CP variants (1bp sub/indel) must never collide across libraries.
            for my $s (@mis) {
                if (exists $CP_VARIANT_OWNER{$s} && $CP_VARIANT_OWNER{$s} ne $r->{lib}) {
                    die "CP collision across libs (variant=$s): $CP_VARIANT_OWNER{$s} vs $r->{lib}\n";
                }
                $CP_VARIANT_OWNER{$s} = $r->{lib};
            }
        }
        for my $s (@mis){ _insert_variant_safe($tries{$type}, $type, $s, $r->{info_m}, $r->{lib}, $r->{tag_id}); }
    }

    # sanity: something must be decodable
    if (%ROWS_SKIPPED) {
        log_msg("WARN", "BB rows skipped: ".join(", ", map { "$_=$ROWS_SKIPPED{$_}" } sort keys %ROWS_SKIPPED));
    }
    die "No CP rows loaded from $file (check the 'type' column: CODON/HP/OP/CP)\n" unless $ANCHOR_ROWS{CP};
    for my $t (qw(HP OP)) { log_msg("WARN", "No $t rows loaded from $file; no read can be decoded") unless $ANCHOR_ROWS{$t}; }

    # lib별 expected cycle 확정 (Codon 기반)
    my %final;
    for my $lib (keys %LIB_EXPECTED_CYC){
        my @cyc = sort {$a<=>$b} keys %{ $LIB_EXPECTED_CYC{$lib} };
        my $maxc = @cyc ? $cyc[-1] : 3;
        $final{$lib} = ($maxc >= 4) ? 4 : 3;
    }
    for my $lib (keys %LIB_IDS){ $final{$lib} //= 3; }  # CP만 있는 경우 3
    %LIB_EXPECTED_CYC = %final;

    my @cycle_str = map { "$_:$LIB_EXPECTED_CYC{$_}" } sort keys %LIB_EXPECTED_CYC;
    my $skipped_str = %VARIANT_SKIPPED ? join(",", map { "$_:$VARIANT_SKIPPED{$_}" } sort keys %VARIANT_SKIPPED) : "none";
    log_msg("INFO","Trie build complete. Libs=".scalar(keys %LIB_IDS)."; cycles=".join(", ", @cycle_str)
                   ."; anchors=".join(",", map { "$_:".($ANCHOR_ROWS{$_}//0) } qw(HP OP CP))."; variants_dropped=$skipped_str");
}

# --------- Decoding helpers ----------
sub parse_info {
    my ($info)=@_;
    my @f = split /\Q$INFO_SEP\E/, $info;
    my %h = ( type=>$f[0], tag_id=>$f[1], bb_id=>$f[2], smiles=>$f[3], qual=>$f[4] );
    $h{cycle} = $f[5] if defined $f[5];
    return \%h;
}
sub pick_best_cp_candidates {
    my ($cps_ref,$max_n)=@_;
    my @cps=@$cps_ref;
    for my $m (@cps){ $m->{meta}=parse_info($m->{info}); $m->{match_len}=length($m->{sequence}); }
    my @sorted = sort {
        my $qa = ($a->{meta}{qual} eq 'perf') ? 1 : 0;
        my $qb = ($b->{meta}{qual} eq 'perf') ? 1 : 0;
        return $qb <=> $qa if $qa != $qb;
        return ($b->{position} <=> $a->{position});
    } @cps;
    if (defined $max_n && $max_n > 0 && @sorted > $max_n){
        splice(@sorted,$max_n);
    }
    return \@sorted;
}
sub pick_best_anchor_before {
    my ($matches_ref,$limit_pos,$max_n)=@_;
    my @cand = grep { $_->{position} < $limit_pos } @$matches_ref;
    return undef unless @cand;
    for my $m (@cand){ $m->{meta}=parse_info($m->{info}); }
    my @sorted = sort {
        my $qa = ($a->{meta}{qual} eq 'perf') ? 1 : 0;
        my $qb = ($b->{meta}{qual} eq 'perf') ? 1 : 0;
        return $qb <=> $qa if $qa != $qb;
        return $b->{position} <=> $a->{position};
    } @cand;
    if (defined $max_n && $max_n > 0 && @sorted > $max_n){
        splice(@sorted,$max_n);
    }
    return \@sorted;
}
sub _gap_ok {
    my ($next_start,$prev_end,$tol)=@_;
    my $gap = $next_start - ($prev_end + 1);
    return (abs($gap) <= $tol);
}
sub _pick_codons_dfs {
    my ($by_cycle,$cycle,$expected,$prev_end,$cp_pos,$tol,$memo)=@_;
    my $key = join(":", $cycle, $prev_end);
    return $memo->{$key} if exists $memo->{$key};
    if ($cycle > $expected){
        return $memo->{$key} = { penalty => 0, sel => {} };
    }

    my $best;
    for my $cand (@{ $by_cycle->{$cycle} || [] }){
        next unless _gap_ok($cand->{position}, $prev_end, $tol);
        if ($cycle == $expected){
            my $gap_cp = $cp_pos - ($cand->{end} + 1);
            next if abs($gap_cp) > $tol;
            my $pen = abs($cand->{position} - ($prev_end + 1)) + abs($gap_cp);
            if (!$best || $pen < $best->{penalty}){
                $best = { penalty => $pen, sel => { $cycle => $cand } };
            }
            next;
        }
        my $next = _pick_codons_dfs($by_cycle, $cycle+1, $expected, $cand->{end}, $cp_pos, $tol, $memo);
        next unless $next;
        my $pen = abs($cand->{position} - ($prev_end + 1)) + $next->{penalty};
        if (!$best || $pen < $best->{penalty}){
            my %sel = %{ $next->{sel} };
            $sel{$cycle} = $cand;
            $best = { penalty => $pen, sel => \%sel };
        }
    }
    return $memo->{$key} = $best;
}
sub pick_codons_by_adjacency {
    my ($lib_id,$seq,$op_end,$cp_pos,$expected_cycles,$adj_tol)=@_;
    my $codon_hits = $CODON_TRIE->search_substrings_by_lib_id($seq,$lib_id);
    return undef unless @$codon_hits;

    my %by_cycle;
    for my $h (@$codon_hits){
        my $pos=$h->{position};
        my $end=$pos + length($h->{sequence}) - 1;
        next if $pos < $op_end - $adj_tol;
        next if $end > $cp_pos + $adj_tol;
        my $meta = parse_info($h->{info});
        next unless defined $meta->{cycle} && $meta->{cycle} =~ /^[1-4]$/;
        my $cycle = $meta->{cycle}+0;
        $h->{meta} = $meta;
        $h->{end}  = $end;
        push @{ $by_cycle{$cycle} }, $h;
    }
    for my $cycle (1..$expected_cycles){
        return undef unless $by_cycle{$cycle} && @{$by_cycle{$cycle}};
        @{ $by_cycle{$cycle} } = sort { $a->{position} <=> $b->{position} } @{ $by_cycle{$cycle} };
    }

    my %memo;
    my $best = _pick_codons_dfs(\%by_cycle, 1, $expected_cycles, $op_end, $cp_pos, $adj_tol, \%memo);
    return undef unless $best && $best->{sel};

    my %out;
    for my $c (1..$expected_cycles){
        my $h=$best->{sel}{$c};
        return undef unless $h;
        $out{$c} = { bb_id=>$h->{meta}{bb_id}, pos=>$h->{position}, seq=>$h->{sequence} };
    }
    return \%out;
}
sub try_decode_direction {
    my ($seq,$direction)=@_;
    my $cp_hits = $CP_TRIE->search_substrings($seq);
    return undef unless @$cp_hits;
    my $cp_cands = pick_best_cp_candidates($cp_hits, $MAX_CP_CANDS);
    CP: for my $cp (@$cp_cands){
        my $lib_id = $cp->{lib_id};
        next CP unless exists $LIB_IDS{$lib_id};
        my $cp_pos  = $cp->{position};

        # Use lib-specific OP hits (and allow global anchors if present)
        my $op_hits = $OP_TRIE->search_substrings_by_lib_id($seq, $lib_id);
        if ($lib_id ne ''){
            my $op_glob = $OP_TRIE->search_substrings_by_lib_id($seq, '');
            push @$op_hits, @$op_glob if $op_glob && @$op_glob;
        }
        next CP unless @$op_hits;
        my $op_cands = pick_best_anchor_before($op_hits,$cp_pos,$MAX_ANCHOR_CANDS);
        next CP unless $op_cands && @$op_cands;

        # HP hits do not depend on the OP candidate: search once per CP candidate (was inside the OP loop)
        # Use lib-specific HP hits (and allow global anchors if present)
        my $hp_hits = $HP_TRIE->search_substrings_by_lib_id($seq, $lib_id);
        if ($lib_id ne ''){
            my $hp_glob = $HP_TRIE->search_substrings_by_lib_id($seq, '');
            push @$hp_hits, @$hp_glob if $hp_glob && @$hp_glob;
        }
        next CP unless @$hp_hits;
        my $expected_cycles = $LIB_EXPECTED_CYC{$lib_id} // 3;

        OP: for my $op (@$op_cands){
            my $op_pos=$op->{position};
            my $op_end = $op_pos + length($op->{sequence}) - 1;
            my $hp_cands = pick_best_anchor_before($hp_hits,$op_pos,$MAX_ANCHOR_CANDS);
            next OP unless $hp_cands && @$hp_cands;

            # Codon selection depends only on (lib, OP end, CP pos): compute once per OP candidate (was inside the HP loop)
            my $codon_sel = pick_codons_by_adjacency($lib_id,$seq,$op_end,$cp_pos,$expected_cycles,$ADJ_TOL);
            next OP unless $codon_sel;

            HP: for my $hp (@$hp_cands){
                my $hp_pos=$hp->{position};
                my $hp_end = $hp_pos + length($hp->{sequence}) - 1;
                next HP unless _gap_ok($op_pos, $hp_end, $ADJ_TOL);

                my @bbs; for my $c (1..$expected_cycles){ push @bbs, $codon_sel->{$c}{bb_id}; }
                push @bbs, "NA" if $expected_cycles==3;

                return {
                    lib_id => $lib_id, cycles=>$expected_cycles, bb_list=>\@bbs, direction=>$direction,
                    anchors => {
                        HP => { pos=>$hp_pos, seq=>$hp->{sequence}, qual=>(parse_info($hp->{info})->{qual}) },
                        OP => { pos=>$op->{position}, seq=>$op->{sequence}, qual=>(parse_info($op->{info})->{qual}) },
                        CP => { pos=>$cp->{position}, seq=>$cp->{sequence}, qual=>(parse_info($cp->{info})->{qual}) },
                    },
                    codons => $codon_sel,
                };
            }
        }
    }
    return undef;
}
sub decode_one_read {
    my ($seq)=@_;
    my $res = try_decode_direction($seq,'forward'); return $res if $res;
    $res = try_decode_direction(rc($seq),'revcomp'); return $res;
}

# --------- Counting / Matrices ----------
my %COUNTS_BY_LIB_ID;  # {lib}{id}{sample} = count
my %LIB_TOTALS;        # {lib}{sample}     = decoded total
my %SAMPLES;           # set of samples
my %SUMMARY;           # {sample}{metric}
my %ANCHOR_QUAL;       # {sample}{lib}{HP|OP|CP}{perf|miss|unknown}
my %DIR_COUNTS;        # {sample}{lib}{forward|revcomp}
my %CYCLE_DIST;        # {sample}{lib}{3|4}
my %BB_FREQ;           # {sample}{lib}{cycle}{bb_id}
my %GAP_STATS;         # {sample}{lib}{metric}->{n,sum,sumsq,min,max}
my %LEN_WIN;           # {sample}{within_3/within_4/within_both/within_multi/outside}
my %LEN_CONF;          # {sample}{... config ...}
my %QC;                # {sample}{codon_len_not9,order_violation,cp_len_out_of_27_29}

sub chain_id_str { my ($cycles,$bbs_ref)=@_; my @bb=@$bbs_ref; join("_",$cycles,@bb); }
sub add_count_matrix {
    my ($sample,$lib,$cycles,$bbs_ref)=@_;
    my $id=chain_id_str($cycles,$bbs_ref);
    $COUNTS_BY_LIB_ID{$lib}{$id}{$sample}++;
    $LIB_TOTALS{$lib}{$sample}++;
}
sub acc_qual  { my ($s,$l,$tag,$q)=@_; $q = (defined $q && length $q) ? $q : 'unknown'; $ANCHOR_QUAL{$s}{$l}{$tag}{$q}++; }
sub acc_dir   { my ($s,$l,$d)=@_;       $DIR_COUNTS{$s}{$l}{$d}++; }
sub acc_cycle { my ($s,$l,$c)=@_;       $CYCLE_DIST{$s}{$l}{$c}++; }
sub acc_bb    { my ($s,$l,$cy,$bb)=@_;  $BB_FREQ{$s}{$l}{$cy}{$bb}++; }
sub acc_gap {
    my ($s,$l,$m,$v)=@_;
    my $h = ($GAP_STATS{$s}{$l}{$m} //= { n=>0,sum=>0,sumsq=>0,min=>undef,max=>undef });
    $h->{n}++; $h->{sum}+=$v; $h->{sumsq}+=$v*$v;
    $h->{min} = defined $h->{min} ? ($v<$h->{min}?$v:$h->{min}) : $v;
    $h->{max} = defined $h->{max} ? ($v>$h->{max}?$v:$h->{max}) : $v;
}

# --------- 길이 윈도우/파일/JSON 유틸 ----------
my $AUTO_FASTQ_REGEX = qr/(?:\.fpmerged\.fq(?:\.gz)?|_merged\.(?:fq|fastq)(?:\.gz)?)$/i;
sub _compile_user_regex { my ($s)=@_; return undef unless defined $s && $s ne ''; my $re; eval { $re = qr/$s/; 1 } or die "Invalid --fastq-regex: $s\n$@\n"; return $re; }
my $FASTQ_REGEX = _compile_user_regex($OPT{fastq_regex}) // $AUTO_FASTQ_REGEX;

sub list_merged_fastqs {
    my ($dir)=@_;
    opendir(my $D,$dir) or die "Cannot open $dir: $!";
    my @all = readdir($D);
    closedir($D);
    my @cand = grep { $_ =~ $FASTQ_REGEX } @all;

    # Deduplicate per sample: prefer uncompressed over .gz, and .fastq over .fq
    my %pick; # sample => { fn => ..., rank => ... }
    for my $fn (@cand){
        my $sample = sample_name_from_filename($fn, $dir);
        my $path   = File::Spec->catfile($dir, $fn);
        my $is_gz     = ($fn =~ /\.gz$/i) ? 1 : 0;
        my $is_fastq  = ($fn =~ /\.fastq(?:\.gz)?\z/i) ? 1 : 0;
        # Lower rank is preferred
        my $rank = ($is_gz ? 10 : 0) + ($is_fastq ? 0 : 1);
        my $mtime = (stat($path))[9] // 0;
        if (!exists $pick{$sample} || $rank < $pick{$sample}{rank} ||
            ($rank == $pick{$sample}{rank} && $mtime > ($pick{$sample}{mtime} // 0))){
            $pick{$sample} = { fn => $fn, rank => $rank, mtime => $mtime };
        }
    }
    my @sel = sort map { $pick{$_}{fn} } keys %pick;
    return map { File::Spec->catfile($dir,$_) } @sel;
}
my %SAMPLE_NAME_CACHE;   # "$dir\0$fn" -> sample (avoids repeating the JSON lookup and its WARN for the same file)
sub sample_name_from_filename {
    my ($fn, $dir)=@_;
    my $ck = ($dir // '')."\0".$fn;
    return $SAMPLE_NAME_CACHE{$ck} if exists $SAMPLE_NAME_CACHE{$ck};
    my $s = $fn;
    $s =~ s/\.gz$//i;
    $s =~ s/\.(?:fastq|fq)$//i;

    my $base = $s;
    my $stripped = $s;
    $stripped =~ s/(?:\.fp)?merged$//i;
    $stripped =~ s/_merged$//i;
    $stripped =~ s/[._-]$//;

    if ($stripped ne $base){
        my $json_stripped = ($dir && find_fastp_json_for_sample($dir, $stripped)) ? 1 : 0;
        my $json_base     = ($dir && find_fastp_json_for_sample($dir, $base)) ? 1 : 0;
        if ($json_stripped && !$json_base){
            $base = $stripped;
        } elsif ($json_stripped && $json_base){
            log_msg("WARN","Ambiguous sample name for $fn (fastp JSON matches '$base' and '$stripped'); keeping '$base'");
        } elsif (!$json_stripped && !$json_base){
            # The file-name regex only admits merge suffixes (_merged / .fpmerged), so stripping them is unambiguous;
            # this keeps '<SAMPLE>' as the column name in the --skip-fastp path (README: FASTQ naming), as documented.
            log_msg("INFO","No fastp JSON for $fn; using sample name '$stripped' (merge suffix stripped)");
            $base = $stripped;
        }
    }
    $base =~ s/[._-]$//;
    $SAMPLE_NAME_CACHE{$ck} = $base;
    return $base;
}
sub open_maybe_gzip {
    my ($path)=@_;
    if ($path =~ /\.gz$/i){
        # Prefer system gzip -cd for robust handling of concatenated (multi-member) streams.
        if (open(my $FH, "-|", "gzip","-cd","--",$path)){
            log_msg("INFO","Using system gzip -cd for $path (multi-stream safe)");
            return $FH;
        } else {
            # Fallback to Perl IO::Uncompress with MultiStream enabled
            my $use_perl = eval { require IO::Uncompress::Gunzip; 1; };
            die "Cannot run gzip -cd on $path and IO::Uncompress::Gunzip not available: $!\n" unless $use_perl;
            my $z = IO::Uncompress::Gunzip->new(
                $path,
                MultiStream => 1,   # ensure all concatenated members are read
                Transparent => 0,
                Append      => 0,
                Strict      => 0,
                BinModeIn   => 1,
            );
            unless ($z) {
                no warnings 'once'; # GunzipError used only once 경고 억제
                die "Cannot gunzip $path (IO::Uncompress::Gunzip failed: $IO::Uncompress::Gunzip::GunzipError)\n";
            }
            log_msg("INFO","Using IO::Uncompress::Gunzip(MultiStream=1) for $path");
            return $z;
        }
    } else {
        open(my $FH, "<", $path) or die "Cannot open $path: $!\n";
        binmode($FH);
        return $FH;
    }
}
sub close_input_or_die {
    # A truncated/corrupt .gz makes `gzip -cd` exit non-zero while Perl only sees EOF; without this check the sample
    # would be counted as complete (exit 0, .decode_done written) with a silently reduced read count.
    my ($fh, $path, $sample) = @_;
    if (ref($fh) && eval { $fh->isa('IO::Uncompress::Gunzip') }) {
        my $err = $IO::Uncompress::Gunzip::GunzipError;
        $fh->close();
        die "Decompression of $path ($sample) failed: $err. Input may be truncated; aborting to avoid partial counts.\n"
            if defined $err && length $err;
        return;
    }
    local $! = 0;
    unless (close($fh)) {
        my $status = $?;
        my $rc  = $status >> 8;
        my $sig = $status & 127;
        die sprintf("Input stream for %s (%s) ended abnormally (decompressor exit=%d%s%s). Input may be truncated; aborting to avoid partial counts.\n",
                    $sample, $path, $rc, ($sig ? ", signal=$sig" : ""), ($! ? ", err=$!" : ""))
            if $status != 0 || $!;
    }
}
sub find_fastp_json_for_sample {
    my ($dir,$sample)=@_;
    my $p1 = File::Spec->catfile($dir, "$sample.fastp.json");
    return $p1 if -e $p1;
    opendir(my $D,$dir) or return undef;
    my @c = grep { /\.fastp\.json$/ && /^\Q$sample\E[._-]/ } readdir($D);
    closedir($D);
    return @c ? File::Spec->catfile($dir,$c[0]) : undef;
}
sub _slurp_text { my ($path)=@_; open(my $F,"<",$path) or die "Cannot open $path: $!\n"; local $/; my $txt=<$F>; close($F); return $txt; }

# histogram normalization
sub _normalize_hist_to_hash {
    my ($h) = @_;
    my %out; return %out unless defined $h;
    if (ref($h) eq 'HASH') {
        my $ok_pairs=0;
        for my $k (keys %$h){
            if (looks_like_number($k) && looks_like_number($h->{$k})) { $out{int($k)} += $h->{$k}; $ok_pairs++ }
        }
        return %out if $ok_pairs;
        if (ref($h->{x}) eq 'ARRAY' && ref($h->{y}) eq 'ARRAY'){
            my $n = @{$h->{x}} < @{$h->{y}} ? @{$h->{x}} : @{$h->{y}};
            for (my $i=0; $i<$n; $i++){
                my ($x,$y)=($h->{x}[$i], $h->{y}[$i]);
                next unless defined $x && defined $y;
                next unless looks_like_number($x) && looks_like_number($y);
                $out{int($x)} += $y+0;
            } return %out;
        }
    } elsif (ref($h) eq 'ARRAY') {
        my $are_pairs = (ref($h->[0]) eq 'ARRAY' && @{$h->[0]}>=2);
        if ($are_pairs){
            for my $row (@$h){
                my ($x,$y) = @$row[0,1];
                next unless defined $x && defined $y;
                next unless looks_like_number($x) && looks_like_number($y);
                $out{int($x)} += $y+0;
            }
        } else {
            for (my $i=0; $i<@$h; $i++){
                my $y=$h->[$i]; next unless defined $y && looks_like_number($y);
                $out{$i} += $y+0;
            }
        } return %out;
    }
    return %out;
}
sub _sum_range {
    my ($href,$lo,$hi)=@_; my $sum=0; for my $k (keys %$href){ next if $k < $lo || $k > $hi; $sum += $href->{$k}; } return $sum;
}
sub _argmax_key {
    my ($href)=@_; my ($best_k,$best_v)=(undef,-1e99);
    for my $k (sort { $a <=> $b } keys %$href){ my $v=$href->{$k}; next unless defined $v; if ($v>$best_v){ $best_v=$v; $best_k=$k; } }   # ties -> smallest length (deterministic)
    return $best_k;
}
sub _top_k_peaks {
    my ($href,$k,$min_sep,$min_frac,$total) = @_;
    my @pairs = sort { $href->{$b} <=> $href->{$a} || $a <=> $b } keys %$href; # by count desc, ties by length asc (deterministic across runs)
    my (@sel,@sel_counts);
    CAND: for my $x (@pairs){
        my $cnt = $href->{$x}+0;
        my $frac = $total>0 ? $cnt/$total : 0;
        next if $frac < $min_frac;
        for my $c (@sel){
            next CAND if abs($x - $c) < $min_sep;
        }
        push @sel, $x; push @sel_counts, $cnt;
        last if @sel >= $k;
    }
    return (\@sel,\@sel_counts);
}

# fastp에서 insert size 피크/지지도 감지 (Top-K)
sub detect_insert_size_for_sample {
    my ($dir,$sample,$exp3_base,$exp4_base,$tol_base,$k,$min_sep,$min_frac)=@_;
    my $json_path = find_fastp_json_for_sample($dir,$sample) or return undef;
    my $txt = eval { _slurp_text($json_path) }; if ($@){ log_msg("WARN","Cannot read fastp json for $sample: $@"); return undef; }
    my $J = eval { JSON::PP->new->utf8->decode($txt) }; if ($@){ log_msg("WARN","JSON parse failed for $json_path: $@"); return undef; }

    my ($peak, $width, %hist);
    if (exists $J->{insert_size} && ref($J->{insert_size}) eq 'HASH'){
        my $ins = $J->{insert_size};
        $peak  = $ins->{peak} // $ins->{mode} // $ins->{peak_insert_size};
        $width = $ins->{width};
        %hist  = _normalize_hist_to_hash($ins->{histogram} // $ins->{hist} // $ins->{distribution});
    } else {
        for my $k2 (qw(insertSize insert-size)){
            next unless exists $J->{$k2} && ref($J->{$k2}) eq 'HASH';
            my $ins = $J->{$k2};
            $peak  //= $ins->{peak} // $ins->{mode};
            $width //= $ins->{width};
            my %h2   = _normalize_hist_to_hash($ins->{histogram} // $ins->{hist} // $ins->{distribution});
            %hist = (%hist, %h2) if %h2;
        }
    }
    if (!defined $peak && %hist){ $peak = _argmax_key(\%hist); }

    my ($support3,$support4,$total_hist)=(undef,undef,0);
    if (%hist){
        $total_hist += $_ for values %hist;
        $support3 = _sum_range(\%hist, $exp3_base-$tol_base, $exp3_base+$tol_base);
        $support4 = _sum_range(\%hist, $exp4_base-$tol_base, $exp4_base+$tol_base);
    }
    my ($peaks,$counts) = (%hist ? _top_k_peaks(\%hist,$k,$min_sep,$min_frac,$total_hist) : ([],[]));

    return {
        peak        => $peak,
        peaks       => $peaks,        # ARRAY ref (Top-K)
        peak_counts => $counts,       # ARRAY ref
        width       => $width,
        hist_n      => scalar(keys %hist),
        total       => $total_hist,
        support3    => $support3,
        support4    => $support4,
        json_path   => $json_path,
    };
}

# --------- 길이 윈도우 판정 ----------
my $ACTIVE_MODE = 'both';        # '3' | '4' | 'both' | 'multi'
my @ACTIVE_CENTERS = ();         # multi 모드에서 사용하는 길이 중심들(Top-K)
sub classify_len {
    my ($L) = @_;
    my $in3 = (abs($L - $EXP_LEN_3) <= $LENGTH_TOL) ? 1 : 0;
    my $in4 = (abs($L - $EXP_LEN_4) <= $LENGTH_TOL) ? 1 : 0;

    if    ($ACTIVE_MODE eq '3')   { $in4 = 0; }
    elsif ($ACTIVE_MODE eq '4')   { $in3 = 0; }

    if ($ACTIVE_MODE ne 'multi'){
        my $cls = $in3 && $in4 ? 'within_both' : ($in3 ? 'within_3' : ($in4 ? 'within_4' : 'outside'));
        return ($cls, $in3, $in4);
    } else {
        my $in_multi = 0;
        for my $c (@ACTIVE_CENTERS){ if (abs($L - $c) <= $LENGTH_TOL){ $in_multi=1; last; } }
        my $cls;
        if ($in_multi){
            $cls = $in3 && $in4 ? 'within_both'
                : ($in3 ? 'within_3'
                : ($in4 ? 'within_4' : 'within_multi'));
        } else {
            $cls = 'outside';
        }
        return ($cls, $in3, $in4);
    }
}
sub len_filter_reject {
    my ($cls) = @_;
    if ($LENGTH_FILTER_MODE eq 'reject_outside') { return ($cls eq 'outside') ? 1 : 0; }
    if ($LENGTH_FILTER_MODE eq 'reject_inside')  { return ($cls =~ /^within_/) ? 1 : 0; }
    return ($cls eq 'outside') ? 1 : 0;
}

# --------- merged fastq 처리 ----------
sub process_merged_dir {
    my ($dir)=@_;
    my @fq_paths = list_merged_fastqs($dir);
    die "No merged fastq found in $dir (regex=".($OPT{fastq_regex}//'AUTO').")\n"
        unless @fq_paths;

    log_msg("INFO","Length policy = $LENGTH_FILTER_MODE; TOL=$LENGTH_TOL_BASE; base exp3=$EXP_LEN_3_BASE; base exp4=$EXP_LEN_4_BASE");
    log_msg("INFO","Scaling denominator policy = $SCALING_DENOM");
    log_msg("INFO","FASTQ match regex = ".($OPT{fastq_regex}//'AUTO')."; files=".scalar(@fq_paths));
    for my $p (@fq_paths){ log_msg("INFO","  found: ".basename($p)); }
    dump_effective_config();

    for my $path (@fq_paths){
        my $fn = basename($path);
        my $sample = sample_name_from_filename($fn, $dir);
        $SAMPLES{$sample}=1;

        # --- fastp 기반 자동 길이 감지/조정 (Top-K => multi) ---
        # NOTE: start from the CLI base values every sample; $LENGTH_TOL is overwritten per sample below (-U 1)
        my ($exp3_local,$exp4_local,$tol_local,$mode_local) = ($EXP_LEN_3_BASE,$EXP_LEN_4_BASE,$LENGTH_TOL_BASE,'both');
        my (@centers_local, @peaks_local, @peak_counts_local);
        if ($AUTO_INSERT_DETECT){
            my $det = detect_insert_size_for_sample($dir,$sample,$EXP_LEN_3_BASE,$EXP_LEN_4_BASE,$LENGTH_TOL_BASE,
                                                    $PEAK_TOP_K,$PEAK_MIN_SEP,$PEAK_MIN_FRAC);
            if ($det && ref($det->{peaks}) eq 'ARRAY' && @{$det->{peaks}}){
                @peaks_local       = @{$det->{peaks}};
                @peak_counts_local = @{$det->{peak_counts} // []};
                $mode_local = 'multi';
                @centers_local = @peaks_local;

                if (defined $det->{support3} && defined $det->{support4} && $det->{total}){
                    my $r3 = $det->{support3} / $det->{total};
                    my $r4 = $det->{support4} / $det->{total};
                    log_msg("INFO", sprintf("fastp[%s] peaks(top%d)=%s; counts=%s; support3=%.3f, support4=%.3f",
                        $det->{json_path}, scalar(@peaks_local),
                        join(",",@peaks_local), ( @peak_counts_local ? join(",",@peak_counts_local) : "NA"),
                        $r3, $r4));
                } else {
                    log_msg("INFO", sprintf("fastp[%s] peaks(top%d)=%s; counts=%s",
                        $det->{json_path}, scalar(@peaks_local),
                        join(",",@peaks_local), ( @peak_counts_local ? join(",",@peak_counts_local) : "NA")));
                }

                # 선택/센터링 로직(원하면 멀티 대신 단일 선택)
                # NOTE: previously `$AUTO_SELECT_WINDOW && !$AUTO_CENTER_WINDOW` combined with the centering
                # guard `$mode_local ne 'multi'` made both branches unreachable under the shipped defaults.
                if ($AUTO_SELECT_WINDOW){
                    if (defined $det->{support3} && defined $det->{support4} && $det->{total}){
                        my $r3 = $det->{support3} / $det->{total};
                        my $r4 = $det->{support4} / $det->{total};
                        if ($r3 <= 0 && $r4 <= 0) {
                            # no fastp reads near either design length: a 3/4 choice would be arbitrary -> keep Top-K 'multi'
                            log_msg("WARN", sprintf("%s: -W 1 requested but fastp shows no support for the 3- or 4-cycle window (exp3=%d, exp4=%d, tol=%d); keeping 'multi'",
                                                    $sample, $EXP_LEN_3_BASE, $EXP_LEN_4_BASE, $LENGTH_TOL_BASE));
                        } else {
                            if ($r3 >= 0.20 && $r4 >= 0.20) { $mode_local = 'both'; }
                            else { $mode_local = ($r3 >= $r4) ? '3' : '4'; }
                            @centers_local = ();   # centres are only meaningful in 'multi' mode
                        }
                    }
                }
                if ($AUTO_CENTER_WINDOW){
                    if ($mode_local ne 'multi'){
                        my ($target, $nearest_peak);
                        if     ($mode_local eq '3'){ $target = $EXP_LEN_3_BASE; }
                        elsif  ($mode_local eq '4'){ $target = $EXP_LEN_4_BASE; }
                        if (defined $target){
                            my $best_d=1e9;
                            for my $p (@peaks_local){ my $d=abs($p-$target); if ($d<$best_d){$best_d=$d;$nearest_peak=$p;} }
                            # Only re-centre within a bounded shift; a far-away Top-K peak (adapter dimer, off-design
                            # product) must not drag the design window away from every real read.
                            my $max_shift = 2*$LENGTH_TOL_BASE; $max_shift = 6 if $max_shift < 6;
                            if (defined $nearest_peak && $best_d <= $max_shift){
                                if ($mode_local eq '3'){ $exp3_local = $nearest_peak; }
                                if ($mode_local eq '4'){ $exp4_local = $nearest_peak; }
                            } elsif (defined $nearest_peak) {
                                log_msg("WARN", sprintf("%s: nearest fastp peak %d is %d nt from design length %d (> max shift %d); window not re-centred",
                                                        $sample, $nearest_peak, $best_d, $target, $max_shift));
                            }
                        }
                    }
                }
                if ($AUTO_TOL_FROM_FASTP && defined $det->{width} && looks_like_number($det->{width})){
                    my $w = $det->{width}+0; my $tol = int($w/2 + 0.5); $tol = 2 if $tol < 2; $tol = 20 if $tol > 20;
                    $tol_local = $tol;
                }
            } else {
                log_msg("WARN","No usable insert-size information for sample=$sample; using base windows.");
            }
        }

        ($EXP_LEN_3,$EXP_LEN_4,$LENGTH_TOL,$ACTIVE_MODE) = ($exp3_local,$exp4_local,$tol_local,$mode_local);
        @ACTIVE_CENTERS = @centers_local;
        $LEN_CONF{$sample} = {
            policy        => $LENGTH_FILTER_MODE,
            tol           => $LENGTH_TOL,
            exp3          => $EXP_LEN_3,
            exp4          => $EXP_LEN_4,
            active_mode   => $ACTIVE_MODE,
            fastp_peak    => (@peaks_local ? $peaks_local[0] : undef),
            fastp_peaks   => [ @peaks_local ],
            active_centers=> [ @ACTIVE_CENTERS ],
            fastp_json    => find_fastp_json_for_sample($dir,$sample),
        };

        my $mode_disp = ($ACTIVE_MODE eq 'multi') ? "multi(".scalar(@ACTIVE_CENTERS)." centers: ".join(",",@ACTIVE_CENTERS).")" : $ACTIVE_MODE;
        log_msg("INFO", sprintf("Decoding sample: %s (%s) | mode=%s, exp3=%d, exp4=%d, tol=%d",
                                $sample,$path,$mode_disp,$EXP_LEN_3,$EXP_LEN_4,$LENGTH_TOL));

        my $dec_path = File::Spec->catfile($OUT_DIR,"decoded_reads_$sample.tsv");
        my $und_path = File::Spec->catfile($OUT_DIR,"undecoded_reads_$sample.tsv");
        open(my $FDEC, ">", $dec_path) or die "Cannot write $dec_path: $!\n";
        my @DEC_HEADER = build_decoded_header();
        print_tsv_header($FDEC,\@DEC_HEADER,'decoded');

        open(my $FUND, ">", $und_path) or die "Cannot write $und_path: $!\n";
        my @UND_HEADER = build_undecoded_header();
        print_tsv_header($FUND,\@UND_HEADER,'undecoded');

        my $und_written=0;
        my $IN = open_maybe_gzip($path);
        my ($total_reads,$decoded,$length_passed)=(0,0,0);
        my %fail=( no_cp=>0, no_op=>0, no_hp=>0, codon_fail=>0, len_out_of_range=>0, order_violation=>0 );

        while (1){
            my $hdr = <$IN>; last unless defined $hdr;
            my $seq = <$IN>; my $plus = <$IN>; my $qual = <$IN>;
            if (!defined $qual) {
                log_msg("WARN", "$sample: truncated FASTQ record after read $total_reads (header present, <4 lines); dropped");
                last;
            }
            chomp($hdr); chomp($seq);
            $hdr =~ s/\r\z//; $seq =~ s/\r\z//; $seq = uc($seq);   # tolerate CRLF (header + sequence) and lower-case FASTQ
            $total_reads++;
            my $read_id=$hdr; $read_id =~ s/^\@//;
            my $read_len=length($seq);

            # 길이 필터
            my ($cls,$in3,$in4) = classify_len($read_len);
            $LEN_WIN{$sample}{$cls}++;
            if ( len_filter_reject($cls) ){
                $fail{len_out_of_range}++;
                if ($und_written < $MAX_FAILED_DUMP){
                    my $range_str = sprintf("exp3=%d±%d[%d..%d],exp4=%d±%d[%d..%d],mode=%s",
                                             $EXP_LEN_3,$LENGTH_TOL,$EXP_LEN_3-$LENGTH_TOL,$EXP_LEN_3+$LENGTH_TOL,
                                             $EXP_LEN_4,$LENGTH_TOL,$EXP_LEN_4-$LENGTH_TOL,$EXP_LEN_4+$LENGTH_TOL,$ACTIVE_MODE);
                    my $reason = "len_out_of_range(len=$read_len,policy=$LENGTH_FILTER_MODE,$range_str)";
                    print $FUND join("\t",$read_id,$read_len,$seq,$reason,0,0,0), "\n";
                    $und_written++;
                }
                next;
            }

            $length_passed++;

            my $res = decode_one_read($seq);
            if ($res){
                my $hp_pos=$res->{anchors}{HP}{pos}; my $hp_len=length($res->{anchors}{HP}{seq});
                my $op_pos=$res->{anchors}{OP}{pos}; my $op_len=length($res->{anchors}{OP}{seq});
                my $cp_pos=$res->{anchors}{CP}{pos}; my $cp_len=length($res->{anchors}{CP}{seq});

                my (@cpos,@clen);
                for my $c (1..$res->{cycles}){
                    my $p=$res->{codons}{$c}{pos};
                    my $l=length($res->{codons}{$c}{seq});
                    push @cpos,$p; push @clen,$l;
                }

                # QC thresholds follow the canonical lengths (-c / -Q); column names kept for compatibility
                for my $l (@clen){ $QC{$sample}{codon_len_not9}++ if ($l != $CODON_CAN_LEN); }
                my $order_ok = 1;
                $order_ok &&= ($hp_pos <= $op_pos);
                my $prev = $op_pos;
                for my $p (@cpos){ $order_ok &&= ($prev <= $p); $prev = $p; }
                $order_ok &&= ($cpos[-1] <= $cp_pos);
                $QC{$sample}{order_violation}++ unless $order_ok;
                $QC{$sample}{cp_len_out_of_27_29}++ unless ($cp_len >= $CP_CAN_LEN-1 && $cp_len <= $CP_CAN_LEN+1);

                if (!$order_ok){
                    $fail{order_violation}++;
                    if ($und_written < $MAX_FAILED_DUMP){
                        print $FUND join("\t",$read_id,$read_len,$seq,"order_violation",1,1,1), "\n";
                        $und_written++;
                    }
                    next;
                }

                add_count_matrix($sample,$res->{lib_id},$res->{cycles},$res->{bb_list});
                $decoded++;
                acc_dir($sample,$res->{lib_id},$res->{direction});
                acc_cycle($sample,$res->{lib_id},$res->{cycles});
                for my $tag (qw(HP OP CP)){ acc_qual($sample,$res->{lib_id},$tag,$res->{anchors}{$tag}{qual}); }

                my $prev_end = $hp_pos+$hp_len-1; acc_gap($sample,$res->{lib_id},"gap_HP_OP",$op_pos-($prev_end+1));
                $prev_end = $op_pos+$op_len-1;

                for my $c (1..$res->{cycles}){
                    acc_bb($sample,$res->{lib_id},$c,$res->{codons}{$c}{bb_id});
                }

                acc_gap($sample,$res->{lib_id},"gap_OP_C1",$cpos[0]-($prev_end+1));
                for my $i (1..$#cpos){ my $pe=$cpos[$i-1]+$clen[$i-1]-1; acc_gap($sample,$res->{lib_id},"gap_C".$i."_C".($i+1), $cpos[$i]-($pe+1)); }
                my $last_c_end = $cpos[-1]+$clen[-1]-1;
                acc_gap($sample,$res->{lib_id},"gap_C_last_CP",$cp_pos-($last_c_end+1));

                my @bb=@{ $res->{bb_list} }; my $idstr=chain_id_str($res->{cycles},\@bb);
                my ($c1p,$c1l,$c1s,$c1b)=("NA","NA","NA","NA");
                my ($c2p,$c2l,$c2s,$c2b)=("NA","NA","NA","NA");
                my ($c3p,$c3l,$c3s,$c3b)=("NA","NA","NA","NA");
                my ($c4p,$c4l,$c4s,$c4b)=("NA","NA","NA","NA");
                if ($res->{cycles}>=1){ $c1p=$res->{codons}{1}{pos}; $c1s=$res->{codons}{1}{seq}; $c1l=length($c1s); $c1b=$res->{codons}{1}{bb_id}; }
                if ($res->{cycles}>=2){ $c2p=$res->{codons}{2}{pos}; $c2s=$res->{codons}{2}{seq}; $c2l=length($c2s); $c2b=$res->{codons}{2}{bb_id}; }
                if ($res->{cycles}>=3){ $c3p=$res->{codons}{3}{pos}; $c3s=$res->{codons}{3}{seq}; $c3l=length($c3s); $c3b=$res->{codons}{3}{bb_id}; }
                if ($res->{cycles}==4){ $c4p=$res->{codons}{4}{pos}; $c4s=$res->{codons}{4}{seq}; $c4l=length($c4s); $c4b=$res->{codons}{4}{bb_id}; }

                print $FDEC join("\t",
                    $read_id,$read_len,$seq,$res->{direction},$res->{lib_id},$res->{cycles},$idstr,
                    $res->{anchors}{HP}{pos},$hp_len,$res->{anchors}{HP}{qual},$res->{anchors}{HP}{seq},
                    $res->{anchors}{OP}{pos},$op_len,$res->{anchors}{OP}{qual},$res->{anchors}{OP}{seq},
                    $cp_pos,$cp_len,$res->{anchors}{CP}{qual},$res->{anchors}{CP}{seq},
                    $c1p,$c1l,$c1s,$c1b,
                    $c2p,$c2l,$c2s,$c2b,
                    $c3p,$c3l,$c3s,$c3b,
                    $c4p,$c4l,$c4s,$c4b
                ), "\n";

            } else {
                # Failure reason: probe the orientation in which a CP exists (decode_one_read tries forward then
                # revcomp), so revcomp reads are classified by their real missing element instead of always 'no_cp'.
                # Each trie is searched once; the flags are reused for the undecoded dump.
                my ($probe,$cp_found) = ($seq,0);
                for my $s ($seq, rc($seq)){
                    if (@{ $CP_TRIE->search_substrings($s) }){ $cp_found=1; $probe=$s; last; }
                }
                my $op_found = $cp_found ? (@{ $OP_TRIE->search_substrings($probe) } ? 1 : 0) : 0;
                my $hp_found = $op_found ? (@{ $HP_TRIE->search_substrings($probe) } ? 1 : 0) : 0;
                my $reason = (!$cp_found) ? "no_cp" : (!$op_found ? "no_op" : (!$hp_found ? "no_hp" : "codon_fail"));
                $fail{$reason}++;
                if ($und_written < $MAX_FAILED_DUMP){
                    print $FUND join("\t",$read_id,$read_len,$seq,$reason,$cp_found,$op_found,$hp_found), "\n";
                    $und_written++;
                }
            }

            if ($total_reads % 200000 == 0){
                log_msg("INFO","  $sample progress: $total_reads reads, decoded=$decoded");
            }
        }
        close_input_or_die($IN, $path, $sample);

        $SUMMARY{$sample}{total_reads}         = $total_reads;
        $SUMMARY{$sample}{length_passed_reads} = $length_passed;
        $SUMMARY{$sample}{decoded_reads}       = $decoded;
        $SUMMARY{$sample}{no_cp}               = $fail{no_cp};
        $SUMMARY{$sample}{no_op}               = $fail{no_op};
        $SUMMARY{$sample}{no_hp}               = $fail{no_hp};
        $SUMMARY{$sample}{codon_fail}          = $fail{codon_fail};
        $SUMMARY{$sample}{len_out_of_range}    = $fail{len_out_of_range};
        $SUMMARY{$sample}{order_violation}     = $fail{order_violation};

        log_msg("INFO","  $sample done: total=$total_reads, length_passed=$length_passed, decoded=$decoded, ".
                       "fails(no_cp=$fail{no_cp}, no_op=$fail{no_op}, no_hp=$fail{no_hp}, codon_fail=$fail{codon_fail}, len_out_of_range=$fail{len_out_of_range}, order_violation=$fail{order_violation})");
        close($FDEC); close($FUND);
    }
}

# --------- 헤더/출력 ----------
sub build_decoded_header {
    my @base    = qw(read_id read_len read_seq direction lib_id cycles id);
    my @anchors = map { ("${_}_pos","${_}_len","${_}_qual","${_}_seq") } qw(HP OP CP);
    my @codons  = map { ("C${_}_pos","C${_}_len","C${_}_seq","C${_}_bb_id") } (1..4);
    return (@base, @anchors, @codons);  # 35
}
sub build_undecoded_header { return qw(read_id read_len read_seq reason cp_found op_found hp_found); } # 7
sub print_tsv_header {
    my ($fh,$cols_ref,$label)=@_;
    my $line = join("\t", @$cols_ref);
    my $expected = ($label eq 'decoded') ? 35 : 7;
    my $n = scalar(@$cols_ref);
    # die() alone is enough: $SIG{__DIE__} logs the message once as [DIE]
    if ($n != $expected){ die "Header[$label] invalid: has $n cols (exp $expected).\n"; }
    if ($label eq 'decoded' && $line =~ /OP_seqCP_pos/){ die "Header concat detected.\n"; }
    print $fh $line, "\n";
    log_msg("INFO","Header[$label] columns = $n");
}

# --------- Matrix/통계 출력 ----------
sub _sort_ids_nicely { return sort { my ($ac,@ar)=split /_/,$a; my ($bc,@br)=split /_/,$b; $ac <=> $bc || $a cmp $b } @_; }

sub _denominator_for {
    my ($lib,$sample)=@_;
    if ($SCALING_DENOM eq 'sample_decoded'){
        return $SUMMARY{$sample}{decoded_reads}//0;
    } elsif ($SCALING_DENOM eq 'sample_length_passed'){
        return $SUMMARY{$sample}{length_passed_reads}//0;
    } elsif ($SCALING_DENOM eq 'sample_total_raw'){
        return $SUMMARY{$sample}{total_reads}//0;
    } elsif ($SCALING_DENOM eq 'library_decoded'){
        return $LIB_TOTALS{$lib}{$sample}//0;
    } else {
        return $SUMMARY{$sample}{decoded_reads}//0;
    }
}
sub write_raw_matrix {
    my ($f)=@_; my @samples=sort keys %SAMPLES;
    open(my $OUT,">",$f) or die "Cannot write $f: $!\n";
    print $OUT join("\t","lib_id","id",@samples), "\n";
    for my $lib (sort keys %COUNTS_BY_LIB_ID){
        for my $id (_sort_ids_nicely(keys %{ $COUNTS_BY_LIB_ID{$lib} })){
            print $OUT join("\t",$lib,$id, map { $COUNTS_BY_LIB_ID{$lib}{$id}{$_}//0 } @samples), "\n";
        }
    }
    close($OUT); log_msg("INFO","Wrote raw matrix: $f");
}
sub write_library_totals {
    my ($f)=@_; my @samples=sort keys %SAMPLES;
    open(my $OUT,">",$f) or die "Cannot write $f: $!\n";
    print $OUT join("\t","lib_id",@samples), "\n";
    for my $lib (sort keys %LIB_IDS){
        print $OUT join("\t",$lib, map { $LIB_TOTALS{$lib}{$_}//0 } @samples), "\n";
    }
    close($OUT); log_msg("INFO","Wrote library totals: $f");
}
sub write_scaled_matrix {
    my ($f)=@_; my @samples=sort keys %SAMPLES;
    open(my $OUT,">",$f) or die "Cannot write $f: $!\n";
    print $OUT join("\t","lib_id","id",@samples), "\n";
    for my $lib (sort keys %COUNTS_BY_LIB_ID){
        for my $id (_sort_ids_nicely(keys %{ $COUNTS_BY_LIB_ID{$lib} })){
            my @vals;
            for my $s (@samples){
                my $raw = $COUNTS_BY_LIB_ID{$lib}{$id}{$s}//0;
                my $den = _denominator_for($lib,$s);
                my $scaled = ($den>0) ? ($raw*1_000_000.0/$den) : 0;
                push @vals, sprintf("%.6f",$scaled);
            }
            print $OUT join("\t",$lib,$id,@vals), "\n";
        }
    }
    close($OUT); log_msg("INFO","Wrote scaled(CPM) matrix (denom_policy=$SCALING_DENOM): $f");
}
sub write_summary_basic {
    my ($f)=@_; open(my $OUT,">",$f) or die "Cannot write $f: $!\n";
    print $OUT join("\t",qw(sample total_reads decoded_reads no_cp no_op no_hp codon_fail len_out_of_range order_violation)), "\n";
    for my $s (sort keys %SUMMARY){
        my $h=$SUMMARY{$s};
        print $OUT join("\t",$s, map { $h->{$_}//0 } qw(total_reads decoded_reads no_cp no_op no_hp codon_fail len_out_of_range order_violation)), "\n";
    }
    close($OUT); log_msg("INFO","Wrote decoding summary: $f");
}
sub write_sample_stats {
    my ($f)=@_; open(my $OUT,">",$f) or die "Cannot write $f: $!\n";
    print $OUT join("\t",qw(sample total_reads length_passed_reads decoded_reads decode_rate_pct)), "\n";
    for my $s (sort keys %SUMMARY){
        my $tot=$SUMMARY{$s}{total_reads}//0; 
        my $lenp=$SUMMARY{$s}{length_passed_reads}//0;
        my $dec=$SUMMARY{$s}{decoded_reads}//0;
        my $rate = $tot ? sprintf("%.3f",100.0*$dec/$tot) : "0.000";
        print $OUT join("\t",$s,$tot,$lenp,$dec,$rate), "\n";
    }
    close($OUT); log_msg("INFO","Wrote sample_stats: $f");
}
sub write_lib_sample_stats {
    my ($f)=@_; my @samples=sort keys %SAMPLES;
    open(my $OUT,">",$f) or die "Cannot write $f: $!\n";
    print $OUT join("\t",qw(sample lib_id decoded_in_lib ratio_among_decoded_pct ratio_among_total_pct)), "\n";
    for my $s (@samples){
        my $dec_total=$SUMMARY{$s}{decoded_reads}//0;
        my $tot_reads=$SUMMARY{$s}{total_reads}//0;
        for my $l (sort keys %LIB_IDS){
            my $x=$LIB_TOTALS{$l}{$s}//0;
            my $r1 = $dec_total ? sprintf("%.3f", 100.0*$x/$dec_total) : "0.000";
            my $r2 = $tot_reads ? sprintf("%.3f", 100.0*$x/$tot_reads) : "0.000";
            print $OUT join("\t",$s,$l,$x,$r1,$r2), "\n";
        }
    }
    close($OUT); log_msg("INFO","Wrote lib_sample_stats: $f");
}
sub write_anchor_qual_counts {
    my ($f)=@_; open(my $OUT,">",$f) or die "Cannot write $f: $!\n";
    print $OUT join("\t",qw(sample lib_id HP_perf HP_miss OP_perf OP_miss CP_perf CP_miss forward_reads revcomp_reads)), "\n";
    for my $s (sort keys %ANCHOR_QUAL){
        for my $l (sort keys %{ $ANCHOR_QUAL{$s} }){
            next unless exists $LIB_IDS{$l};
            my $HPp=$ANCHOR_QUAL{$s}{$l}{HP}{perf}//0; my $HPm=$ANCHOR_QUAL{$s}{$l}{HP}{miss}//0;
            my $OPp=$ANCHOR_QUAL{$s}{$l}{OP}{perf}//0; my $OPm=$ANCHOR_QUAL{$s}{$l}{OP}{miss}//0;
            my $CPp=$ANCHOR_QUAL{$s}{$l}{CP}{perf}//0; my $CPm=$ANCHOR_QUAL{$s}{$l}{CP}{miss}//0;
            my $dfw=$DIR_COUNTS{$s}{$l}{forward}//0;  my $drc=$DIR_COUNTS{$s}{$l}{revcomp}//0;
            print $OUT join("\t",$s,$l,$HPp,$HPm,$OPp,$OPm,$CPp,$CPm,$dfw,$drc), "\n";
        }
    }
    close($OUT); log_msg("INFO","Wrote anchor_quality_counts: $f");
}
sub write_gap_stats {
    my ($f)=@_; open(my $OUT,">",$f) or die "Cannot write $f: $!\n";
    print $OUT join("\t",qw(sample lib_id metric n mean sd min max)), "\n";
    for my $s (sort keys %GAP_STATS){
        for my $l (sort keys %{ $GAP_STATS{$s} }){
            next unless exists $LIB_IDS{$l};
            for my $m (sort keys %{ $GAP_STATS{$s}{$l} }){
                my $h=$GAP_STATS{$s}{$l}{$m};
                my ($n,$sum,$sumsq,$min,$max)=@{$h}{qw(n sum sumsq min max)};
                my $mean = $n ? $sum/$n : 0;
                my $var  = ($n>1) ? ($sumsq/$n - $mean*$mean) : 0; $var=0 if $var<0;
                my $sd   = sqrt($var);
                printf $OUT "%s\t%s\t%s\t%d\t%.6f\t%.6f\t%d\t%d\n", $s,$l,$m,$n,$mean,$sd,$min,$max;
            }
        }
    }
    close($OUT); log_msg("INFO","Wrote gap_stats: $f");
}
sub write_cycle_distribution {
    my ($f)=@_; open(my $OUT,">",$f) or die "Cannot write $f: $!\n";
    print $OUT join("\t",qw(sample lib_id cycles count)), "\n";
    for my $s (sort keys %CYCLE_DIST){
        for my $l (sort keys %{ $CYCLE_DIST{$s} }){
            next unless exists $LIB_IDS{$l};
            for my $c (sort {$a<=>$b} keys %{ $CYCLE_DIST{$s}{$l} }){
                print $OUT join("\t",$s,$l,$c,$CYCLE_DIST{$s}{$l}{$c}), "\n";
            }
        }
    }
    close($OUT); log_msg("INFO","Wrote cycle_distribution: $f");
}
sub write_bb_frequency_by_cycle {
    my ($f)=@_; open(my $OUT,">",$f) or die "Cannot write $f: $!\n";
    print $OUT join("\t",qw(sample lib_id cycle bb_id count)), "\n";
    for my $s (sort keys %BB_FREQ){
        for my $l (sort keys %{ $BB_FREQ{$s} }){
            next unless exists $LIB_IDS{$l};
            for my $c (sort {$a<=>$b} keys %{ $BB_FREQ{$s}{$l} }){
                for my $bb (sort keys %{ $BB_FREQ{$s}{$l}{$c} }){
                    print $OUT join("\t",$s,$l,$c,$bb,$BB_FREQ{$s}{$l}{$c}{$bb}), "\n";
                }
            }
        }
    }
    close($OUT); log_msg("INFO","Wrote bb_frequency_by_cycle: $f");
}
sub write_length_window_stats {
    my ($f)=@_;
    open(my $OUT, ">", $f) or die "Cannot write $f: $!\n";
    # 멀티 모드 대응 컬럼 추가(fastp_peaks, active_centers, within_multi)
    print $OUT join("\t", qw(sample within_3 within_4 within_both within_multi outside policy tol exp3 exp4 active_mode fastp_peak fastp_peaks active_centers fastp_json)), "\n";
    for my $s (sort keys %LEN_WIN){
        my $w3  = $LEN_WIN{$s}{within_3}    // 0;
        my $w4  = $LEN_WIN{$s}{within_4}    // 0;
        my $wb  = $LEN_WIN{$s}{within_both} // 0;
        my $wm  = $LEN_WIN{$s}{within_multi}// 0;
        my $out = $LEN_WIN{$s}{outside}     // 0;
        my $conf = $LEN_CONF{$s} // { policy=>$LENGTH_FILTER_MODE, tol=>$LENGTH_TOL_BASE, exp3=>$EXP_LEN_3_BASE, exp4=>$EXP_LEN_4_BASE, active_mode=>'both' };
        my $peaks_str   = (exists $conf->{fastp_peaks} && ref($conf->{fastp_peaks}) eq 'ARRAY') ? join(",", @{$conf->{fastp_peaks}}) : "NA";
        my $centers_str = (exists $conf->{active_centers} && ref($conf->{active_centers}) eq 'ARRAY') ? join(",", @{$conf->{active_centers}}) : "NA";
        print $OUT join("\t", $s, $w3, $w4, $wb, $wm, $out,
            $conf->{policy}, $conf->{tol}, $conf->{exp3}, $conf->{exp4}, $conf->{active_mode},
            (defined $conf->{fastp_peak} ? $conf->{fastp_peak} : "NA"),
            $peaks_str, $centers_str,
            ($conf->{fastp_json} // "NA")
        ), "\n";
    }
    close($OUT);
    log_msg("INFO","Wrote length_window_stats: $f");
}
sub write_qc_checks {
    my ($f)=@_;
    open(my $OUT, ">", $f) or die "Cannot write $f: $!\n";
    print $OUT join("\t", qw(sample codon_len_not9 order_violation cp_len_out_of_27_29)), "\n";
    for my $s (sort keys %SAMPLES){
        my $c1 = $QC{$s}{codon_len_not9}      // 0;
        my $c2 = $QC{$s}{order_violation}     // 0;
        my $c3 = $QC{$s}{cp_len_out_of_27_29} // 0;
        print $OUT join("\t", $s, $c1, $c2, $c3), "\n";
    }
    close($OUT);
    log_msg("INFO","Wrote qc_checks: $f");
}

# --------- Main ---------
log_msg("INFO","Merged dir = $MERGED_DIR");
die "Fixed BB file not found: $FIXED_BB_FILE\n" unless -s $FIXED_BB_FILE;
die "Merged dir not found: $MERGED_DIR\n"    unless -d $MERGED_DIR;

load_fixed_bb_and_build_tries($FIXED_BB_FILE);
process_merged_dir($MERGED_DIR);

write_raw_matrix($RAW_MATRIX_FILE);
write_library_totals($LIB_TOTALS_FILE);
write_scaled_matrix($SCALED_MATRIX_FILE);
write_summary_basic($SUMMARY_FILE);

write_sample_stats($SAMPLE_STATS_FILE);
write_lib_sample_stats($LIB_SAMPLE_STATS_FILE);
write_anchor_qual_counts($ANCHOR_QUAL_FILE);
write_gap_stats($GAP_STATS_FILE);
write_cycle_distribution($CYCLE_DIST_FILE);
write_bb_frequency_by_cycle($BB_FREQ_BY_CYCLE_FILE);

write_length_window_stats($LEN_STATS_FILE);
write_qc_checks($QC_CHECKS_FILE);

if (open(my $DM, ">", $DONE_MARKER)) { print $DM scalar(localtime()), "\n"; close($DM); }
else { log_msg("WARN", "Could not write completion marker $DONE_MARKER: $!"); }
log_msg("INFO","All done.");
exit 0;

############################
# Helpers
############################
sub print_usage_and_exit {
    my ($code) = @_; $code //= 0;   # 0 for -h; 2 when invoked without any argument
    print <<"USAGE";
Usage:
  perl $0 -m <MERGED_DIR> -b <FIXED_BB_TSV> [options]
  perl $0 --merged-dir <DIR> --fixed-bb-file <FILE> [options]

Required inputs
  -m, --merged-dir <DIR>         Directory with merged fastq files
  -b, --fixed-bb-file <FILE>     DEL BB info (fixed) TSV path

I/O / Layout
  -o, --out-dir <DIR>            Output dir (default: \$Bin/<target>/02_decoded)
  -p, --target-protein <NAME>    Target folder under script dir (default: $OPT{target_protein})
  -n, --fastp-outdir-name <NAME> Subfolder name with merged fastq (default: $OPT{fastp_outdir_name})
  -R, --fastq-regex <REGEX>      Perl regex for merged fastq selection
                                 (default: AUTO: /(?:\\.fpmerged\\.fq(?:\\.gz)?|_merged\\.(?:fq|fastq)(?:\\.gz)?)\$/i)

Scaling (CPM)
  -S, --scaling-denom <POLICY>   sample_decoded | sample_length_passed | sample_total_raw | library_decoded
                                 (default: $OPT{scaling_denom})

Length filter (base)
  -M, --length-filter-mode <MODE>  reject_outside | reject_inside (default: $OPT{length_filter_mode})
  -T, --length-tol <INT>           Tolerance for expected length windows (default: $OPT{length_tol})
  -H, --hp-len <INT>               Canonical HP length (default: 8)
  -O, --op-len <INT>               Canonical OP length (default: 20)
  -c, --codon-len <INT>            Canonical Codon length (default: 9)
  -Q, --cp-len <INT>               Canonical CP length (default: 28)
  -J, --adj-tol <INT>              Adjacency tolerance between HP/OP/Codon/CP (default: $OPT{adj_tol})

fastp-based auto adjust
  -A, --auto-insert-detect 0|1     Enable auto detection from fastp JSON (default: $OPT{auto_insert_detect})
  -W, --auto-select-window 0|1     Select active window (3/4) from support (mixed => both; no support => keep multi)
                                   (default: $OPT{auto_select_window})
  -C, --auto-center-window 0|1     Center selected canonical window to nearest peak, requires -W 1; shift capped at
                                   max(2*length-tol, 6) nt (default: $OPT{auto_center_window})
  -U, --auto-tol-from-fastp 0|1    If width available, set TOL to width/2 (clamped 2..20) (default: 0)

Mixed DEL (Top-K peaks)
  -K, --peak-top-k <INT>           Number of peaks to use from fastp histogram (default: 5)
  -E, --peak-min-sep <INT>         Minimal separation between peaks (nt) (default: 2)
  -F, --peak-min-frac <FLOAT>      Minimal fraction per peak (0..1) (default: 0.01)

General
  -x, --max-failed-dump <INT>      Max undecoded dump lines per sample (default: $OPT{max_failed_dump})
  -P, --max-cp-cands <INT>         Max CP candidates to consider (0 = no limit; default: $OPT{max_cp_cands})
  -Y, --max-anchor-cands <INT>     Max HP/OP candidates to consider (0 = no limit; default: $OPT{max_anchor_cands})
      --mismatch <MODE>            1-bp mismatch indexing: hp_op_cp (default) | none
  -h, --help                       Show this help and exit

Notes
- Run without options to see this help (exit code 2); -h exits 0.
- Per-sample policy actually used is saved to length_window_stats.tsv, incl. fastp peaks and active centers.
- Input fastqs auto-detected: *.fpmerged.fq(.gz), *_merged.fq/.fastq(.gz).
USAGE
    exit $code;
}

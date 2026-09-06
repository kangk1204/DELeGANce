#!/usr/bin/env perl
# 01_preprocess_reads.pl (refactored to unified style)
# - 기능/논리/입출력 규격: 변경 없음
# - 형태 통일점:
#   * Embedded LibraryTrie (공통): insert/insert_one_mismatch/search_substrings/_by_lib_id
#   * ensure_dir: Windows/Unix 안전
#   * 로깅(ts/log_msg) + 경고/예외 트랩 공통
#   * 섹션 레이아웃/주석 규약 통일
# - 설계 고정: CODON mismatch는 비활성, HP/OP/CP mismatch는 옵션(--mismatch hp_op_cp)일 때만 치환 기반 1-bp 인덱싱

############################
# LibraryTrie.pm (embedded, unified)
############################
package LibraryTrie;
use strict;
use warnings;

sub new {
    my ($class, %opt) = @_;
    # die_on_dup: 같은 'sequence' 노드에서 동일 lib_id에 대해 2회 이상 insert 시 에러(디코딩 스크립트에서 사용)
    my $die_on_dup = $opt{die_on_dup} // 0;
    return bless { children => {}, is_end => 0, values => {}, die_on_dup => $die_on_dup }, $class;
}

sub insert {
    my ($self, $sequence, $info, $lib_id) = @_;
    die "Error: 'info' must be defined.\n"   unless defined $info;
    die "Error: 'lib_id' must be defined.\n" unless defined $lib_id;

    my $node = $self;
    foreach my $char (split //, $sequence) {
        $node->{children}{$char} //= { children => {}, is_end => 0, values => {} };
        $node = $node->{children}{$char};
    }
    $node->{is_end} = 1;

    # values->{$lib_id}는 배열로 통일(디코딩에서는 die_on_dup로 단건 보장)
    if ($self->{die_on_dup}) {
        die "Error: dup sequence for lib_id in trie: $sequence / $lib_id\n"
            if exists $node->{values}{$lib_id} && @{$node->{values}{$lib_id}};
    }
    $node->{values}{$lib_id} //= [];
    push @{ $node->{values}{$lib_id} }, $info;
}

# 치환 기반 1-bp mismatch(삽입/삭제는 여기서 다루지 않음; 디코딩 스크립트는 자체로 indel 생성)
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

# 읽기 서열 내 substring 매칭(모든 lib_id/모든 info 반환)
sub search_substrings {
    my ($self, $search_seq) = @_;
    my @matches;
    for (my $i = 0; $i < length($search_seq); $i++) {
        my $node = $self; my $j=$i;
        while ($j < length($search_seq)) {
            my $char = substr($search_seq, $j, 1);
            last unless exists $node->{children}{$char};
            $node = $node->{children}{$char};
            if ($node->{is_end}) {
                foreach my $lib_id (keys %{$node->{values}}) {
                    my $aref = $node->{values}{$lib_id} // [];
                    for my $info (@$aref) {
                        push @matches, {
                            sequence => substr($search_seq, $i, $j-$i+1),
                            lib_id   => $lib_id,
                            info     => $info,
                            position => $i
                        };
                    }
                }
            }
            $j++;
        }
    }
    return \@matches;
}

# 특정 lib_id에 한정한 substring 매칭
sub search_substrings_by_lib_id {
    my ($self, $search_seq, $target_lib_id) = @_;
    die "Error: 'target_lib_id' must be defined.\n" unless defined $target_lib_id;
    my @matches;
    for (my $i = 0; $i < length($search_seq); $i++) {
        my $node = $self; my $j=$i;
        while ($j < length($search_seq)) {
            my $char = substr($search_seq, $j, 1);
            last unless exists $node->{children}{$char};
            $node = $node->{children}{$char};
            if ($node->{is_end} && exists $node->{values}{$target_lib_id}) {
                my $aref = $node->{values}{$target_lib_id} // [];
                for my $info (@$aref) {
                    push @matches, {
                        sequence => substr($search_seq, $i, $j-$i+1),
                        lib_id   => $target_lib_id,
                        info     => $info,
                        position => $i
                    };
                }
            }
            $j++;
        }
    }
    return \@matches;
}

############################
# main
############################
package main;
use strict;
use warnings;
use File::Spec;
use FindBin qw($Bin);
use IO::Handle;
use Getopt::Long qw(GetOptions);
use Cwd qw(abs_path);
use File::Spec::Functions qw(file_name_is_absolute);
use IPC::Open3;
use Symbol qw(gensym);
use POSIX qw(:sys_wait_h);

#--------------------------- CLI
my ($OPT_BBINFO, $OPT_FASTQ_DIR, $OPT_OUTDIR, $OPT_THREADS, $OPT_HELP);
my ($OPT_MISMATCH, $OPT_SKIP_FASTP) = ('hp_op_cp', 0);  # default: enable 1-bp mismatch for HP/OP/CP
# Optional: custom FASTQ pairing regex with named groups (?P<base>...) and (?P<read>1|2)
my ($OPT_FASTQ_REGEX);
# mismatch: none | hp_op_cp  (CODON mismatch disabled by design)

sub print_usage {
    my $script = $0;
    # Try to align with repo naming; fall back to auto-detect
    my $default_bb     = File::Spec->catfile($Bin, "00_BB_information.txt");
    my $default_fastq  = File::Spec->catdir($Bin, "00_original_files");
    my $default_base   = File::Spec->catdir($Bin, "DELeGANce_out");
    my $default_thr    = $ENV{FASTP_THREADS} // 4;
    my $regex_example  = q{(?P<base>.+)[._-](?:R)?(?P<read>[12]).*\.(?:fastq|fq)(?:\.gz)?$};

    print <<"USAGE";
Usage:
  perl $script --bbinfo <FILE> [--fastq_dir <DIR>] [--outdir <RUN_DIR>] [--threads <N>]
               [--mismatch <MODE>] [--skip-fastp]
               [--fastq-regex <REGEX>]
  perl $script -b <FILE> [-f <DIR>] [-o <RUN_DIR>] [-t <N>] [--mismatch none|hp_op_cp] [--skip-fastp]
USAGE
    print <<"USAGE";
  perl $script --help

Required:
  -b, --bbinfo FILE         BB information 파일 경로

Optional:
  -f, --fastq_dir DIR       원본 FASTQ 디렉토리 (기본: $default_fastq)
  -o, --outdir RUN_DIR      전체 출력 루트(런 폴더).
                            * 절대경로면 해당 경로 사용
                            * 상대경로면 $default_base/<RUN_DIR>
                            * 미지정 시: $default_base
  -t, --threads N           fastp 스레드 수 (기본: $default_thr)
      --mismatch MODE       1-bp mismatch 인덱싱 모드: hp_op_cp(기본) | none
                             (주의: CODON mismatch는 설계상 비활성화)
      --skip-fastp          fastp 단계 생략(인덱싱/검증만 수행)
      --fastq-regex REGEX   사용자 정의 페어링 정규식(선택). 명명 그룹 (?P<base>...) 과
                            (?P<read>[12])를 포함해야 함. 예) $regex_example
  -h, --help                이 도움말 출력 후 종료

Examples:
  perl $script -b $default_bb -o example_run
  perl $script -b /data/BB_info.tsv -f /data/fastq -o /data/runs/example_20250815 -t 6 --mismatch hp_op_cp
USAGE
}

my $ARGC = @ARGV;
GetOptions(
    'bbinfo|b=s'   => \$OPT_BBINFO,
    'fastq_dir|f=s'=> \$OPT_FASTQ_DIR,
    'outdir|o=s'   => \$OPT_OUTDIR,
    'threads|t=i'  => \$OPT_THREADS,
    'mismatch=s'   => \$OPT_MISMATCH,
    'skip-fastp'   => \$OPT_SKIP_FASTP,
    'fastq-regex|R=s' => \$OPT_FASTQ_REGEX,
    'help|h'       => \$OPT_HELP,
) or do { print STDERR "Invalid options.\n\n"; print_usage(); exit 2; };

if ($OPT_HELP || $ARGC == 0) { print_usage(); exit 0; }
$OPT_MISMATCH = lc($OPT_MISMATCH // 'hp_op_cp');
die "Invalid --mismatch: $OPT_MISMATCH (use none|hp_op_cp)\n"
    unless $OPT_MISMATCH =~ /^(none|hp_op_cp)$/;

#--------------------------- Paths & logging
my $OUTPUT_BASE   = File::Spec->catdir($Bin, "DELeGANce_out");
my $RUN_DIR       = do {
    if (defined $OPT_OUTDIR && $OPT_OUTDIR ne "") {
        file_name_is_absolute($OPT_OUTDIR)
            ? $OPT_OUTDIR
            : File::Spec->catdir($OUTPUT_BASE, $OPT_OUTDIR);
    } else {
        $OUTPUT_BASE;
    }
};
my $FASTQ_DIR     = $OPT_FASTQ_DIR // File::Spec->catdir($Bin, "00_original_files");
my $FASTP_DIR     = File::Spec->catdir($RUN_DIR, "01_fastp_out");
my $DECODE_DIR    = File::Spec->catdir($RUN_DIR, "02_decoded");
my $NORM_DIR      = File::Spec->catdir($RUN_DIR, "03_normalized");

my $BB_INFO_FILE       = $OPT_BBINFO // File::Spec->catfile($Bin, "00_BB_information.txt");
my $BB_INFO_FIXED_FILE = File::Spec->catfile($RUN_DIR, "BB_information_fixed.tsv");
my $BB_ID_MAP_FILE     = File::Spec->catfile($RUN_DIR, "BB_id_resolution_map.tsv");
my $BB_ID_CHGONLY_FILE = File::Spec->catfile($RUN_DIR, "BB_id_changes_only.tsv");
my $BB_ERR_FILE        = File::Spec->catfile($RUN_DIR, "BB_duplicated.txt");
my $LOG_FILE           = File::Spec->catfile($RUN_DIR, "01_preprocess_reads.log");

sub ensure_dir {
    my ($dir) = @_;
    return if -d $dir;
    # File::Path avoids interpolating the path into a shell command (quotes/special chars were unsafe)
    require File::Path;
    File::Path::make_path($dir, { error => \my $err });
    if ($err && @$err) { die "Failed to create directory: $dir\n"; }
}
ensure_dir($RUN_DIR);
ensure_dir($FASTP_DIR);
ensure_dir($DECODE_DIR);
ensure_dir($NORM_DIR);

open(my $LOG_FH, ">>", $LOG_FILE) or die "Cannot open log file $LOG_FILE: $!\n";
$LOG_FH->autoflush(1);
STDOUT->autoflush(1);
# Stage completion marker (the log is append-mode, so an old "All done." line must not be used as a resume marker).
my $DONE_MARKER = File::Spec->catfile($FASTP_DIR, ".preprocess_done");
unlink($DONE_MARKER) if -e $DONE_MARKER;

# Resolve BB information file (explicit or auto-detect)
sub _auto_detect_bbinfo {
    my ($dir) = @_;
    opendir(my $DH, $dir) or return undef;
    my @candidates;
    while (my $f = readdir($DH)) {
        next if $f =~ /^\./;
        next unless $f =~ /^00_.*BB[_-]?information.*\.(?:txt|tsv)(?:\.gz)?$/i;
        push @candidates, File::Spec->catfile($dir, $f);
    }
    closedir($DH);
    return undef unless @candidates;
    return $candidates[0] if @candidates == 1;
    # Prefer a file containing 'DELeGANce' or 'BB_information' if multiple
    my @prefer = grep { /DELeGANce|BB[_-]?information/i } @candidates;
    return $prefer[0] if @prefer;
    return $candidates[0];
}

if (defined $BB_INFO_FILE && -e $BB_INFO_FILE) {
    my $abs = eval { abs_path($BB_INFO_FILE) };
    $BB_INFO_FILE = $abs if defined $abs;
} else {
    my $det = _auto_detect_bbinfo($Bin);
    if (defined $det && -e $det) {
        $BB_INFO_FILE = $det;
        my $abs = eval { abs_path($BB_INFO_FILE) };
        $BB_INFO_FILE = $abs if defined $abs;
        log_msg("INFO", "Auto-detected BB info: $BB_INFO_FILE");
    } else {
        print STDERR "BB information file not found: $BB_INFO_FILE\n";
        print STDERR "Tried auto-detect under: $Bin\n\n";
        print_usage();
        exit 2;
    }
}

#--------------------------- Log helpers
sub ts { scalar localtime(); }
sub log_msg {
    my ($level, $msg) = @_;
    my $line = ts() . " [$level] $msg\n";
    print $line;
    print $LOG_FH $line;
}
$SIG{__WARN__} = sub { my $m = join("", @_); chomp($m); log_msg("WARN", $m); };
$SIG{__DIE__}  = sub { my $m = join("", @_); chomp($m); log_msg("DIE",  $m); CORE::die("$m\n"); };

#--------------------------- Tries & info store
my $CODON_PATTERN = LibraryTrie->new();           # die_on_dup=0 (preprocess에서는 중복 허용/집계)
my $HP_PATTERN    = LibraryTrie->new();
my $OP_PATTERN    = LibraryTrie->new();
my $CP_PATTERN    = LibraryTrie->new();

my @INFO_STORE;  # index -> {type, tag_id, bb_id, cycle, lib_id, smiles, compseq}
my %CNT_EXACT = (HP=>0, OP=>0, CP=>0, CODON=>0);
my %CNT_MM    = (HP=>0, OP=>0, CP=>0, CODON=>0);
my %CNT_ROWS  = (HP=>0, OP=>0, CP=>0, CODON=>0);

#--------------------------- Utils
sub normalize_smiles {
    my ($s) = @_;
    return "" unless defined $s;
    $s =~ s/^\s+|\s+$//g;
    $s =~ s/\s+//g;
    return $s;
}
sub normalize_dna_seq {
    my ($s) = @_;
    return "" unless defined $s;
    $s =~ s/[\x{A0}\s]+//g;  # 공백류/NBSP 제거
    $s = uc($s);
    return $s;
}
sub mismatch_enabled_for_type {
    my ($type) = @_;
    return 0 if $type eq 'CODON';               # hard rule: CODON mismatch disabled
    return 0 if $OPT_MISMATCH eq 'none';
    return 1 if $OPT_MISMATCH eq 'hp_op_cp' && ($type eq 'HP' || $type eq 'OP' || $type eq 'CP');
    return 0;
}
sub index_seq {
    my ($type, $seq, $info_idx, $lib_id) = @_;
    if    ($type eq 'HP') {
        $HP_PATTERN->insert($seq, $info_idx, $lib_id); $CNT_EXACT{HP}++;
        if (mismatch_enabled_for_type('HP')) { $HP_PATTERN->insert_one_mismatch($seq, $info_idx, $lib_id); $CNT_MM{HP} += 3*length($seq); }
    }
    elsif ($type eq 'OP') {
        $OP_PATTERN->insert($seq, $info_idx, $lib_id); $CNT_EXACT{OP}++;
        if (mismatch_enabled_for_type('OP')) { $OP_PATTERN->insert_one_mismatch($seq, $info_idx, $lib_id); $CNT_MM{OP} += 3*length($seq); }
    }
    elsif ($type eq 'CP') {
        $CP_PATTERN->insert($seq, $info_idx, $lib_id); $CNT_EXACT{CP}++;
        if (mismatch_enabled_for_type('CP')) { $CP_PATTERN->insert_one_mismatch($seq, $info_idx, $lib_id); $CNT_MM{CP} += 3*length($seq); }
    }
    elsif ($type eq 'CODON') {
        $CODON_PATTERN->insert($seq, $info_idx, $lib_id); $CNT_EXACT{CODON}++;
    }
}

#--------------------------- Pipeline
log_msg("INFO", "Starting preprocess.");
log_msg("INFO", "  FASTQ_DIR=$FASTQ_DIR");
log_msg("INFO", "  OUTPUT_RUN_DIR=$RUN_DIR");
log_msg("INFO", "  FASTP_DIR=$FASTP_DIR");
log_msg("INFO", "  DECODE_DIR=$DECODE_DIR");
log_msg("INFO", "  NORMALIZED_DIR=$NORM_DIR");
log_msg("INFO", "  MISMATCH_MODE=$OPT_MISMATCH (CODON mismatch disabled)");
log_msg("INFO", "  SKIP_FASTP=$OPT_SKIP_FASTP");

resolve_bb_id_smiles_and_write_fixed($BB_INFO_FILE, $BB_INFO_FIXED_FILE, $BB_ID_MAP_FILE, $BB_ID_CHGONLY_FILE);
check_bb_info_and_build_index($BB_INFO_FIXED_FILE, $BB_ERR_FILE);

if (!$OPT_SKIP_FASTP) {
    run_fastp();
} else {
    # When skipping fastp, ensure merged FASTQ files already exist.
    opendir(my $DH, $FASTP_DIR) or die "Cannot open $FASTP_DIR: $!\n";
    # same pattern as 02_decode_reads.pl AUTO_FASTQ_REGEX (accept *.fpmerged.fq(.gz) as well)
    my @merged = grep { /(?:\.fpmerged\.fq(?:\.gz)?|_merged\.(?:fq|fastq)(?:\.gz)?)$/i } readdir($DH);
    closedir($DH);
    if (!@merged) {
        die "SKIP_FASTP=1 but no merged FASTQ found in $FASTP_DIR. ".
            "Provide pre-merged *_merged.fq.gz files or rerun without --skip-fastp.\n";
    }
}
if (open(my $DM, ">", $DONE_MARKER)) { print $DM scalar(localtime()), "\n"; close($DM); }
else { log_msg("WARN", "Could not write completion marker $DONE_MARKER: $!"); }
log_msg("INFO", "All done. See log: $LOG_FILE");

#--------------------------- Step 1: bb_id-SMILES fix
sub resolve_bb_id_smiles_and_write_fixed {
    my ($in_file, $out_fixed, $out_map, $out_changes_only) = @_;
    log_msg("INFO", "Resolving bb_id ⇄ SMILES and writing fixed file: $out_fixed");

    my %uniq_smiles_by_bb;     # {bb_id}{lib_id}{smiles_norm} = 1
    my %order_by_bb;           # {bb_id}{lib_id} = [smiles_norm first-seen order]
    my %libs_by_bb;            # {bb_id}{lib_id}=1 to detect cross-library reuse

    # 1st pass: bb_id → [SMILES 순서] 수집
    open(my $IN1, "<", $in_file) or die "Cannot open $in_file: $!\n";
    while (my $line = <$IN1>) {
        chomp($line);
        next if $line =~ /^\s*$/ || $line =~ /^\s*#/;

        my @d = split /\t/, $line, -1;
        next if @d < 7;

        my ($bb_type_raw, $seq_raw, $bb_id, $cycle, $tag_id, $lib_id, $smiles_raw) = @d;
        next if (defined $bb_type_raw && $bb_type_raw =~ /^(type|bb[_-]?type)$/i);

        my $smiles = normalize_smiles($smiles_raw);
        next if (!defined $bb_id || $bb_id eq "");
        next if (!defined $smiles || $smiles eq "");
        my $lib_key = defined($lib_id) ? $lib_id : "";

        $uniq_smiles_by_bb{$bb_id}{$lib_key} //= {};
        $order_by_bb{$bb_id}{$lib_key}       //= [];
        $libs_by_bb{$bb_id}{$lib_key} = 1;
        if (!exists $uniq_smiles_by_bb{$bb_id}{$lib_key}{$smiles}) {
            $uniq_smiles_by_bb{$bb_id}{$lib_key}{$smiles} = 1;
            push @{$order_by_bb{$bb_id}{$lib_key}}, $smiles;
        }
    }
    close($IN1);

    # 2nd pass: 고정본 작성 (여기서 DNA 서열 정규화 수행)
    open(my $IN2,  "<", $in_file)   or die "Cannot open $in_file: $!\n";
    open(my $OUTF, ">", $out_fixed) or die "Cannot write $out_fixed: $!\n";

    my $header_seen = 0;

    while (my $line = <$IN2>) {
        chomp($line);
        next if $line =~ /^\s*$/ || $line =~ /^\s*#/;
        my @d = split /\t/, $line, -1;

        # 헤더 처리
        if (!$header_seen && defined($d[0]) && $d[0] =~ /^(type|bb[_-]?type)$/i) {
            $d[2] = "bb_id_fixed";
            print $OUTF join("\t", @d), "\n";
            $header_seen = 1;
            next;
        }

        # 데이터 라인 파싱
        my ($bb_type_raw, $seq_raw, $bb_id, $cycle, $tag_id, $lib_id, $smiles_raw) = @d;

        # ❶ 타입/SMILES 정규화
        my $bb_type = uc($bb_type_raw // "");
        my $smiles  = normalize_smiles($smiles_raw // "");
        my $lib_key = defined($lib_id) ? $lib_id : "";

        # ❷ DNA 서열 정규화 (필수) — 2번째 열을 무조건 정규화하여 출력에 반영
        my $seq_norm = normalize_dna_seq($seq_raw // "");
        $d[1] = $seq_norm;

        # 유효성 체크
        next if (!defined $bb_id || $bb_id eq "");
        next if ($bb_type !~ /^(HP|OP|CP|CODON)$/);
        next if ($seq_norm eq "");

        # ❸ IUPAC-like 문자열만 추가 정규화
        for (my $i = 0; $i <= $#d; $i++) {
            next if $i == 2; # bb_id/fixed는 별도 처리
            my $v = $d[$i];
            next unless defined $v && $v ne "";
            my $strip = $v;
            $strip =~ s/[\x{A0}\s]+//g;
            if ($strip ne "" && $strip =~ /^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+$/) {
                $d[$i] = uc($strip);
            }
        }

        # bb_id_fixed 결정
        my $fixed_bb_id = $bb_id;
        if (exists $order_by_bb{$bb_id}{$lib_key}) {
            my @smiles_list = @{$order_by_bb{$bb_id}{$lib_key}};
            if (@smiles_list > 1) {
                my %idx; my $k = 0;
                for my $s (@smiles_list) { $k++; $idx{$s} = $k; }
                my $kuse = $idx{$smiles} // 1;
                $fixed_bb_id = $bb_id . "." . $kuse;
            }
        }
        # If bb_id is reused across multiple libraries, namespace by lib_id to avoid cross-lib collisions
        if (exists $libs_by_bb{$bb_id} && scalar(keys %{$libs_by_bb{$bb_id}}) > 1) {
            my $lib_sanitized = $lib_key;
            $lib_sanitized =~ s/[^A-Za-z0-9_.-]/_/g;
            $fixed_bb_id = $fixed_bb_id . "_LIB" . $lib_sanitized;
        }

        # 출력 (type 은 대문자로 정규화해 기록 — 02_decode_reads.pl 은 대소문자 무시 비교)
        $d[0] = $bb_type;
        $d[2] = $fixed_bb_id;
        print $OUTF join("\t", @d), "\n";
    }
    close($IN2);
    close($OUTF);

    # 매핑 파일 출력
    open(my $MAP_ALL, ">", $out_map) or die "Cannot write $out_map: $!\n";
    print $MAP_ALL "original_bb_id\tlib_id\tfixed_bb_id\tSMILES_index\n";
    for my $bb_id (sort keys %order_by_bb) {
        for my $lib (sort keys %{$order_by_bb{$bb_id}}) {
            my @slist = @{$order_by_bb{$bb_id}{$lib}};
            my $base = $bb_id;
            if (exists $libs_by_bb{$bb_id} && scalar(keys %{$libs_by_bb{$bb_id}}) > 1) {
                my $lib_sanitized = $lib; $lib_sanitized =~ s/[^A-Za-z0-9_.-]/_/g;
                $base = $base . "_LIB" . $lib_sanitized;
            }
            if (@slist <= 1) {
                print $MAP_ALL "$bb_id\t$lib\t$base\t1\n";
            } else {
                my $i=0; for my $s (@slist) { $i++; print $MAP_ALL "$bb_id\t$lib\t$base.$i\t$i\n"; }
            }
        }
    }
    close($MAP_ALL);

    open(my $MAP_CHG, ">", $out_changes_only) or die "Cannot write $out_changes_only: $!\n";
    print $MAP_CHG "original_bb_id\tlib_id\tfixed_bb_id\tSMILES_index\n";
    for my $bb_id (sort keys %order_by_bb) {
        for my $lib (sort keys %{$order_by_bb{$bb_id}}) {
            my @slist = @{$order_by_bb{$bb_id}{$lib}};
            next unless @slist > 1;
            my $base = $bb_id;
            if (exists $libs_by_bb{$bb_id} && scalar(keys %{$libs_by_bb{$bb_id}}) > 1) {
                my $lib_sanitized = $lib; $lib_sanitized =~ s/[^A-Za-z0-9_.-]/_/g;
                $base = $base . "_LIB" . $lib_sanitized;
            }
            my $i=0; for my $s (@slist) { $i++; print $MAP_CHG "$bb_id\t$lib\t$base.$i\t$i\n"; }
        }
    }
    close($MAP_CHG);

    log_msg("INFO", "Fixed BB information written: $out_fixed");
    log_msg("INFO", "Mapping files written: $out_map, $out_changes_only");
}

#--------------------------- Step 2: read fixed file, build index
sub check_bb_info_and_build_index {
    my ($file, $file_out) = @_;
    log_msg("INFO", "# Reading FIXED BB info: $file");

    my $line_no = 0;
    my $header_checked = 0;
    my @bad_rows;
    my %codon_seen_per_lib;   # enforce codon uniqueness within lib_id
    my %anchor_seen_per_lib;  # HP/OP/CP seq reuse within lib_id/type
    my %dup_anchor_count;

    # reset duplicate report
    if (defined $file_out && $file_out ne '') {
        unlink $file_out if -e $file_out;
        open(my $HDR, '>', $file_out) or die "Cannot write $file_out: $!\n";
        print $HDR "issue\tlib_id\tseq\tbb_id_first\tbb_id_conflict\n";
        close($HDR);
    }

    open(my $FO, "<", $file) or die "Cannot open $file: $!\n";
    while (my $line = <$FO>) {
        $line_no++;
        chomp($line);
        next if $line =~ /^\s*$/;
        my @d = split /\t/, $line, -1;

        if (!$header_checked) { $header_checked = 1; next if (defined($d[0]) && $d[0] =~ /^(type|bb[_-]?type)$/i); }
        next if @d < 7;

        my ($bb_type_raw, $seq_raw, $bb_id_fixed, $cycle, $tag_id, $lib_id, $smiles_raw, $comp_seq_raw) = @d;
        my $bb_type = uc($bb_type_raw // "");
        my $seq     = uc($seq_raw     // "");
        my $smiles  = normalize_smiles($smiles_raw // "");
        my $compseq = $comp_seq_raw // "";

        if ($bb_type !~ /^(HP|OP|CP|CODON)$/) { push @bad_rows, "Line $line_no: invalid type='$bb_type_raw'"; next; }
        if (!$seq)                             { push @bad_rows, "Line $line_no: empty sequence"; next; }
        if (!$bb_id_fixed)                     { push @bad_rows, "Line $line_no: empty fixed bb_id"; next; }
        if (!$lib_id)                          { push @bad_rows, "Line $line_no: empty lib_id"; next; }

        $CNT_ROWS{$bb_type}++;

        my $info_idx = scalar @INFO_STORE;
        push @INFO_STORE, {
            type    => $bb_type,
            tag_id  => $tag_id,
            bb_id   => $bb_id_fixed,
            cycle   => $cycle,
            lib_id  => $lib_id,
            smiles  => $smiles,
            compseq => $compseq,
        };

        if ($bb_type eq 'CODON') {
            my $key = join('|||', $lib_id, $seq);
            if (exists $codon_seen_per_lib{$key} && $codon_seen_per_lib{$key} ne $bb_id_fixed) {
                open(my $ERRF, ">>", $file_out) or die "Cannot write $file_out: $!\n";
                print $ERRF join("\t", 'DUP_CODON', $lib_id, $seq, $codon_seen_per_lib{$key}, $bb_id_fixed), "\n";
                close($ERRF);
            } else {
                $codon_seen_per_lib{$key} = $bb_id_fixed;
            }
        } else {
            my $key = join('|||', $lib_id, $bb_type, $seq);
            if (exists $anchor_seen_per_lib{$key} && $anchor_seen_per_lib{$key} ne $bb_id_fixed) {
                open(my $ERRF, ">>", $file_out) or die "Cannot write $file_out: $!\n";
                print $ERRF join("\t", "DUP_$bb_type", $lib_id, $seq, $anchor_seen_per_lib{$key}, $bb_id_fixed), "\n";
                close($ERRF);
                $dup_anchor_count{$bb_type}++;
            } else {
                $anchor_seen_per_lib{$key} = $bb_id_fixed;
            }
        }

        index_seq($bb_type, $seq, $info_idx, $lib_id);
    }
    close($FO);

    log_msg("INFO", "Indexing completed.");
    log_msg("INFO", sprintf("Rows: HP=%d, OP=%d, CP=%d, CODON=%d", $CNT_ROWS{HP}||0, $CNT_ROWS{OP}||0, $CNT_ROWS{CP}||0, $CNT_ROWS{CODON}||0));
    log_msg("INFO", sprintf("Exact insertions: HP=%d, OP=%d, CP=%d, CODON=%d", $CNT_EXACT{HP}||0, $CNT_EXACT{OP}||0, $CNT_EXACT{CP}||0, $CNT_EXACT{CODON}||0));
    log_msg("INFO", sprintf("1-bp mismatch insertions (approx): HP=%d, OP=%d, CP=%d, CODON=0 (disabled)", $CNT_MM{HP}||0, $CNT_MM{OP}||0, $CNT_MM{CP}||0));

    if (@bad_rows) {
        log_msg("WARN", "Found format issues in FIXED BB info:");
        log_msg("WARN", "  " . $_) for @bad_rows;
    }
    if (%dup_anchor_count) {
        my $msg = join(", ", map { "$_=$dup_anchor_count{$_}" } sort keys %dup_anchor_count);
        log_msg("WARN", "Found duplicate anchor sequences per lib/type: $msg (see $file_out)");
    }
}

#--------------------------- Step 3: fastp
sub _assign_pair {
    my ($paired, $base, $read, $path) = @_;
    my $slot = ($read == 1) ? 'f' : 'r';
    if (exists $paired->{$base}{$slot}) {
        die "FASTQ pairing conflict: read $read of sample '$base' matched twice ".
            "('$paired->{$base}{$slot}' and '$path'). Use --fastq-regex with (?<base>...)(?<read>[12]) to disambiguate.\n";
    }
    $paired->{$base}{$slot} = $path;
}
sub run_fastp {
    my $threads = $OPT_THREADS // $ENV{FASTP_THREADS} // 4;
    my $ver = `fastp --version 2>&1`;
    if ($? != 0) { die "fastp not found in PATH. Please install fastp and try again.\n"; }
    chomp($ver);
    log_msg("INFO", "# fastp detected: $ver (threads=$threads)");
    # Build paired-end list robustly
    opendir(my $OD, $FASTQ_DIR) or die "Cannot open directory $FASTQ_DIR: $!";
    my @file = readdir($OD); closedir($OD);

    my %paired;
    my @failed_samples;
    my $ext_re = qr/(?:fastq|fq)(?:\.gz)?/i;
    # Optional custom pairing regex with named groups (?P<base>...), (?P<read>[12])
    my $PAIR_RE;
    if (defined $OPT_FASTQ_REGEX && $OPT_FASTQ_REGEX ne '') {
        eval { $PAIR_RE = qr/$OPT_FASTQ_REGEX/ };
        if ($@) {
            log_msg("WARN", "Invalid --fastq-regex: $OPT_FASTQ_REGEX ($@). Falling back to auto-detect.");
            undef $PAIR_RE;
        } else {
            log_msg("INFO", "Using custom --fastq-regex for pairing");
        }
    }
    FILE: for my $fn (@file) {
        next unless $fn =~ /\.$ext_re$/i || $fn =~ /\.(?:fastq|fq)(?:\.gz)?$/i;
        my ($base, $read);
        if ($PAIR_RE && $fn =~ $PAIR_RE) {
            $base = exists $+{base} ? $+{base} : undef;
            $read = exists $+{read} ? $+{read} : undef;
            unless (defined $base && defined $read && $read =~ /^[12]$/) {
                next FILE;
            }
            _assign_pair(\%paired, $base, $read, File::Spec->catfile($FASTQ_DIR, $fn));
            next FILE;
        }
        # Common bcl2fastq style: *_R1_*.fastq.gz / *_R2_*.fastq.gz
        # Order matters: the simple *_1.fastq.gz / *_2.fastq.gz style is the most specific (digit directly before
        # the extension) and must be tested FIRST. Testing the bcl2fastq "_R1" pattern first mis-pairs sample names
        # that themselves contain _R1/_R2 (e.g. NEG_R1_2.fastq.gz was taken as read 1 of sample "NEG").
        if ($fn =~ /^(.*?)[._-]([12])\.(?:fastq|fq)(?:\.gz)?$/i) {
            ($base, $read) = ($1, $2);
        }
        # bcl2fastq style: <base>_R1_001.fastq.gz / <base>_R1.fastq.gz — R[12] must be followed by a separator
        # token (lane/chunk) or directly by the extension, never by another word character.
        # The base is GREEDY so that a base containing "_R1_" itself (e.g. X_R1_L001_R1_001.fastq.gz) is split at
        # the LAST _R[12] token: X_R1_L001 / read 1.
        elsif ($fn =~ /^(.*)[._-]R([12])(?:[._-][^.]*)?\.(?:fastq|fq)(?:\.gz)?$/i) {
            ($base, $read) = ($1, $2);
        } else {
            next FILE;
        }
        _assign_pair(\%paired, $base, $read, File::Spec->catfile($FASTQ_DIR, $fn));
    }

    my @samples = sort keys %paired;
    if (!@samples) { die "No paired-end FASTQ files detected in $FASTQ_DIR (expected <sample>_1/_2.fastq.gz or <sample>_R1/_R2*.fastq.gz).\n"; }

    log_msg("INFO", "# Processing FASTQ files sequentially..");
    my $idx = 0;
    for my $sample (@samples) {
        $idx++;
        my $r1 = $paired{$sample}{f} // "";
        my $r2 = $paired{$sample}{r} // "";
        unless (-s $r1 && -s $r2) {
            # A missing mate or an empty file is an input error: it must not silently drop the sample from the run.
            log_msg("ERROR", sprintf("[%d/%d] %s: pair missing or empty (R1=%s, R2=%s)", $idx, scalar(@samples), $sample, $r1 || '-', $r2 || '-'));
            push @failed_samples, $sample;
            next;
        }

        log_msg("INFO", sprintf("[%d/%d] %s", $idx, scalar(@samples), $sample));
        my $out_json = File::Spec->catfile($FASTP_DIR, "${sample}.fastp.json");
        my $out_html = File::Spec->catfile($FASTP_DIR, "${sample}.fastp.html");
        my $out_r1   = File::Spec->catfile($FASTP_DIR, "${sample}_clean_R1.fq.gz");
        my $out_r2   = File::Spec->catfile($FASTP_DIR, "${sample}_clean_R2.fq.gz");
        my $out_mg   = File::Spec->catfile($FASTP_DIR, "${sample}_merged.fq.gz");

        my @cmd = (
            'fastp',
            '-i', $r1,
            '-I', $r2,
            '-o', $out_r1,
            '-O', $out_r2,
            '--detect_adapter_for_pe',
            '--json', $out_json,
            '--html', $out_html,
            '--overrepresentation_analysis',
            '--thread', $threads,
            '-m',
            '--merged_out', $out_mg,
            '-P', 1,
        );
        log_msg("INFO", "    > fastp -i $r1 -I $r2 -o $out_r1 -O $out_r2 -m --merged_out $out_mg --thread $threads");

        # Run and capture output to log file robustly
        my $err = gensym; my $ch_out = gensym; my $ch_err = gensym;
        my $pid = eval { open3(undef, $ch_out, $ch_err, @cmd) };
        if ($@) {
            log_msg("ERROR", "failed to start fastp: $@");
            push @failed_samples, $sample;
            next;
        }
        while (1) {
            my $rin = my $win = '';
            my $fo = fileno($ch_out); my $fe = fileno($ch_err);
            vec($rin, $fo, 1) = 1 if defined $fo;
            vec($rin, $fe, 1) = 1 if defined $fe;
            last unless $rin;
            my $n = select($win=$rin, undef, undef, 0.1);
            if ($n) {
                if (defined $fo && vec($win, $fo, 1)) {
                    my $line = eval { scalar(<$ch_out>) };
                    if (defined $line) { print $LOG_FH $line; print $line; }
                }
                if (defined $fe && vec($win, $fe, 1)) {
                    my $line = eval { scalar(<$ch_err>) };
                    if (defined $line) { print $LOG_FH $line; print $line; }
                }
            }
            my $kid = waitpid($pid, WNOHANG);
            last if $kid > 0;
        }
        my $exit = $?;
        # drain whatever is still buffered in the pipes after the child exited
        while (defined(my $line = <$ch_out>)) { print $LOG_FH $line; print $line; }
        while (defined(my $line = <$ch_err>)) { print $LOG_FH $line; print $line; }
        close($ch_out); close($ch_err);
        if ($exit != 0) {
            log_msg("ERROR", "fastp failed for sample '$sample' (exit=$exit). Check $LOG_FILE");
            push @failed_samples, $sample;
        }
    }
    if (@failed_samples) {
        die "fastp failed for ".scalar(@failed_samples)." sample(s): ".join(", ", @failed_samples).". See $LOG_FILE\n";
    }
}
__END__

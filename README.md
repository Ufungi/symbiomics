<p align="center">
  <img src="https://img.shields.io/badge/platform-Linux-blue?logo=linux&logoColor=white" alt="Platform">
  <img src="https://img.shields.io/badge/nextflow-%E2%89%A524.10-23aa62?logo=nextflow&logoColor=white" alt="Nextflow">
  <img src="https://img.shields.io/badge/run%20with-singularity%20%7C%20docker%20%7C%20conda-blue?logo=singularity&logoColor=white" alt="Containers">
  <img src="https://img.shields.io/badge/language-Nextflow%20%7C%20Python%20%7C%20Bash-informational?logo=gnu-bash&logoColor=white" alt="Language">
  <img src="https://img.shields.io/badge/license-MIT-lightgrey" alt="License">
</p>

<p align="right"><a href="README.en.md">English</a> | <strong>한국어</strong></p>

# symbiomics

의사결정 기반 진핵생물 mRNA-seq 정량 및 게놈 주석 Nextflow 파이프라인. 게놈을 먼저
읽고, 어떻게 처리할지 결정한 뒤, 왜 그렇게 결정했는지 코어 시간을 쓰기 전에
알려줍니다.

```bash
./setup.sh                          # conda 자동 탐지 -> nextflow 환경 구축
symbiomics init -p matsutake        # project/matsutake/ 생성 (input/, output/, config.yaml)
# project/matsutake/input/samplesheet.tsv 편집, genome.fasta를 input/에 추가
symbiomics run --project matsutake  # 파이프라인 전체 실행
```

성격이 전혀 다른 두 타깃에 대해 구축·검증했습니다: *Tricholoma matsutake*
(161 Mb, 소프트마스킹된 균류 게놈)와 *Pinus densiflora*(21.74 Gb,
haplotype-resolved 침엽수 게놈, 2,006개 contig, 최장 contig 2.008 Gb). 동일한
명령이 둘 다에서 돌아가며, 파이프라인이 내부적으로 다른 도구를 선택합니다.

---

## ✨ 핵심 기능

*   **도구 목록이 아니라 정책 엔진** — 게놈을 측정한 뒤
    [`decide_strategy.py`](bin/decide_strategy.py)가 `strategy.yml`을 작성해
    모든 선택에 대해 무엇을 실행할지, 왜 그런지, 무엇을 기각했는지, 어떤
    플래그로 재정의하는지를 기록합니다. `-entry strategy`로 실행 전체 계획을
    해석하고 멈추므로, 몇 주짜리 실행을 시작 전에 검토할 수 있습니다.
*   **취약한 실패가 아니라 우아한 성능 저하** — 무거운 트랙은 시간 예산을
    갖습니다. 예산을 초과한 트랙은 제외되고, 나머지로 실행이 완료되며, 그
    누락은 조용히 사라지는 대신 리포트에 기록됩니다.
*   **정렬은 한 번만** — 하나의 HISAT2 BAM 채널이 주석 증거, StringTie
    조립, read 카운팅으로 갈라지므로, 63개 샘플 프로젝트도 같은 read를 두 번
    정렬하지 않습니다.
*   **정량 엔진 두 가지** — 실제 정렬이 필요할 때(대립유전자 특이적 발현,
    향후 구조 주석 증거)는 [HISAT2](https://github.com/DaehwanKimLab/hisat2),
    빠르고 편향 보정된 decoy-aware 전사체 정량에는
    [Salmon](https://github.com/COMBINE-lab/salmon). `--quant_engine`으로
    실행마다 선택합니다.
*   **구조 주석은 선택이지 전제가 아님** — `-entry functional`은 단백질
    FASTA를 바로 기능주석으로 보냅니다. 이미 좋은 유전자 모델이 있다면(*P.
    densiflora*가 그렇듯) BRAKER를 돌려서 원점으로 돌아갈 필요가 없습니다.
*   **주석 없는 strandedness 추론** — 게놈에서 직접 스플라이스 접합부 모티프
    (`GT..AG` / `CT..AC`)를 읽습니다. 파이프라인의 이 시점에는 Salmon이나
    RSeQC 방식 추론을 돌릴 전사체가 아직 없기 때문입니다.
*   **우연이 아니라 설계로 haplotype을 고려함** — haplotype-resolved diploid
    조립체는 모든 카운터가 전제하는 "read 하나 = 위치 하나" 가정을 깨뜨립니다.
    `docs/haplotypes.md`가 그 함정을 명시하고, 파이프라인은 이를 우회하는
    세 가지 방어적인 방법을 구현합니다. [WASP](https://github.com/bmvdgeijn/WASP)로
    필터링한 read에 [phASER](https://github.com/secastel/phaser)를 적용하는
    실제 대립유전자 특이적 발현(ASE) 워크플로우도 포함됩니다.

---

## 파이프라인 개요

```
samplesheet.tsv ──► INPUT_CHECK (SE/PE 자동 감지) ────────────────────┐
                                                                       │
genome.fasta ──► GENOME_PREP ──► DECIDE_STRATEGY ──► strategy.yml ────┤
                                                                       │
                        ┌──────────────────────────────────────────────┘
                        ▼
              fastp ─► strandedness 프로브 ─► HISAT2 ─► 정렬 ─► CSI/BAI 인덱스
                                                                       │
                 ┌─────────────────────────────┬─────────────────────┼──────────────┐
                 ▼                             ▼                     ▼              ▼
          StringTie 조립                featureCounts / HTSeq    (증거,        phASER ASE
                                       또는 Salmon 정량           v0.3+)        (선택,
                                                                                haplotype-
                                                                                resolved
                                                                                게놈)
                                              │
                                              ▼
                                        MultiQC 리포트

  ── 독립 경로, 게놈 불필요 ──
  단백질 FASTA ──► -entry functional ──► taxon=fungi:  funannotate2 + f2a (미배선)
                                          taxon=other:  DIAMOND/Swiss-Prot + eggNOG-mapper +
                                                         dbCAN v3 (HMMER+DIAMOND+eCAMI)  [Phase A]
                                                         + InterProScan + KofamScan       [Phase B]
                                       ──► product-name 사다리 (AHRD 없음) ──► functional TSV
```

---

## 사전 요구사항

| | 요구사항 |
|---|---|
| OS | Linux |
| 오케스트레이터 | [Nextflow](https://www.nextflow.io/) ≥ 24.10 (Java 17–21 필요) |
| 런타임 | [Singularity](https://sylabs.io/singularity/) 또는 Docker(권장), 또는 Conda(부분 지원 — `conf/conda.config` 참고) |
| 저장공간 | 게놈 크기에 따라 다름; work 디렉터리에 게놈 크기의 3배 이상 확보 |
| GPU | 선택 — 이후 마일스톤 트랙(Helixer, TMbed)에서만 사용 |

---

## 설치

**빠른 시작** (conda가 설치되어 있다면 한 명령으로 끝납니다):

```bash
git clone https://github.com/Ufungi/symbiomics.git
cd symbiomics
./setup.sh          # conda 탐지 -> nextflow 실행 환경 구축 -> 검증 -> PATH 등록
# 새 셸을 열거나 source ~/.bashrc를 한 뒤:
symbiomics init -p matsutake         # project/matsutake/ 생성 (input/, output/, config.yaml)
symbiomics run --project matsutake   # 전체 파이프라인 실행
```

전체 설치 매뉴얼(conda 설치, 런타임/프로파일 선택, DB 프로비저닝, 첫 실행,
문제 해결): [`docs/installation.md`](docs/installation.md).

### 1. Nextflow 실행 환경 구축

`scripts/symbiomics` 런처는 Nextflow(≥ 24.10)와 Java 17–21 JDK를 담은
`nextflow` conda 환경을 자동으로 해석·사용합니다. **그 환경을 만드는 것은
`setup.sh`의 역할**입니다(단 한 번만 하면 됩니다):

```bash
./setup.sh              # 최초: envs/runner.yml로 nextflow 환경 생성
./setup.sh --update     # envs/runner.yml이 바뀐 뒤 재동기화
./setup.sh --check      # 설치 없이 현재 상태만 점검
```

개별 도구 conda 환경(`envs/*.yml`)은 `-profile conda`에서 Nextflow가
실행 시 자동으로 만듭니다 — 수동 `conda env create` 불필요.

### 2. 환경 점검

`setup.sh`가 런처를 PATH에 등록하므로(새 셸 또는 `source ~/.bashrc` 후) `symbiomics`를
어디서든 한 단어로 쓸 수 있습니다. 직접 부르려면 `scripts/symbiomics ...`를 쓰세요.

```bash
symbiomics preflight
```

Java 버전, Nextflow, Singularity/Docker, GPU, 디스크 공간, 데이터베이스
프로비저닝을 확인하며 — 실패마다 정확한 해결 명령을 출력합니다.
`symbiomics`(즉 `scripts/symbiomics`)가 `nextflow`를 직접 부르는 대신 지원되는
진입점입니다: 프로젝트 전용 conda 환경(`nextflow`)에 번들된 JDK로 `JAVA_CMD`를
지정하는데, Nextflow는 Java 17–21가 필요하지만 많은 호스트의 기본값은 Java 11이기
때문입니다.

```bash
which nextflow    # `conda activate nextflow` 이후에만 resolve 됨 — 또는
                   # 그냥 symbiomics을 쓰면 자동으로 처리됩니다
```

### 3. (선택) 기능주석용 데이터베이스 프로비저닝

```bash
symbiomics run . -entry download_dbs --db_dir /data/db/eukannot
```

mRNA-seq arm에는 필요하지 않습니다. `-entry functional`의 Phase B
모듈(InterProScan 6, dbCAN v5, KofamScan)에만 필요하며 — Phase A(DIAMOND/
Swiss-Prot + eggNOG-mapper + dbCAN v3)는 이 서버에 이미 있는 데이터베이스로
바로 돌아갑니다. `docs/functional_annotation.md` 참고.

---

## 입력 파일

| 파일 | 필요한 경우 | 설명 |
|---|---|---|
| `samplesheet.tsv` | mRNA-seq arm | 샘플당 한 줄: `sample_id`, `fastq_1`/`fastq_2` 또는 `bam` (나머지 열은 전부 선택/자동). SE/PE는 자동 감지됩니다. `assets/samplesheet.example.tsv` 참고. |
| `genome.fasta` | mRNA-seq arm, 구조 주석 | 참조 게놈. 소프트마스킹 여부는 무관 — 필요하면 파이프라인이 마스킹하고, 소문자가 전혀 없는 `--premasked` 주장은 신뢰하지 않고 거부합니다. |
| `proteome.fasta` | `-entry functional` | 예측된 단백질 서열. 이 경로에는 게놈도 samplesheet도 필요 없습니다. |
| `genomes.tsv` | 다중 게놈 배치 (예약) | 현재 파이프라인에서 읽히지 않음. 단일 유전체 실행은 config의 `genome:` 키만 사용합니다. `assets/genomes.example.tsv`는 향후 스키마 참조용. |
| `reference_proteomes.tsv` | product-name 사다리 (예약) | 현재 파이프라인에서 읽히지 않음. 기능주석 Phase B 배선 시 사용 예정. `assets/reference_proteomes.example.tsv`는 스키마 참조용. |

---

## 사용법

**가장 간단한 방법** — 프로젝트 단위로 시작합니다. `init -p <name>`이
`project/<name>/` 아래 시작 파일(input/)과 결과 폴더(output/), 설정을 만들고,
각 키는 Nextflow 파라미터에 매핑되며, `profile` 키는 `-profile`이 되고,
CLI 인자가 config보다 우선합니다:

```bash
symbiomics init -p matsutake          # project/matsutake/{input,output,config.yaml}
# project/matsutake/input/samplesheet.tsv 편집, genome.fasta를 input/에 추가
symbiomics run --project matsutake    # 전체 실행 (정렬, 조립, 카운트)
symbiomics run --project matsutake --genome other.fasta   # 게놈만 바꿔 재실행
```

**전체 실행 — 플래그로 직접** (정렬, 조립, 카운트):

```bash
symbiomics run . -profile singularity,local64 \
    --input samplesheet.tsv --genome genome.fasta \
    --taxon plant --project plant
```

**정량만, HISAT2 대신 Salmon** (BAM 생성 없음, decoy-aware 편향 보정):

```bash
symbiomics run . -profile singularity,local64 \
    --input samplesheet.tsv --quant_engine salmon \
    --transcript_fasta HA.CDS.fa --project plant
```

**기능주석만 — 구조 주석 없음, 게놈 없음:**

```bash
symbiomics run . -entry functional -profile singularity,local64 \
    --proteome HA.PEP.fa --genome_id Pinde_HA --taxon plant \
    --project Pinde_HA
```

**며칠짜리 실행을 시작하기 전에 실행 계획 검토:**

```bash
symbiomics run . -entry strategy \
    --genome genome.fasta --taxon plant --clade gymnosperm
```

더 보기: `docs/usage.md`(파라미터, 출력 구조), `docs/decisions.md`(정책
엔진이 어떻게 판단하는지), `docs/haplotypes.md`(diploid/haplotype-resolved
게놈), `docs/functional_annotation.md`.

---

## 단계 참조

| 단계 | 하는 일 | 실행 조건 |
|---|---|---|
| `INPUT_CHECK` | samplesheet 파싱, SE/PE 추론, 검증 | 항상 |
| `GENOME_PREP` / `DECIDE_STRATEGY` | 게놈 통계, 크기 등급, `strategy.yml` | 게놈을 읽는 모든 entry에서 항상 |
| `RNASEQ_ALIGN` | fastp, strandedness 추론, HISAT2, 정렬, 인덱싱 | `--steps`에 `align` 포함 시 |
| `RNASEQ_ASSEMBLE` | StringTie 샘플별 조립 + 병합 | `assemble` |
| `QUANTIFY` | featureCounts / HTSeq / Salmon, 병합된 매트릭스 | `quantify` |
| `FUNCTIONAL` | DIAMOND-vs-Swiss-Prot / eggNOG-mapper / dbCAN v3(HMMER+DIAMOND+eCAMI) — Phase A, 배포됨; InterProScan / KofamScan / funannotate2 — Phase B, 아직 미배선 | `-entry functional` |
| `HAPLOTYPE_PAIRING` | minimap2 + SyRI: HA↔HB phased VCF와 allele-pairing 테이블 | `-entry pairing`, 독립 실행(그 출력이 아래 두 행에 공급됨) |
| `QUANTIFY` (Salmon, diploid 모드) | 합산 + haplotype별 allele 매트릭스 | `--quant_engine salmon --transcript_fasta HA.CDS.fa,HB.CDS.fa --haplotype_pairs <pairing 출력>` |
| `ALLELE_SPECIFIC_EXPRESSION` | WASP로 필터링한 HISAT2 read → phASER Gene AE | `--run_ase true --ase_phased_vcf <pairing 출력>` |
| `MULTIQC` | 통합 리포트, `versions.yml`, `strategy.yml` | 항상 |

---

## 사용 도구

| 도구 | 역할 | 인용 |
|---|---|---|
| [Nextflow](https://www.nextflow.io/) | 워크플로우 엔진 | Di Tommaso et al., *Nat. Biotechnol.* 2017 |
| [HISAT2](https://github.com/DaehwanKimLab/hisat2) | 스플라이스 인식 short-read 정렬 | Kim et al., *Nat. Biotechnol.* 2019 |
| [Salmon](https://github.com/COMBINE-lab/salmon) | decoy-aware, 편향 보정된 전사체 정량 | Patro et al., *Nat. Methods* 2017 |
| [StringTie](https://github.com/gpertea/stringtie) | 전사체 조립 | Pertea et al., *Nat. Biotechnol.* 2015 |
| [Subread/featureCounts](https://subread.sourceforge.net/) | exon 단위 read 카운팅 | Liao et al., *Bioinformatics* 2014 |
| [HTSeq](https://htseq.readthedocs.io/) | exon 단위 read 카운팅(교차검증) | Anders et al., *Bioinformatics* 2015 |
| [samtools](https://www.htslib.org/) | BAM 정렬/인덱싱/통계 | Danecek et al., *GigaScience* 2021 |
| [fastp](https://github.com/OpenGene/fastp) | read 트리밍/QC | Chen et al., *Bioinformatics* 2018 |
| [DIAMOND](https://github.com/bbuchfink/diamond) | Swiss-Prot best-hit product 전이(Phase A; UPIMAPI 자체의 ID-mapping 단계가 이 자리의 Phase B 업그레이드) | Buchfink et al., *Nat. Methods* 2021 |
| [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) | orthology 기반 GO/KEGG/description 전이 | Cantalapiedra et al., *Mol. Biol. Evol.* 2021 |
| [InterProScan](https://github.com/ebi-pf-team/interproscan6) | 도메인/family/GO 주석(Phase B) | Jones et al., *Bioinformatics* 2014 |
| [dbCAN / run_dbcan](https://github.com/bcb-unl/run_dbcan)(v3, HMMER+DIAMOND+eCAMI 합의) | CAZyme 주석 | Zheng et al., *Nucleic Acids Res.* 2023 |
| [funannotate2](https://github.com/nextgenusfs/funannotate2) | 균류 구조+기능 주석 | Palmer & Stajich |
| [minimap2](https://github.com/lh3/minimap2) | 조립체 간 정렬(HA↔HB) | Li, *Bioinformatics* 2018 |
| [SyRI](https://github.com/schneebergerlab/syri) | 두 조립체 간 synteny 및 구조변이 콜링 | Goel et al., *Genome Biol.* 2019 |
| [WASP](https://github.com/bmvdgeijn/WASP) | 대립유전자 특이적 read의 매핑 편향 필터링 | van de Geijn et al., *Nat. Methods* 2015 |
| [phASER](https://github.com/secastel/phaser) | read 기반 haplotype phasing 및 대립유전자 발현 | Castel et al., *Nat. Commun.* 2016 |
| [MultiQC](https://multiqc.info/) | 통합 QC 리포트 | Ewels et al., *Bioinformatics* 2016 |

DOI가 포함된 전체 인용 목록: `docs/citations.md`(각 도구의 모듈이 검증될
때마다 채워짐 — `scripts/check_containers.sh` 참고).

---

## 현재 상태

**v0.1.0** — mRNA-seq arm(samplesheet, strandedness 추론, HISAT2,
StringTie, 카운팅, 정책 엔진)이 구축되어 태그되었습니다.

**v0.2** — 두 번째 정량 엔진 Salmon, `-entry functional`(구조 주석과
분리된 단백질 입력 기능주석), haplotype-resolved 게놈 지원(HA↔HB pairing
`-entry pairing`, diploid Salmon 정량, phASER 대립유전자 특이적 발현
`--run_ase`)이 추가되었습니다. Salmon과 기능주석 Phase A(DIAMOND/Swiss-Prot
전이, eggNOG-mapper, dbCAN v3)는 이번 패스에서 실제 데이터로 검증되었고,
haplotype pairing과 ASE 워크플로우는 `-stub-run` 검증 및 synthetic
데이터 단위 테스트를 마쳤으며, 알려진 gap(WASP 컨테이너를 아직 빌드해야
함 — `docs/allele_specific_expression.md` 참고)이 route B의 실제 실행
전에 남아 있습니다. `docs/roadmap.md` 참고.

구조 주석(repeat, BRAKER, Helixer, consensus)은 의도적으로 뒤로 미뤘습니다
— 이 파이프라인의 참조 타깃에는 이미 좋은 외부 주석이 존재하므로, 기능주석과
정량을 먼저 우선했습니다.

---

## 라이선스

MIT. `LICENSE` 참고.

---

*진핵생물 게놈 주석 · mRNA-seq · Nextflow · haplotype-resolved 게놈 · 대립유전자 특이적 발현 · Pinus densiflora · Tricholoma matsutake*

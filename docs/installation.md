# 설치 (Installation)

symbiomics를 처음부터 실행 가능한 상태까지 갖추는 단계별 안내입니다. 두 런타임
경로(`-profile singularity` / `-profile docker` / `-profile conda`) 중 하나만
갖추면 되며, 어느 쪽이든 **conda 하나와 Nextflow를 실행하는 프로젝트 전용
`nextflow` 환경**(이름은 `SYMBIOMICS_NF_ENV`로 변경 가능)이 필요합니다.

대부분의 설치는 **한 명령**으로 끝납니다:

```bash
git clone https://github.com/Ufungi/symbiomics.git
cd symbiomics
./setup.sh          # conda 자동 탐지 -> nextflow 환경 구축 -> 검증 -> PATH 등록
# 새 셸을 열고(또는 source ~/.bashrc) 준비가 끝났으면:
symbiomics init -p matsutake         # project/matsutake/ 생성 (input/, output/, config.yaml)
symbiomics run --project matsutake   # 파이프라인 전체 실행
```

`setup.sh`가 끝나면 런처가 PATH에 등록되어, `scripts/symbiomics` 대신
`symbiomics`를 어디서든 한 단어로 부를 수 있습니다. 직접 부르고 싶다면
`scripts/symbiomics ...`를 그대로 써도 됩니다.

## 1. 사전 요구사항

| 항목 | 요구사항 | 비고 |
|---|---|---|
| OS | Linux | macOS는 미검증 |
| conda | Miniforge 또는 Miniconda | 파이프라인 자체에 conda는 없어도 되나, **Nextflow 실행 환경을 만드는 `setup.sh`가 conda를 필요로 함** |
| 런타임 | Singularity/Apptainer 또는 Docker(권장), 또는 Conda | `-profile`로 선택 (아래) |
| 디스크 | 게놈 크기의 3배 이상 (work 디렉터리) | preflight가 확인 |
| 저장공간 | 기능주석 DB는 별도 (선택) | `download_dbs`로 프로비저닝 |

## 2. conda 설치 (아직 없다면)

```bash
# Miniforge (conda-forge 기본 채널)
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
# 셸 재시작 후:
conda --version   # 확인
```

`setup.sh`는 conda를 자동으로 탐지합니다(`CONDA_EXE` → PATH의 `conda`).
탐지에 실패하면 `CONDA_EXE`를 conda 바이너리 경로로 지정해 주세요.

## 3. Nextflow 실행 환경 구축 (`setup.sh`)

`scripts/symbiomics` 런처는 Nextflow(≥ 24.10)와 Java 17–21 JDK를 담은
`nextflow` conda 환경을 요구합니다. (많은 호스트의 기본 Java는 11이라
순수 `nextflow` 바이너리가 기동을 거부합니다.) **이 환경을 만드는 것이
`setup.sh`의 역할이며, 저장소에는 그 외에 아무것도 미리 만들 필요가 없습니다.**

```bash
./setup.sh              # 최초: nextflow 환경 생성 (envs/runner.yml 기반)
./setup.sh --update     # envs/runner.yml이 바뀐 뒤 동기화
./setup.sh --check      # 설치 없이 현재 상태만 점검
./setup.sh --pull-images # (선택) singularity 이미지 미리 당겨두기
```

`setup.sh`가 하는 일:

1. conda 기반 경로를 탐지하고, 이 호스트의 `conda.sh` 경로를
   `conf/host.sh`에 기록합니다 — 런처(`scripts/symbiomics`)가 이 파일을 읽어
   더 이상 서버 전용 기본값(`/home/genome/...`)에 의존하지 않습니다.
2. `envs/runner.yml`에서 `nextflow` 환경을 생성/갱신합니다
   (`nextflow=24.10.4`, `openjdk=21`).
3. 환경 안의 Java 버전(17–21)과 nextflow 존재를 **검증**하고, 실패 시
   [FAIL] 줄을 출력하고 종료 코드 1로 실패합니다.

환경 이름은 기본값 `nextflow`이며 `SYMBIOMICS_NF_ENV`로
바꿀 수 있습니다. 예전 기본값인 `symbiomics` 환경이 남아 있으면
setup.sh가 실행 시 한 번 알려주며, 이 프로젝트는 그 환경을 사용하지
않습니다.

### 참고: `-profile conda`의 툴 환경은?

symbiomics의 개별 도구 conda 환경(`envs/*.yml` 22개)은 **Nextflow가 실행 시
자동으로 생성**합니다. 수동으로 `conda env create` 할 필요가 없습니다. 다만
일부 도구(Helixer, InterProScan 등)는 conda 레시피가 없어 이 프로파일에서
건너뛰므로, `conf/conda.config`의 "HONEST GAP LIST"를 확인하세요. 전체 툴
커버리지를 원하면 Singularity/Docker를 쓰는 것이 권장 경로입니다.

## 4. 런타임 선택 (`-profile`)

| profile | 의미 | 툴 커버리지 |
|---|---|---|
| `singularity` | 각 툴은 별도 컨테이너로 실행 (권장) | 전체 |
| `docker` | 동일하되 Docker 사용 | 전체 |
| `conda` | 각 툴은 conda 환경으로 실행 | 부분(일부 툴 건너뜀) |

실행은 항상 로컬 실행기 프로파일과 함께 결합하는 것이 일반적입니다:

```bash
scripts/symbiomics run . -profile singularity,local64 \
    --input samplesheet.tsv --genome genome.fasta --taxon plant --project plant
```

`conf/local64.config`는 실험실 서버(64코어/503GB) 기준이며, 다른 호스트는
자체 `.config` 또는 `-profile standard,conda` 등을 사용하세요.

## 4.5 프로젝트 시작: init -p (권장)

매 실행마다 플래그를 나열하는 대신, 프로젝트 단위로 시작합니다.
`init -p <name>`이 `project/<name>/` 아래 시작 파일(input/)과
결과 폴더(output/), 설정(config.yaml)을 한 번에 만듭니다:

```bash
symbiomics init -p matsutake   # project/matsutake/{input,output,config.yaml}
# project/matsutake/input/samplesheet.tsv 편집, genome.fasta를 input/에 추가
symbiomics run --project matsutake   # project/matsutake/config.yaml 자동 사용
```

`project/<name>/input/`에는 시작 파일이 복사되어 있습니다:
`samplesheet.tsv`, `genomes.tsv`, `reference_proteomes.tsv` — 필요한 것만
편집하면 됩니다.

수동 설정을 원하면 템플릿을 복사해 쓰는 방법도 그대로 유효합니다:

```bash
cp config.example.yaml config.yaml   # 주석 달린 템플릿을 복사해 편집
```

`config.yaml`의 각 키는 Nextflow 파라미터(`nextflow.config`)에 매핑되며,
`profile` 키는 `-profile` 플래그가 됩니다. 파일 경로는 **이 config 파일이
있는 디렉터리** 기준으로 해석되며, 파이프라인이 모르는 키는 무시됩니다.

```yaml
# config.yaml
project: plant               # output -> project/<project>/output/
input: input/samplesheet.tsv       # config 파일 기준 경로
genome: input/genome.fasta
taxon: plant
profile: "singularity,local64"
```

그리고 한 줄로 실행:

```bash
symbiomics run --project plant
```

명령줄에 직접 준 인자(플래그)가 config보다 우선합니다. 예를 들어 다음은
게놈만 바꿔 다시 실행합니다:

```bash
symbiomics run --project plant --genome other.fasta
```

`-entry`/`-profile` 같은 Nextflow 플래그도 그대로 함께 쓸 수 있습니다:

```bash
symbiomics run --project plant -entry strategy --clade gymnosperm
```

## 5. 환경 점검

```bash
symbiomics preflight
```

Java, Nextflow, 컨테이너 런타임, GPU, 디스크, `ulimit`, 로케일, 데이터베이스
디렉터리를 점검하며 **실패마다 정확한 해결 명령을 출력**합니다. 하드 실패
(Java/Nextflow)만 종료 코드 1을 냅니다.

## 6. (선택) 기능주석 데이터베이스 프로비저닝

`-entry functional`(Phase B: InterProScan 6, dbCAN v5, KofamScan)만 별도 DB가
필요합니다. mRNA-seq arm과 Phase A(DIAMOND/Swiss-Prot, eggNOG-mapper,
dbCAN v3)는 필요 없습니다.

```bash
symbiomics run . -entry download_dbs --db_dir /data/db/eukannot
```

`db_dir`의 기본값은 `SYMBIOMICS_DB_DIR`(기본 `/data/db/eukannot`)입니다.
상세: `docs/functional_annotation.md`.

## 7. 첫 실행

```bash
# 입력 파일 준비 (assets/samplesheet.example.tsv 참고)
cp assets/samplesheet.example.tsv samplesheet.tsv
# 전체 실행 -- config.yaml 경로
symbiomics run -c config.yaml
# 또는 플래그로 직접:
# symbiomics run . -profile singularity,local64 \
#     --input samplesheet.tsv --genome genome.fasta --taxon plant --project plant
```

며칠 걸리는 실행을 시작하기 전에 계획부터 검토하세요.

```bash
symbiomics run . -entry strategy \
    --genome genome.fasta --taxon plant --clade gymnosperm
```

## 8. 문제 해결

- `setup.sh`가 conda를 못 찾음 → conda 설치 후 PATH에 `conda`가 있는 셸에서 재실행,
  또는 `CONDA_EXE` 지정.
- `symbiomics: command not found` → `setup.sh` 후 새 셸을 열지 않았거나
  `source ~/.bashrc`를 안 했기 때문. 직접 `scripts/symbiomics ...`로 부르면 됩니다.
- `scripts/symbiomics`가 conda 프로파일을 못 찾음 → `./setup.sh`를
  먼저 실행해 `conf/host.sh`를 생성하거나, `SYMBIOMICS_CONDA_SH`를 직접 지정.
- `Cannot find Java or it's a wrong version` → `scripts/symbiomics`가
  `nextflow` 환경의 JDK를 사용하도록 되어 있으니, 환경이 올바른지
  `setup.sh --check`로 확인.

종합 문제 해결: `docs/troubleshooting.md`.

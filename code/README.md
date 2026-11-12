# GSE47460 COPD 잔여 DLCO 연구 — 분석 코드

기류제한(FEV1)과 CT 폐기종으로 설명되지 않는 COPD 확산능(DLCO) 저하의 폐조직 발현 특성을
두 Agilent 마이크로어레이 플랫폼(GPL14550 탐색 / GPL6480 재현)에서 분석한 코드다.
원고(Manuscript)의 모든 표·그림 수치는 이 코드의 출력에서 나온다.

## 폴더 구조

```
code/
├── README.md               이 파일
├── requirements.txt        Python 패키지
├── R_packages.txt          R 패키지
├── run_all.sh              전체 실행 (macOS/Linux, Git Bash)
├── run_all.bat             전체 실행 (Windows)
├── python/
│   ├── make_gmt.py                 유전자군 파일(travaglini_lung_c8.gmt) 생성
│   ├── 01_cohort_table1.py         Step 1  코호트(186/188), Table 1(QC 전), 임상변수 그림
│   ├── 03b_array_qc_metrics.py     Step 3b 배열 QC 지표(qc_array_*.csv) 재생성 ※
│   ├── 05_raw_qc_batch.py          Step 5  원시 Agilent FE 파일 헤더 → 스캔 날짜·STATS
│   ├── 06_table1_qc177.py          Step 6  주분석 177명 Table 1 (원고 표 1)
│   └── 07_manuscript_figures.py    Step 7  원고 그림 1–5, 부록 그림 S5
└── R/
    ├── 00_common.R                 공통 로더: 공변량 코딩, 유전자 집합, limma 모형, 배열 QC 규칙
    ├── 02_h1_limma.R               Step 2  H1 1차 분석(QC 전), VIF, 민감도
    ├── 02b_h2_preview_scores.R     (탐색적 미리보기, 확증 분석 아님)
    ├── 03_replication_diagnostics.R Step 3 재현 실패 진단(split-half, 검정력, 설계 일치, SVA)
    ├── 04_main_h1_to_h5.R          Step 4  주분석 H1–H5 (QC 177명)
    └── 05b_batch_tissue.R          Step 5b 조직 구성 점수, 스캔 배치 민감도
```

※ `03b_array_qc_metrics.py`: 원래 `qc_array_*.csv`를 만든 스크립트는 보존되지 않아 새로 작성했다.
저장된 파일과 두 플랫폼 모두 소수 셋째 자리까지 같은 값을 낸다(medcor, z_medcor, PC1, PC2, PC_maha 최대 차이 0).
`diag5_array_qc.csv`의 요약값도 생성 코드가 남아 있지 않으나, 같은 수치가 `04_main`의 `H1_summary.csv`
(noQC_186, primary_QC)와 `H3_replication.csv`에 재현된다.

## 입력 자료

자료 폴더 `GEO_GSE47460_COPD/`(NCBI GEO GSE47460에서 받은 원본 + `prepared/` 파생표)를 그대로 쓴다.
코드 폴더를 자료 폴더 안(`GEO_GSE47460_COPD/code/`)에 두면 `run_all`이 상대경로로 동작한다.
`python VERIFY_FILES.py`(자료 폴더 최상위)로 입력 파일 무결성(SHA256)을 먼저 확인한다.

## 실행 순서

| 단계 | 명령 (자료 폴더에서 실행) | 주요 출력 |
|---|---|---|
| 0 | `python code/python/make_gmt.py` (analysis/genesets에서) | travaglini_lung_c8.gmt |
| 1 | `python code/python/01_cohort_table1.py prepared/all_sample_metadata.tsv analysis/01_cohort` | cohort_main_186.tsv, cohort_sensitivity_188.tsv, cohort_flow.csv, table1_main_186.csv |
| 2 | `Rscript code/R/02_h1_limma.R prepared analysis/01_cohort analysis/02_h1` | h1_*.csv, fig2/fig3 (QC 전) |
| 3 | `Rscript code/R/03_replication_diagnostics.R prepared analysis/01_cohort analysis/03_diagnostics` | diag1–4, diag6, fig4 |
| 3b | `python code/python/03b_array_qc_metrics.py prepared analysis/01_cohort analysis/03_diagnostics` | qc_array_GPL14550.csv, qc_array_GPL6480.csv |
| 4 | `Rscript code/R/04_main_h1_to_h5.R prepared analysis/01_cohort analysis/genesets/travaglini_lung_c8.gmt analysis/04_main` | H1–H5 결과 |
| 5 | `python code/python/05_raw_qc_batch.py expansion/RAW analysis/01_cohort/cohort_main_186.tsv analysis/03_diagnostics analysis/05_raw_qc` | raw_fe_headers_cohort186.csv |
| 5b | `Rscript code/R/05b_batch_tissue.R prepared analysis/01_cohort analysis/genesets/travaglini_lung_c8.gmt analysis/05_raw_qc/raw_fe_headers_cohort186.csv analysis/05_raw_qc` | 조직 점수, 배치 민감도 |
| 6 | `python code/python/06_table1_qc177.py analysis/01_cohort analysis/03_diagnostics analysis/06_manuscript` | table1_qc177.csv |
| 7 | `python code/python/07_manuscript_figures.py analysis analysis/03_diagnostics analysis/06_manuscript/figures` | Figure1–5, FigureS5 |

Step 5는 `expansion/GSE47460_RAW.tar`를 `expansion/RAW/`에 풀어 둔 뒤 실행한다(1.76 GB).

## 핵심 분석 설정

- 공변량 코딩(`00_common.R`): sqrtEmph = √(%LAA-950), Male = 남성 1, Current = 현재 흡연 1(과거·비흡연 0)
- 유전자 집합: 심볼 있는 프로브 중 탐색 코호트 평균발현 상위 75%, 심볼당 평균발현 최고 프로브 1개 → 11,394개
- 연관 모형: 발현 ~ DLCO + FEV1 + sqrtEmph + Age + Male + Current, `limma::eBayes(trend = TRUE, robust = TRUE)`, 플랫폼별 별도 적합
- 배열 QC: 플랫폼별 다른 배열과의 상관 중앙값 < Q1 − 3×IQR → 제외 (GPL14550 7, GPL6480 2)
- H2 점수: GSVA(minSize = 10), 플랫폼 안에서 z 표준화. 유전자군은 검정 전에 고정
- 난수 시드: 20261020 (순열검정 1,000회, 부트스트랩 2,000회, split-half 200회)

## 실행 환경

분석은 R 4.x(limma, GSVA, sva)와 Python 3.10+ (pandas, numpy, scipy, matplotlib)에서 수행했다.
패키지 목록은 `R_packages.txt`, `requirements.txt` 참조.
그림의 한글 글꼴은 맑은 고딕 → Noto Sans CJK KR → AppleGothic 순으로 찾는다.

## 데이터 출처·인용

Tedrow J, Kaminski N, Guardela BJ, Schwartz DA. Gene expression profiling of chronic lung disease for the
Lung Genomics Research Consortium. NCBI GEO, GSE47460 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47460).
NCBI는 GEO 자료의 이용·배포를 제한하지 않는다(https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html).

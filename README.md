# FRPP 패키지 개발판

`emhawkes`의 지수 커널과 별도 잔차분포를 결합하는 R 패키지입니다.
별도 폴더와 Git 저장소에서 관리하며, 의존 방향은 `frpp → emhawkes`입니다.
CRAN 제출본은 그대로 유지하며, 공통 커널 연결은 `emhawkes` 개발판에서 검증합니다.

## 설치와 개발 API

```r
install.packages("remotes")
remotes::install_github("ksublee/frpp", build_vignettes = FALSE)
```

기본 `bridge` 우도와 `fsim()`·`ffit()`은 emhawkes 1.0.0 이상을 사용합니다.
`backend = "state"`는 `exp_kernel_replay()`를 제공하는 별도 emhawkes 개발판이
필요합니다. 이 API는 현재 로컬 개발 변경이며 emhawkes의 GitHub main이나
CRAN 배포본에는 아직 포함되지 않았습니다. API가 없으면 명시적으로 안내합니다.
현재 전체 테스트와 비네트는 해당 API가 설치된 개발 환경에서 검증했습니다.
공통 API 배포 전에는 기본 bridge를 사용하고 비네트 빌드는 생략해 주세요.

## 모형과 시뮬레이션

`fspec`은 `hspec`을 상속하지 않고 `kernel`과 `residual`을 보유합니다.
이전 시제품의 `frpp_spec()`과 저장 객체도 계속 사용할 수 있습니다.
`as_fspec(old_hspec)`으로 기존 FRPP 모형을 변환합니다.

```r
library(frpp)
model <- fspec(
  methods::new("hspec", mu = 0.3, alpha = 0.2, beta = 1),
  residual_law(c(shape = 2, rate = 2), stats::dgamma, stats::pgamma,
    stats::qgamma, stats::rgamma, name = "Gamma(2,2)"))
path <- fsim(model, size = 1000, seed = 73, lambda_component0 = 0)
# stats::simulate(model, size = 1000, seed = 73, lambda_component0 = 0)도 동일합니다.
logLik(model, path$inter_arrival, path$type, lambda_component0 = 0)
```

`size`는 초기 행을 포함합니다. 모든 글로벌 사건에서 모든 차원의 후보 잔차를 새로 생성합니다.
FRPP의 `auto` 초기화는 zero이며 Hawkes 정상 평균을 사용하지 않습니다.
구동량과 적분값을 일반 FRPP의 조건부 강도라고 바로 해석해서는 안 됩니다.

## 추정과 복원

```r
fit <- ffit(model, path$inter_arrival, path$type, lambda_component0 = 0)
summary(fit)
fdiagnostics(fit)
estimated <- fitted(fit)  # fitted_fspec(fit)도 같은 결과
do.call(logLik, c(list(object = estimated,
  inter_arrival = path$inter_arrival, type = path$type), finitial_state(fit)))
```

잔차분포 모수는 기본적으로 고정합니다. 함께 추정하려면 `estimate_residual = TRUE`를 지정합니다.
이때 커널과 잔차분포 척도가 식별되는지 확인하고 필요한 모수를 `fixed`로 고정해야 합니다.
임의 분포의 평균이나 정규화 규칙을 자동으로 추측하지 않습니다.
행렬 안의 같은 시작값은 기본적으로 같은 모수이며, `reduced = FALSE`로 풀 수 있습니다.
함수의 같은 모수 이름은 공유됩니다. 모수 공유·고정·매핑은 복원한 모형에도 유지합니다.

기본 `optimizer = "maxLik"` 외에 `optimizer = "stats"`를 사용할 수 있습니다.
`constraints`와 `constraint` 중 하나로 선형 제약을 전달합니다.
Stats 경로의 부등식은 내부 시작점이 필요하며 등식은 지원하지 않습니다.
`valid = function(param)`으로 분포별 유효 영역을 검사할 수 있습니다.
FALSE는 해당 후보를 배제하지만, 검사 함수와 사용자 callback의 오류는 숨기지 않습니다.

`coef()`, `logLik()`, `nobs()`, `vcov()`, `summary()`, `fitted()`와 저장·복원을 지원합니다.
곡률이 부적절하거나 경계 근처인 해에서는 자유 모수의 Wald 공분산을 NA로 반환합니다.
양의 공분산도 `available_unvalidated`이며 신뢰구간 포함률 검증을 뜻하지 않습니다.
수치 Hessian의 간격 민감성, 약한 식별성, tzexp의 비매끄러운 우도는 별도로 검토해야 합니다.

## 계산 경계

시뮬레이션과 기본 `bridge` 우도는 공개 emhawkes FRPP 계산에 위임합니다.
추정도 공개 `hfit()`의 모수 매핑·optimizer·복원을 재사용하지만 목적함수는 FRPP 우도입니다.
내부 함수 호출이나 Hawkes Hessian 안정화를 적용하지 않습니다.
`backend = "state"`는 공개 커널 replay의 적분·사건 직전 성분으로 FRPP 밀도·생존확률을 조합합니다.
이 경로에는 `exp_kernel_replay()`를 제공하는 emhawkes 개발판이 필요합니다.
이전 emhawkes에서는 기본 `bridge` 경로를 사용합니다.
숫자형 커널·행내 동일 감쇠율·숫자형 마크 효과·매핑을 지원하며 impact와 마크 밀도는 지원하지 않습니다.
사용자 callback 모형은 `bridge`에서 사용할 수 있습니다. 시뮬레이션은 행내 동일 감쇠율이 필요합니다.

상태 우도는 Hawkes 우도를 호출하지 않으며, `infer = FALSE`에서는 전체 커널 성분 이력을 저장하지 않습니다.
구간 적분과 성분 합의 메모리는 경로 길이와 사건 유형 수에 비례합니다. `infer = TRUE`는 기존 성분 이력을 반환합니다.
현재 replay는 R 기준 구현입니다. 10,000행·4유형·16커널의 3회 중앙값에서 이력 없는 우도는 누적 R 할당이 78.28→31.92MB로 줄었지만, 전체 이력은 79.33→88.01MB로 늘었습니다.
본체의 FRPP 분기를 제거하려면 별도 시뮬레이션 경로와 추정·복원 경계 정리가 더 필요합니다.
사용자 함수의 외부 환경은 참조로 유지되며 저장 시 외부 상태까지 자동 고정되지 않습니다.

## 예제와 검증

[비네트](vignettes/frpp-examples.Rmd)는 기존 FRPP 예제 3개를 포함합니다.
1차원 tzexp 3,000행, 2차원 tzexp 10,000행, 2차원 균등분포 5,000행을 비교합니다.
모든 비교에 seed 1107을 명시하며, 원래 seed가 없던 두 구간에도 비교용 seed를 지정했습니다.
시뮬레이션 출력과 최종 난수 상태, 호환 우도는 기존 표현과 정확히 일치했습니다.
독립 우도의 최대 절대 차이는 약 4.95e-10이었습니다.

```r
rmarkdown::render("vignettes/frpp-examples.Rmd")
# 10,000행 자료의 긴 추정까지 실행:
Sys.setenv(FRPP_RUN_LONG_EXAMPLES = "true")
rmarkdown::render("vignettes/frpp-examples.Rmd")
```

긴 추정은 기본 패키지 빌드에서 생략합니다. 회귀 검사는 기존 36개 시뮬레이션·난수 조합,
독립 적분, 해석해가 있는 추정, 고정·공유·제약·초기 상태·매핑·callback과 새 세션 복원을 포함합니다.
초기 시제품의 진단과 측정 기록은 [tools/prototype-results.md](tools/prototype-results.md)에 있습니다.

2026-10-17 추가 검증에서 `ffit()`의 10,000행 예제는 이전의 6개 추정값과 최대 로그우도를
정확히 재현했습니다. 기존의 부적절한 곡률은 새 진단에서 Wald 공분산 unavailable로 표시합니다.
테스트 5개 파일, 비네트 재실행, PDF 매뉴얼을 포함한 패키지 검사는 `Status: OK`였습니다.

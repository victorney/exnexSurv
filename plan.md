# exnexSurv — Nova função de utilidade: posterior & fit tools

**Objetivo:** adicionar, **sem tocar no C++ nem nos métodos existentes**, um conjunto de
funções R puras que transformam o ajuste (`exnex_surv`) em inferência útil de sobrevivência
e seleção de modelo, e uma função de simulação para testes/vignettes.

**Todas as funções operam apenas sobre o que já está retido no objeto:** `$draws`
(colunas `theta_1..theta_K`, `beta_1..beta_P`, `sigma2` — NOTA: na ausência de
covariáveis, `beta_*` **não existem**, logo P = 0) e `$data` (`time`, `event`, `group`,
`X`, `n`, `n_groups`, `n_covariates`, `cov_names`, `chain_seeds`).

## Escopo (após decisões com o usuário)

- NÃO modificar `src/exnex_mcmc.cpp` — Z não será retido (matriz `z_draws`).
  `prob_exchangeability()` fica para uma versão futura.
- Entregas:
  1. `survival_curves(fit, newdata, times, ...)`
  2. `plot.survival_exnex(...)` (S3 para o resultado de `survival_curves`)
  3. `median_survival(fit, newdata, ...)`
  4. `rmst(fit, tmax, newdata, ...)`
  5. `compute_waic(fit, ...)` e `compare_waic(...)`
  6. `probability_superiority(fit, a, b, newdata, function_of, times, tmax, ...)`
  7. `simulate_data(n_each, beta, sigma, ...)` (geração de dados de basket trial)

## Arquivos

- NOVO `R/survival.R` — `survival_curves`, `plot.survival_exnex`, `median_survival`,
  `rmst`
- NOVO `R/waic.R` — `compute_waic`, `compare_waic`
- NOVO `R/probability_superiority.R` — `probability_superiority`
- NOVO `R/simulate_data.R` — `simulate_data`
- NOVO `tests/testthat/test-survival.R`
- NOVO `tests/testthat/test-waic.R`
- NOVO `tests/testthat/test-probability_superiority.R`
- NOVO `tests/testthat/test-simulate_data.R`
- NOVO `vignettes/tools-for-experimenters.Rmd`
- EDITAR `DESCRIPTION` — adicionar `posterior` a Suggests
- EDITAR `_pkgdown.yml` — seção "Inference and model selection" na `reference`
- EDITAR `inst/WORDLIST` — palavras novas
- EDITAR `NEWS.md` — entrada `1.3.0`

## Restrições estruturais (ler ANTES de codar)

- `fit$draws` é `data.frame`; colunas θ, β, `sigma2`. Em fit sem covariável,
  `n_covariates == 0` e **não** há colunas `beta_*` → use `n_covariates` em vez de
  `colnames` hardcoded, e o vetor de covariates é `numeric(0)`.
- `fit$data$X` é a matriz de design (com `n` linhas, `n_covariates` colunas);
  `fit$data$cov_names` dá os nomes; `fit$data$time`/`event`/`group` são os vetores
  originais. Para projeções em uma grade de tempos, use as faixas de `time`.
- O vetor de efeito de cesta θ é **considerado centrado** de forma que o "efeito médio"
  das cestas está nos θ_j (o artigo usa `theta` em log-time). Equação:
  \[
  \log T_i = \theta_{g[i]} + X_i^\top\beta + \varepsilon,\quad
  \varepsilon\sim\mathcal N(0,\sigma^2).
  \]
- Em `compute_waic`, use o pacote `posterior::log_lin_pmf` e `posterior::rfun`? NÃO —
  use `loo` se desejado, mas implementação independente é mais simples: transforme cada
  desenho em `log_t = theta[g] + rowSums(X * beta)`, compute densidade da log-normal
  condicional em `log(time)`, para eventos use `dlnorm`, para censura use
  `log(1 - plnorm(...))`. Pontuação: `lpd = log(colMeans(exp(log_pd)))`,
  `p_waic = colVars(log_pd)` (amostral, com correção `n/(n-1)`), `waic = -2*(lpd - p_waic)`.
- `compare_waic(x, y, ...)` recebe **ajustes** (`exnex_surv`), computa `compute_waic` de
  cada um e retorna uma `data.frame` com colunas `model`, `waic`, `se_waic`,
  `lpd`, `p_waic`, `elpd_waic`. Ordenar por `waic` crescente. `se = sd(lpd)*sqrt(n)`.
- `plot.survival_exnex` deve retornar `invisible(x)`.
- Use `params(newdata)`? — `median_survival` precisa de X; use:
  `Xp <- as.matrix(newdata[, fit$data$cov_names])` se `n_covariates>0`, senão `numeric(0)`.
- Em todas as funções: checar `inherits(fit, "exnex_surv")`, checar `newdata`
  (data.frame opcional), e valores de `alpha`.

## Assinaturas e validações

- `survival_curves(fit, newdata = NULL, times = NULL, level = 0.95, ...)`.
  - Desenho: para cada desenho `d`, para cada linha de newdata `i`
    (ou 1 linha, se newdata NULL, com covariates `0`), calcule
    `logt = theta[g_i] + X_i dot beta`, `S(t) = 1 - pnorm((log(t) - logt) / sigma)`.
  - Exige covariáveis? Não — mas `theta` referencia a cesta. Se newdata tem grupo,
    usa como cesta; senão usa a primeira cesta e um aviso.
  - Retorna objeto `survival_exnex`: `data.frame` long com colunas
    `time`, `median` (mediana por tempo), `lower`, `upper` (quantis),
    `group` (rótulo da linha), e `level`. Múltiplas linhassão empilhadas.
- `plot.survival_exnex(x, ...)` desenha `S(t)` com faixa de credibilidade (ribbon) por
  grupo usando ggplot2. Retorna `invisible(x)`.
- `median_survival(fit, newdata = NULL, level = 0.95, ...)`.
  - Mediana de sobrevivência = menor `t` tal que `S(t)=0.5`: resolva
    `0.5 = 1 - pnorm((log t - logt)/sigma)` → `log t = logt` → `t = exp(logt)`.
    Ou seja, mediana por desenho `t_med = exp(logt_d)`. Para cada linha/sujeito agregue
    `median = exp(mean or quantile of logt))`. Retorna `data.frame` com
    `group`, `median`, `lower`, `upper` (quantis do `t_med` amostral).
  - Exige covariável? Não; usa média dos `logt` na linha (ou `0` sem novas covariáveis).
- `rmst(fit, tmax, newdata = NULL, level = 0.95, ...)`.
  - RMST até `tmax`: `RMST = \int_0^tmax S(t) dt`. Por desenho, **estime por quadratura
    média** que é barato: `S(t)` linha a linha, somatória de trapézio.
    Para cada desenho, `rmst_d = integrate S(t) dt` (aprox. trapezoidal sobre grade com,
    digamos, 500 pontos). Retorna `data.frame` com `group`, `rmst`, `lower`, `upper`.
  - Checar `tmax > 0` e `length(tmax) == 1`.
- `compute_waic(fit, ...)`.
  - Retorna `list(waic, se_waic, lpd, p_waic, elpd_waic)` e uma matriz `pointwise`
    (`n` linhas, `lpd`, `p_waic`, `waic` por observação).
  - Pontuação como descrito acima com `posterior` se disponível, senão com pure R.
    Como é Suggests, não exija `posterior` em runtime; implemente `log_dlnorm`/`log_1_pnorm`
    e `variance` à mão.
- `compare_waic(x, y, ...)`.
  - `...` pode incluir mais ajustes. Retorna `data.frame`.
- `probability_superiority(fit, a, b, newdata = NULL, function_of = "median",
  times = NULL, tmax = NULL, level = 0.95, ...)`.
  - Liga `survival_curves` (para `function_of %in% c("Surv_t","S_t","surv")` usa S(t) em
    `times`), `median_survival`, e `rmst`.
  - Retorna `list(prob, target, level)`.
- `simulate_data(n_each = 30, beta = c(0.5, -0.2), sigma = 1.2,
  resistant = NULL, resist_delta = -0.8, censoring_rate = 0.0,
  seed = NULL)`.
  - Gera K cestas com `n_each` sujeitos. Verdadeiros `theta_j` = execução lógica:
    baseline 0, margem por cesta; `resistant` indica índices de cesta cujas `theta`
    saem da "população" → `theta_j = theta_j + resist_delta`.
  - Retorna `data.frame` com `time`, `event`, `group` (fator), `X` (matriz design? não:
    retorne colunas `x1`, ... como parte do data.frame), `cov_names`, `n_each` etc.
  - Na verdade a API do fit espera um data.frame com `time`, `event`, `group` e covs.
  - Adicione `truth`? Não: guarde nos atributos: `attr(d, "true_theta")`,
    `attr(d, "true_beta")`, `attr(d, "true_sigma")`.
  - Valideção: `n_each` inteiro >0; `sigma>0`; `length(resistant)` <= K; `seed` escalar.

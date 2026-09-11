# Pointwise log-likelihood and WAIC from an pooling_surv fit

Computes the Watanabe-Akaike Information Criterion (WAIC) for a fitted
\`pooling_surv\` model. WAIC is a fully Bayesian information criterion
that uses the posterior draws of the log-likelihood and adds a penalty
for effective number of parameters \\p\_{waic}\\.

## Usage

``` r
compute_waic(fit, ...)
```

## Arguments

- fit:

  A fitted \`pooling_surv\` object.

- ...:

  Unused.

## Value

A named list with elements \`waic\`, \`se_elpd_waic\`, \`lpd\`,
\`p_waic\`, \`elpd_waic\`, and \`pointwise\` (an \`n\` by 3 matrix with
columns \`lpd\`, \`p_waic\`, \`waic\`).

## Details

For each observation \\i\\, with posterior draws indexed by
\\s=1,\ldots,S\\, we compute the pointwise log-likelihood \\\log p(y_i
\mid \theta^{(s)})\\, where the log-normal AFT model gives for an
observed event (\\\delta_i=1\\) \$\$\log p(y_i \mid \theta^{(s)}) =
\log\phi\\\left(\frac{\log t_i-\eta_i^{(s)}}{\sigma^{(s)}}\right) -
\log(\sigma^{(s)} t_i),\$\$ and for a censored observation
(\\\delta_i=0\\) \$\$\log p(y_i \mid \theta^{(s)}) =
\log\\\left\[1-\Phi\\\left(\frac{\log
t_i-\eta_i^{(s)}}{\sigma^{(s)}}\right)\right\],\$\$ where
\\\eta_i^{(s)}\\ is the linear predictor built from the draws.

WAIC is then \$\$lpd = \sum_i \log\\\left(\frac{1}{S}\sum_s \exp\log
p(y_i\mid\theta^{(s)})\right),\$\$ \$\$p\_{waic} = \sum_i
\mathrm{var}\_s\big(\log p(y_i\mid\theta^{(s)})\big),\$\$
\$\$\mathrm{WAIC} = -2(lpd - p\_{waic}).\$\$ Rows with extreme leverage
(\`p_waic\` above a large threshold) may flag issues; the function
returns a \`pointwise\` matrix so users can investigate.

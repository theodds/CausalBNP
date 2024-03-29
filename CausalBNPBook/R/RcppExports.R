# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

CCMVMI <- function(Y, R, omega, log_omega, gamma, log_gamma, log_1_m_gamma, beta, log_beta, log_1_m_beta) {
    .Call(`_CausalBNPBook_CCMVMI`, Y, R, omega, log_omega, gamma, log_gamma, log_1_m_gamma, beta, log_beta, log_1_m_beta)
}

CCMVGcomp <- function(N_sim, J, j_0, omega, log_omega, gamma, log_gamma, log_1_m_gamma, beta, log_beta, log_1_m_beta, xi) {
    .Call(`_CausalBNPBook_CCMVGcomp`, N_sim, J, j_0, omega, log_omega, gamma, log_gamma, log_1_m_gamma, beta, log_beta, log_1_m_beta, xi)
}

cluster <- function(y, X, matX, Sy, Sx, betaY, xPiPars, xMuPars, xSigPars, alphapsi, alphatheta, h0y, h0i, uniqueS, c0, mu0, nu0, tau, a0, b0, betainit, diagbetacov0, p1, ptx, p2) {
    .Call(`_CausalBNPBook_cluster`, y, X, matX, Sy, Sx, betaY, xPiPars, xMuPars, xSigPars, alphapsi, alphatheta, h0y, h0i, uniqueS, c0, mu0, nu0, tau, a0, b0, betainit, diagbetacov0, p1, ptx, p2)
}

cluster_continuous <- function(y, X, matX, Sy, Sx, betaY, sig2, xPiPars, xMuPars, xSigPars, alphapsi, alphatheta, h0y, h0i, uniqueS, c0, mu0, nu0, tau, a0, b0, betaa0, betab0, betainit, diagbetacov0, p1, ptx, p2) {
    .Call(`_CausalBNPBook_cluster_continuous`, y, X, matX, Sy, Sx, betaY, sig2, xPiPars, xMuPars, xSigPars, alphapsi, alphatheta, h0y, h0i, uniqueS, c0, mu0, nu0, tau, a0, b0, betaa0, betab0, betainit, diagbetacov0, p1, ptx, p2)
}

MARMI <- function(Y, R, omega, log_omega, beta, log_beta, log_1_m_beta) {
    .Call(`_CausalBNPBook_MARMI`, Y, R, omega, log_omega, beta, log_beta, log_1_m_beta)
}

LogSumExp <- function(x) {
    .Call(`_CausalBNPBook_LogSumExp`, x)
}

NIPMI <- function(Y, R, omega, log_omega, gamma, log_gamma, log_1_m_gamma, beta, log_beta, log_1_m_beta, j_0) {
    .Call(`_CausalBNPBook_NIPMI`, Y, R, omega, log_omega, gamma, log_gamma, log_1_m_gamma, beta, log_beta, log_1_m_beta, j_0)
}

ParafacGcomp <- function(N_sim, J, j_0, omega, log_omega, gamma, log_gamma, log_1_m_gamma, beta, log_beta, log_1_m_beta, xi) {
    .Call(`_CausalBNPBook_ParafacGcomp`, N_sim, J, j_0, omega, log_omega, gamma, log_gamma, log_1_m_gamma, beta, log_beta, log_1_m_beta, xi)
}

UpdateClassMARCpp <- function(Y, R, log_beta, log_1_m_beta, log_omega) {
    .Call(`_CausalBNPBook_UpdateClassMARCpp`, Y, R, log_beta, log_1_m_beta, log_omega)
}

TLOMI <- function(Y, R, omega, log_omega, gamma, log_gamma, log_1_m_gamma, beta, log_beta, log_1_m_beta, xi, j_0) {
    .Call(`_CausalBNPBook_TLOMI`, Y, R, omega, log_omega, gamma, log_gamma, log_1_m_gamma, beta, log_beta, log_1_m_beta, xi, j_0)
}

UpdateBetaCpp <- function(success_counts, failure_counts, col_shape_1, col_shape_2) {
    .Call(`_CausalBNPBook_UpdateBetaCpp`, success_counts, failure_counts, col_shape_1, col_shape_2)
}

UpdateSufficient <- function(Y, R, C, K) {
    .Call(`_CausalBNPBook_UpdateSufficient`, Y, R, C, K)
}

UpdateClassCpp <- function(Y, R, log_beta, log_1_m_beta, log_gamma, log_1_m_gamma, log_omega) {
    .Call(`_CausalBNPBook_UpdateClassCpp`, Y, R, log_beta, log_1_m_beta, log_gamma, log_1_m_gamma, log_omega)
}


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define N_SYMBOLS 4  // A,C,G,T
#define MAX_SEQ_LEN 5000000
#define MAX_ITER 10

typedef struct {
    int n_states;
    double *pi;       // 初始機率 (size n_states)
    double **A;       // 轉移矩陣 (n_states x n_states)
    double **B;       // 發射矩陣 (n_states x n_symbols)
} HMM;

//---------------- 初始化 HMM ----------------
HMM *init_hmm(int n_states) {
    HMM *hmm = malloc(sizeof(HMM));
    hmm->n_states = n_states;
    hmm->pi = malloc(n_states * sizeof(double));
    hmm->A = malloc(n_states * sizeof(double*));
    hmm->B = malloc(n_states * sizeof(double*));

    for (int i=0;i<n_states;i++) {
        hmm->A[i] = malloc(n_states * sizeof(double));
        hmm->B[i] = malloc(N_SYMBOLS * sizeof(double));
        hmm->pi[i] = 1.0/n_states;
        for (int j=0;j<n_states;j++) hmm->A[i][j] = 1.0/n_states;
        for (int j=0;j<N_SYMBOLS;j++) hmm->B[i][j] = 1.0/N_SYMBOLS;
    }
    return hmm;
}

//---------------- 字元映射 ----------------
int symbol_idx(char c) {
    switch(c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

//---------------- Forward演算法(log-space) ----------------
void forward(HMM *hmm, char *seq, int len, double **alpha, double *scale) {
    int N = hmm->n_states;
    scale[0] = 0.0;

    // 初始化
    for (int i=0;i<N;i++) {
        int s = symbol_idx(seq[0]);
        alpha[i][0] = log(hmm->pi[i]) + log(hmm->B[i][s]);
    }

    // 遍歷序列
    for (int t=1;t<len;t++) {
        int s = symbol_idx(seq[t]);
        for (int j=0;j<N;j++) {
            double sum = -INFINITY;
            for (int i=0;i<N;i++) {
                double val = alpha[i][t-1] + log(hmm->A[i][j]);
                if (sum==-INFINITY) sum=val;
                else sum = log(exp(sum)+exp(val)); // log-sum-exp
            }
            alpha[j][t] = sum + log(hmm->B[j][s]);
        }
    }
}

//---------------- Backward演算法(log-space) ----------------
void backward(HMM *hmm, char *seq, int len, double **beta) {
    int N = hmm->n_states;
    for (int i=0;i<N;i++) beta[i][len-1]=0.0;

    for (int t=len-2;t>=0;t--) {
        int s = symbol_idx(seq[t+1]);
        for (int i=0;i<N;i++) {
            double sum = -INFINITY;
            for (int j=0;j<N;j++) {
                double val = log(hmm->A[i][j]) + log(hmm->B[j][s]) + beta[j][t+1];
                if (sum==-INFINITY) sum=val;
                else sum = log(exp(sum)+exp(val));
            }
            beta[i][t] = sum;
        }
    }
}

//---------------- Baum-Welch(EM) ----------------
void baum_welch(HMM *hmm, char *seq, int len, int max_iter) {
    int N = hmm->n_states;
    int T = len;

    double **alpha = malloc(N*sizeof(double*));
    double **beta  = malloc(N*sizeof(double*));
    double *scale = malloc(T*sizeof(double));
    for (int i=0;i<N;i++) { alpha[i]=malloc(T*sizeof(double)); beta[i]=malloc(T*sizeof(double)); }

    for (int iter=0; iter<max_iter; iter++) {
        forward(hmm, seq, T, alpha, scale);
        backward(hmm, seq, T, beta);

        // 計算 γ[i][t] = P(state i at t)
        double **gamma = malloc(N*sizeof(double*));
        for (int i=0;i<N;i++) gamma[i]=malloc(T*sizeof(double));

        for (int t=0;t<T;t++) {
            double sum = -INFINITY;
            for (int i=0;i<N;i++)
                if (sum==-INFINITY) sum=alpha[i][t]+beta[i][t];
                else sum=log(exp(sum)+exp(alpha[i][t]+beta[i][t]));
            for (int i=0;i<N;i++) gamma[i][t] = exp(alpha[i][t]+beta[i][t]-sum);
        }

        // 更新 π
        for (int i=0;i<N;i++) hmm->pi[i]=gamma[i][0];

        // 更新 A
        for (int i=0;i<N;i++) {
            double denom = 0.0;
            for (int t=0;t<T-1;t++) denom += gamma[i][t];
            for (int j=0;j<N;j++) {
                double numer = 0.0;
                for (int t=0;t<T-1;t++) {
                    int s = symbol_idx(seq[t+1]);
                    numer += gamma[i][t]*hmm->A[i][j]*hmm->B[j][s]; // approximation
                }
                hmm->A[i][j] = numer/denom;
            }
        }

        // 更新 B
        for (int i=0;i<N;i++) {
            double denom = 0.0;
            for (int t=0;t<T;t++) denom += gamma[i][t];
            for (int k=0;k<N_SYMBOLS;k++) {
                double numer = 0.0;
                for (int t=0;t<T;t++)
                    if (symbol_idx(seq[t])==k) numer+=gamma[i][t];
                hmm->B[i][k] = numer/denom;
            }
        }

        for (int i=0;i<N;i++) free(gamma[i]);
        free(gamma);
    }

    for (int i=0;i<N;i++) { free(alpha[i]); free(beta[i]); }
    free(alpha); free(beta); free(scale);
}

//---------------- 計算序列的 log probability ----------------
double log_prob_seq(HMM *hmm, char *seq, int len) {
    int N=hmm->n_states;
    double **alpha = malloc(N*sizeof(double*));
    double *scale = malloc(len*sizeof(double));
    for (int i=0;i<N;i++) alpha[i]=malloc(len*sizeof(double));
    forward(hmm, seq, len, alpha, scale);
    double logp=-INFINITY;
    for (int i=0;i<N;i++)
        if (logp==-INFINITY) logp=alpha[i][len-1];
        else logp=log(exp(logp)+exp(alpha[i][len-1]));
    for (int i=0;i<N;i++) free(alpha[i]);
    free(alpha); free(scale);
    return logp;
}

//---------------- 主程式 ----------------
int main() {
    char *seq_train = "ACGTACGTACGTACGTACGT"; // 替換為實際抓取序列
    char *seq_test1 = "TGCATGCATGCATGCA"; 
    char *seq_test2 = "GATTACA"; 
    int len_train = strlen(seq_train);
    int len_test1 = strlen(seq_test1);
    int len_test2 = strlen(seq_test2);

    // HMM1: 2 states
    HMM *hmm1 = init_hmm(2);
    baum_welch(hmm1, seq_train, len_train, MAX_ITER);

    // HMM2: 3 states
    HMM *hmm2 = init_hmm(3);
    baum_welch(hmm2, seq_train, len_train, MAX_ITER);

    printf("HMM1 (2 states) log probabilities:\n");
    printf("  training: %.3f\n", log_prob_seq(hmm1, seq_train, len_train));
    printf("  test1:    %.3f\n", log_prob_seq(hmm1, seq_test1, len_test1));
    printf("  test2:    %.3f\n", log_prob_seq(hmm1, seq_test2, len_test2));

    printf("HMM2 (3 states) log probabilities:\n");
    printf("  training: %.3f\n", log_prob_seq(hmm2, seq_train, len_train));
    printf("  test1:    %.3f\n", log_prob_seq(hmm2, seq_test1, len_test1));
    printf("  test2:    %.3f\n", log_prob_seq(hmm2, seq_test2, len_test2));

    return 0;
}

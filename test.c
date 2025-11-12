#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define MAX_LEN 5000000
#define MAX_ITER 5
#define N_SYMBOLS 4
#define LOG_ZERO -1e100

typedef struct {
    int n_states;
    double *pi;
    double **A;
    double **B;
} HMM;

// ---------------- 基本函式 ----------------
char to_base(char c) {
    c = toupper(c);
    if (c=='A'||c=='C'||c=='G'||c=='T') return c;
    return 0;
}

int symbol_idx(char c) {
    switch(c){
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

// ---------------- 抓序列 ----------------
long fetch_sequence_safe(const char *filename, char *chrom, long start_pos, long length, char *buffer) {
    FILE *fp = fopen(filename,"r");
    if (!fp) { perror("Cannot open file"); return -1; }

    char line[1024];
    int in_chrom = 0;
    long seq_pos = 0;
    long idx = 0;

    while (fgets(line,sizeof(line),fp)) {
        if (line[0]=='>') {
            in_chrom = (strncmp(line+1, chrom, strlen(chrom))==0);
            continue;
        }
        if (!in_chrom) continue;

        for (int i=0; line[i] && line[i]!='\n'; i++) {
            char base = to_base(line[i]);
            seq_pos++;
            if (base != 0 && seq_pos >= start_pos && idx < length)
                buffer[idx++] = base;
            if (idx >= length) break;
        }
        if (idx >= length) break;
    }
    fclose(fp);
    return idx;
}

// ---------------- HMM 初始化 ----------------
HMM *init_hmm(int n_states) {
    HMM *hmm = malloc(sizeof(HMM));
    hmm->n_states = n_states;
    hmm->pi = malloc(n_states * sizeof(double));
    hmm->A = malloc(n_states * sizeof(double*));
    hmm->B = malloc(n_states * sizeof(double*));

    for (int i=0;i<n_states;i++) {
        hmm->A[i] = malloc(n_states*sizeof(double));
        hmm->B[i] = malloc(N_SYMBOLS*sizeof(double));
        hmm->pi[i] = 1.0/n_states;
        for (int j=0;j<n_states;j++) hmm->A[i][j] = 1.0/n_states;
        for (int j=0;j<N_SYMBOLS;j++) hmm->B[i][j] = 1.0/N_SYMBOLS;
    }
    return hmm;
}

// ---------------- log-sum-exp ----------------
double log_sum_exp(double a, double b) {
    if (a < LOG_ZERO) return b;
    if (b < LOG_ZERO) return a;
    if (a > b) return a + log1p(exp(b - a));
    else return b + log1p(exp(a - b));
}

// ---------------- Forward ----------------
void forward(HMM *hmm, char *seq, int len, double **alpha) {
    int N=hmm->n_states;
    for (int i=0;i<N;i++) {
        int s = symbol_idx(seq[0]);
        alpha[i][0] = log(hmm->pi[i]) + log(hmm->B[i][s]);
    }

    for (int t=1;t<len;t++) {
        int s = symbol_idx(seq[t]);
        for (int j=0;j<N;j++) {
            double sum = LOG_ZERO;
            for (int i=0;i<N;i++) {
                double val = alpha[i][t-1] + log(hmm->A[i][j]);
                sum = log_sum_exp(sum,val);
            }
            alpha[j][t] = sum + log(hmm->B[j][s]);
        }
    }
}

// ---------------- Backward ----------------
void backward(HMM *hmm, char *seq, int len, double **beta) {
    int N=hmm->n_states;
    for (int i=0;i<N;i++) beta[i][len-1] = 0.0;

    for (int t=len-2;t>=0;t--) {
        int s = symbol_idx(seq[t+1]);
        for (int i=0;i<N;i++) {
            double sum = LOG_ZERO;
            for (int j=0;j<N;j++) {
                double val = log(hmm->A[i][j]) + log(hmm->B[j][s]) + beta[j][t+1];
                sum = log_sum_exp(sum,val);
            }
            beta[i][t] = sum;
        }
    }
}

// ---------------- Baum-Welch ----------------
void baum_welch(HMM *hmm, char *seq, int len, int max_iter) {
    int N=hmm->n_states;
    int T=len;
    double **alpha = malloc(N*sizeof(double*));
    double **beta = malloc(N*sizeof(double*));
    for (int i=0;i<N;i++) {
        alpha[i] = malloc(T*sizeof(double));
        beta[i] = malloc(T*sizeof(double));
    }

    for (int iter=0;iter<max_iter;iter++) {
        forward(hmm, seq, T, alpha);
        backward(hmm, seq, T, beta);

        double *gamma_sum = calloc(N,sizeof(double));
        double **gamma = malloc(N*sizeof(double*));
        for (int i=0;i<N;i++) gamma[i] = calloc(T,sizeof(double));

        // 計算 gamma[i][t] = P(state i at t | seq)
        for (int t=0;t<T;t++) {
            double denom = LOG_ZERO;
            for (int i=0;i<N;i++)
                denom = log_sum_exp(denom, alpha[i][t]+beta[i][t]);
            for (int i=0;i<N;i++) {
                gamma[i][t] = exp(alpha[i][t]+beta[i][t]-denom);
                gamma_sum[i] += gamma[i][t];
            }
        }

        // 更新 pi
        for (int i=0;i<N;i++) hmm->pi[i] = gamma[i][0];

        // 更新 A
        for (int i=0;i<N;i++) {
            double denom = 0.0;
            for (int t=0;t<T-1;t++) denom += gamma[i][t];
            for (int j=0;j<N;j++) {
                double numer=0.0;
                for (int t=0;t<T-1;t++) {
                    int s = symbol_idx(seq[t+1]);
                    numer += gamma[i][t]*hmm->A[i][j]*hmm->B[j][s]; // 近似
                }
                hmm->A[i][j] = numer/denom;
            }
        }

        // 更新 B
        for (int i=0;i<N;i++) {
            for (int k=0;k<N_SYMBOLS;k++) {
                double numer=0.0;
                for (int t=0;t<T;t++)
                    if (symbol_idx(seq[t])==k) numer+=gamma[i][t];
                hmm->B[i][k] = numer/gamma_sum[i];
            }
        }

        for (int i=0;i<N;i++) free(gamma[i]);
        free(gamma);
        free(gamma_sum);
    }

    for (int i=0;i<N;i++) { free(alpha[i]); free(beta[i]); }
    free(alpha); free(beta);
}

// ---------------- log probability ----------------
double log_prob_seq(HMM *hmm, char *seq, int len) {
    int N=hmm->n_states;
    double **alpha = malloc(N*sizeof(double*));
    for (int i=0;i<N;i++) alpha[i]=malloc(len*sizeof(double));
    forward(hmm, seq, len, alpha);
    double logp = LOG_ZERO;
    for (int i=0;i<N;i++) logp = log_sum_exp(logp, alpha[i][len-1]);
    for (int i=0;i<N;i++) free(alpha[i]);
    free(alpha);
    return logp;
}

// ---------------- 主程式 ----------------
int main() {
    const char *fasta_file = "C:/Users/katya/Desktop/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna";

    char *buffer_train = malloc(MAX_LEN);
    char *buffer_test1 = malloc(MAX_LEN);
    char *buffer_test2 = malloc(MAX_LEN);

    long len_train = fetch_sequence_safe(fasta_file,"chr1",100000000,MAX_LEN,buffer_train);
    long len_test1 = fetch_sequence_safe(fasta_file,"chr1",145000000,MAX_LEN,buffer_test1);
    long len_test2 = fetch_sequence_safe(fasta_file,"chr2",100000000,MAX_LEN,buffer_test2);

    if (len_train<=0 || len_test1<=0 || len_test2<=0) {
        printf("Failed to fetch sequences\n");
        return 1;
    }

    HMM *hmm1 = init_hmm(2);
    baum_welch(hmm1, buffer_train, len_train, MAX_ITER);

    HMM *hmm2 = init_hmm(3);
    baum_welch(hmm2, buffer_train, len_train, MAX_ITER);

    printf("HMM1 (2 states) log probabilities:\n");
    printf("  training: %.3f\n", log_prob_seq(hmm1, buffer_train, len_train));
    printf("  test1:    %.3f\n", log_prob_seq(hmm1, buffer_test1, len_test1));
    printf("  test2:    %.3f\n", log_prob_seq(hmm1, buffer_test2, len_test2));

    printf("HMM2 (3 states) log probabilities:\n");
    printf("  training: %.3f\n", log_prob_seq(hmm2, buffer_train, len_train));
    printf("  test1:    %.3f\n", log_prob_seq(hmm2, buffer_test1, len_test1));
    printf("  test2:    %.3f\n", log_prob_seq(hmm2, buffer_test2, len_test2));

    free(buffer_train); free(buffer_test1); free(buffer_test2);
    return 0;
}

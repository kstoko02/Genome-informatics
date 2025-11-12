#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define MAX_LEN 5000000
#define MAX_ORDER 9       // 支援 0~9 階
#define CALC_ORDER 7      // 只計算到 7 階機率
#define ALPHABET "ACGT"
#define ALPHABET_SIZE 4
#define PSEUDO 0.5

// ---------------- 基本函式 ----------------
char to_base(char c) {
    c = toupper(c);
    if (c=='A'||c=='C'||c=='G'||c=='T') return c;
    return 0;
}

char complement(char c) {
    switch(c) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default:  return 0;
    }
}

long pow4(int k) {
    long r = 1;
    for (int i=0;i<k;i++) r*=4;
    return r;
}

typedef struct {
    int order;
    long *counts;
    double *prob;
} MarkovModel;

// ---------------- 索引編碼 ----------------
long idx_from_sequence(char *seq, long pos, int k) {
    long idx = 0;
    for (int i=0;i<=k;i++) {
        int c = to_base(seq[pos - k + i]);
        if (c==0) return -1;
        switch(c){
            case 'A': c=0; break;
            case 'C': c=1; break;
            case 'G': c=2; break;
            case 'T': c=3; break;
        }
        idx = idx * 4 + c;
    }
    return idx;
}

void revcomp_seq(char *seq, char *rev, long len) {
    for (long i = 0; i < len; i++) {
        char c = seq[len - 1 - i];
        switch(c) {
            case 'A': rev[i]='T'; break;
            case 'T': rev[i]='A'; break;
            case 'C': rev[i]='G'; break;
            case 'G': rev[i]='C'; break;
            default:  rev[i]='N'; break;
        }
    }
    rev[len]='\0';
}


// ---------------- 初始化模型 ----------------
MarkovModel init_model(int k) {
    MarkovModel m;
    m.order = k;
    long size = pow4(k+1);
    m.counts = calloc(size,sizeof(long));
    if (k <= CALC_ORDER)
        m.prob   = calloc(size,sizeof(double));
    else
        m.prob = NULL; // 8~9階不計算機率
    if (!m.counts || (k<=CALC_ORDER && !m.prob)) {
        fprintf(stderr,"Memory allocation failed for order %d\n", k);
        exit(1);
    }
    return m;
}

// ---------------- 抓取序列 ----------------
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
            if (base != 0) {
                if (seq_pos >= start_pos && idx < length)
                    buffer[idx++] = base;
            }
            if (idx >= length) break;
        }
        if (idx >= length) break;
    }
    fclose(fp);
    return idx;
}

// ---------------- 同時計算正向 + 反向互補 ----------------
void train_model(MarkovModel *m, char *seq, long len) {
    int k = m->order;

    // 建立反向互補序列
    char *rev = malloc((len + 1) * sizeof(char));
    if (!rev) { fprintf(stderr,"Memory allocation failed in reverse complement\n"); exit(1); }
    for (long i = 0; i < len; i++) {
        rev[len - 1 - i] = complement(seq[i]);
    }
    rev[len] = '\0';

    // 計數正向與反向互補
    for (long i = k; i < len; i++) {
        long idx_f = idx_from_sequence(seq, i, k);   // 正向
        long idx_r = idx_from_sequence(rev, i, k);   // 反向互補
        if (idx_f != -1) m->counts[idx_f]++;
        if (idx_r != -1) m->counts[idx_r]++;
    }

    free(rev);
}

// ---------------- 機率計算 ----------------
// ---------------- 機率計算（依題目公式） ----------------
void compute_probability(MarkovModel *m, long total_kmer_count) {
    if (m->order > CALC_ORDER) return;

    int k = m->order;
    long total_size = pow4(k+1);
    double denom = (PSEUDO * pow4(k+1)) + total_kmer_count;

    for (long idx = 0; idx < total_size; idx++) {
        m->prob[idx] = (PSEUDO + m->counts[idx]) / denom;
    }
}


// ---------------- 列印模型 ----------------
void print_model(MarkovModel *m) {
    if (m->order > CALC_ORDER) return; // 只印 0~7 階
    int k = m->order;
    long total_size = pow4(k+1);
    char map[ALPHABET_SIZE] = {'A','C','G','T'};

    printf("\nOrder %d Markov model (前 16 個 k-mer):\n", k);
    for (long idx=0; idx<total_size && idx<16; idx++) {
        if (m->prob[idx]>0) {
            long tmp = idx;
            char kmer[12];
            for (int i=k;i>=0;i--) {
                kmer[i] = map[tmp%4];
                tmp/=4;
            }
            kmer[k+1]='\0';
            printf("%s: %.4f\n", kmer, m->prob[idx]);
        }
    }
}

// ---------------- 計算單序列的 log probability ----------------

// 計算單序列 log probability
double log_probability_single(MarkovModel *m, char *seq, long len) {
    int k = m->order;
    long total_size = pow4(k+1);
    double logp = 0.0;

    if (len <= k) {
        long idx = idx_from_sequence(seq, len-1, len-1);
        if (idx == -1) return 0.0;
        double p = (m->prob && idx < total_size) ? m->prob[idx] : 1.0 / total_size;
        return log(p);
    }

    for (long i = k; i < len; i++) {
        long idx = idx_from_sequence(seq, i, k);
        if (idx == -1) {
            double denom = (PSEUDO * pow4(k+1)) + (len - k);
            double p = PSEUDO / denom;
            logp += log(p);
        } else {
            double p = (m->prob && idx < total_size) ? m->prob[idx] : 1.0 / total_size;
            logp += log(p);
        }
    }
    return logp;
}

// 計算雙股 log probability（取正向與反向互補最大值）
double log_probability_ds(MarkovModel *m, char *seq, long len) {
    char *rev = malloc((len + 1) * sizeof(char));
    if (!rev) { fprintf(stderr,"Memory allocation failed for revcomp\n"); exit(1); }

    revcomp_seq(seq, rev, len);

    double lp_fwd = log_probability_single(m, seq, len);
    double lp_rev = log_probability_single(m, rev, len);

    free(rev);
    return fmax(lp_fwd, lp_rev);
}



// ---------------- 主程式 ----------------
int main() {
    char *seq_train = malloc(MAX_LEN * sizeof(char));
    char *seq_test1 = malloc(MAX_LEN * sizeof(char));
    char *seq_test2 = malloc(MAX_LEN * sizeof(char));
    if (!seq_train || !seq_test1 || !seq_test2) {
        fprintf(stderr,"Memory allocation failed\n");
        return 1;
    }

    // 抓取 train
    long len_train = fetch_sequence_safe(
        "C:/Users/katya/Desktop/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        "chr1", 100000000, MAX_LEN, seq_train);
    if (len_train <= 0) { 
        free(seq_train); free(seq_test1); free(seq_test2);
        return 1; 
    }
    seq_train[len_train] = '\0';

    // 抓取 test1
    long len_test1 = fetch_sequence_safe(
        "C:/Users/katya/Desktop/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        "chr1", 145000000, MAX_LEN, seq_test1);
    if (len_test1 <= 0) { 
        free(seq_train); free(seq_test1); free(seq_test2);
        return 1; 
    }
    seq_test1[len_test1] = '\0';

    // 抓取 test2
    long len_test2 = fetch_sequence_safe(
        "C:/Users/katya/Desktop/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        "chr2", 100000000, MAX_LEN, seq_test2);
    if (len_test2 <= 0) { 
        free(seq_train); free(seq_test1); free(seq_test2);
        return 1; 
    }
    seq_test2[len_test2] = '\0';

    // 初始化模型
    MarkovModel models[MAX_ORDER+1];
    for (int k=0; k<=MAX_ORDER; k++) {
        models[k] = init_model(k);
        train_model(&models[k], seq_train, len_train);
        long total_k = len_train - k;
        compute_probability(&models[k], total_k);
    }

    // 計算 log probability (double-stranded)
    for (int k=0; k<=CALC_ORDER; k++) {
        double lp_train = log_probability_ds(&models[k], seq_train, len_train);
        double lp_test1 = log_probability_ds(&models[k], seq_test1, len_test1);
        double lp_test2 = log_probability_ds(&models[k], seq_test2, len_test2);

        printf("\nOrder %d log probabilities (double-stranded):\n", k);
        printf("  training: %.3f\n", lp_train);
        printf("  test1:    %.3f\n", lp_test1);
        printf("  test2:    %.3f\n", lp_test2);
    }

    // 釋放記憶體
    for (int k=0; k<=MAX_ORDER; k++) {
        free(models[k].counts);
        if (models[k].prob) free(models[k].prob);
    }
    free(seq_train);
    free(seq_test1);
    free(seq_test2);

    return 0;
}



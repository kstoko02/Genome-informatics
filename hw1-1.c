#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#define PSEUDO 0.5
#define BASES 4
#define MAX_ORDER 9
#define MAX_LEN 5000000

// ---------------- Base utilities ----------------
char to_base(char c){ 
    c=toupper(c); 
    return (c=='A'||c=='C'||c=='G'||c=='T')?c:0; 
}
int base_index(char b){ 
    switch(b){ 
        case 'A':return 0; 
        case 'C':return 1; 
        case 'G':return 2; 
        case 'T':return 3; 
        default:return -1; 
    } 
}
char complement(char b){ 
    switch(b){ 
        case 'A':return 'T'; 
        case 'T':return 'A'; 
        case 'C':return 'G'; 
        case 'G':return 'C'; 
        default:return 'N'; 
    } 
}

// ---------------- Fetch sequence ----------------
long fetch_sequence_safe(const char *filename, char *chrom, long start_pos, long length, char *buffer){
    FILE *fp=fopen(filename,"r"); if(!fp){ perror("Cannot open"); return -1; }
    char line[1024]; int in_chrom=0; long seq_pos=0, idx=0;
    while(fgets(line,sizeof(line),fp)){
        if(line[0]=='>'){ in_chrom=(strncmp(line+1,chrom,strlen(chrom))==0); continue; }
        if(!in_chrom) continue;
        for(int i=0; line[i]&&line[i]!='\n'; i++){
            char b=to_base(line[i]); seq_pos++;
            if(b && seq_pos>=start_pos && idx<length) buffer[idx++]=b;
            if(idx>=length) break;
        }
        if(idx>=length) break;
    }
    fclose(fp);
    return idx;
}

// ---------------- k-mer encoding ----------------
unsigned long long kmer_index(const char *seq,int k){
    unsigned long long idx=0;
    for(int i=0;i<k;i++){
        int b=base_index(seq[i]); if(b<0) return (unsigned long long)-1;
        idx=(idx<<2)|b;
    }
    return idx;
}

// ---------------- Build double strand ----------------
char* build_double_strand(const char *seq,long len){
    char *double_seq=malloc(len*2+1);
    memcpy(double_seq,seq,len);
    for(long i=0;i<len;i++) double_seq[len+i]=complement(seq[len-1-i]);
    double_seq[len*2]='\0';
    return double_seq;
}

// ---------------- Markov model struct ----------------
typedef struct{
    int order;
    unsigned long long *count_kmer;
    unsigned long long *count_prev;
    unsigned long long n_kmers;
    unsigned long long n_prev;
    unsigned long long total_kmers;
    unsigned long long total_prev;
} MarkovModel;

// ---------------- Train Markov model ----------------
MarkovModel* train_markov(const char *seq,long len,int order){
    int k=order+1;
    unsigned long long n_kmers=1ULL<<(2*k);
    unsigned long long n_prev=1ULL<<(2*order);
    unsigned long long *count_kmer=calloc(n_kmers,sizeof(unsigned long long));
    unsigned long long *count_prev=calloc(n_prev,sizeof(unsigned long long));
    if(!count_kmer || !count_prev){ fprintf(stderr,"malloc failed\n"); exit(1); }

    unsigned long long total_kmers=0,total_prev=0;
    for(long i=0;i+k<=len;i++){ unsigned long long idx=kmer_index(&seq[i],k); if(idx!=(unsigned long long)-1){ count_kmer[idx]++; total_kmers++; } }
    for(long i=0;i+order<=len;i++){ unsigned long long idx=kmer_index(&seq[i],order); if(idx!=(unsigned long long)-1){ count_prev[idx]++; total_prev++; } }

    MarkovModel *mm=malloc(sizeof(MarkovModel));
    mm->order=order; mm->count_kmer=count_kmer; mm->count_prev=count_prev;
    mm->n_kmers=n_kmers; mm->n_prev=n_prev;
    mm->total_kmers=total_kmers; mm->total_prev=total_prev;
    return mm;
}

// ---------------- Compute logP using trained model ----------------
double logp_markov(const MarkovModel *mm,const char *seq,long len){
    int k=mm->order+1; double logp=0.0;
    for(long i=0;i+k<=len;i++){
        unsigned long long idx=kmer_index(&seq[i],k);
        unsigned long long idx_prev=kmer_index(&seq[i],mm->order);
        if(idx==(unsigned long long)-1 || idx_prev==(unsigned long long)-1) continue;
        double p_k=(mm->count_kmer[idx]+PSEUDO)/(mm->total_kmers+PSEUDO*mm->n_kmers);
        double p_prev=(mm->count_prev[idx_prev]+PSEUDO)/(mm->total_prev+PSEUDO*mm->n_prev);
        logp+=log(p_k)-log(p_prev);
    }
    return logp;
}

// ---------------- Main ----------------
int main(){
    char *seq_train=malloc(MAX_LEN+1), *seq_test1=malloc(MAX_LEN+1), *seq_test2=malloc(MAX_LEN+1);
    if(!seq_train || !seq_test1 || !seq_test2){ fprintf(stderr,"malloc failed\n"); return 1; }

    long len_train=fetch_sequence_safe("C:/Users/katya/Desktop/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna","chr1",100000000,MAX_LEN,seq_train);
    long len_test1=fetch_sequence_safe("C:/Users/katya/Desktop/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna","chr1",145000000,MAX_LEN,seq_test1);
    long len_test2=fetch_sequence_safe("C:/Users/katya/Desktop/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna","chr2",100000000,MAX_LEN,seq_test2);

    char *train_double=build_double_strand(seq_train,len_train);
    char *test1_double=build_double_strand(seq_test1,len_test1);
    char *test2_double=build_double_strand(seq_test2,len_test2);

    for(int order=0;order<=MAX_ORDER;order++){
        MarkovModel *mm=train_markov(train_double,len_train*2,order);
        double logp_train=logp_markov(mm,train_double,len_train*2);
        double logp_test1=logp_markov(mm,test1_double,len_test1*2);
        double logp_test2=logp_markov(mm,test2_double,len_test2*2);
        printf("Order %d: logP(train)=%.6f, logP(test1)=%.6f, logP(test2)=%.6f\n",order,logp_train,logp_test1,logp_test2);
        free(mm->count_kmer); free(mm->count_prev); free(mm);
    }

    free(seq_train); free(seq_test1); free(seq_test2);
    free(train_double); free(test1_double); free(test2_double);
    return 0;
}

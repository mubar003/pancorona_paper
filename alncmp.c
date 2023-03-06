#define VERSION "alncmp v0.8"
#define USAGE "Usage: alncmp x.lin.msa y.lin.msa outputs.txt [<cutoff>] [LENIENT] [M <matchscore>]\n"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

int main(int argc, char* argv[]) {
    // Parse input parameters, validate files
    if (argc < 4) {puts(VERSION"\n"USAGE); exit(1);}
    int lenient = 0; double cutoff = 10., match = 5;
    if (argc > 5 && !strcmp(argv[argc-2],"M")) match = atof(argv[argc-1]), printf("Match %.4f\n",match), argc-=2;
    if (argc > 4 && !strcmp(argv[argc-1],"LENIENT")) puts("LENIENT mode."),lenient = 1, --argc;
    if (argc > 4) cutoff = atof(argv[argc-1]);
    FILE *f1 = fopen(argv[1],"rb");
    FILE *f2 = fopen(argv[2],"rb");
    FILE *outf = fopen(argv[3],"wb");
    if (!f1 || !f2 || !outf) {puts("ERROR in file handling"); exit(2);}

    // Allocate and read in sequences, validate lengths
    char *s1 = calloc(INT32_MAX,1), *s2 = calloc(INT32_MAX,1);
    fgets(s1,INT32_MAX,f1); fgets(s2,INT32_MAX,f2); // ignore header
    fgets(s1,INT32_MAX,f1); fgets(s2,INT32_MAX,f2);
    int len = strlen(s1), len2=strlen(s2);
    if (len != len2) {printf("Error: MSAs are not equal length: %d and %d\n",len,len2); exit(4);}

    // Convert sequences to uniform 32-based numerical format
    for (int i = 0; i < len; ++i) // Cast all letters (upper and lower) to 32 bins
        s1[i] = s1[i] >= 64 ? s1[i] & 31 : 31,  // cast letters, treat non-letters as gaps (31)
        s2[i] = s2[i] >= 64 ? s2[i] & 31 : 31, 
        s1[i] = s1[i] == 21 ? 20 : s1[i], // cast all U to T in seq. 1
        s2[i] = s2[i] == 21 ? 20 : s2[i], // cast all U to T in seq. 2
        s1[i] = s1[i] == 24 ? 14 : s1[i], // Cast X's to N's
        s2[i] = s2[i] == 24 ? 14 : s2[i]; // Cast X's to N's

    // Define scoring criteria using IUPAC convention
    typedef struct { char L, *M; } S_t;
    S_t S[32] = {0};
    S['A' & 31] = (S_t){1,(char[]){'A'&31}}; // unambiguous
    S['C' & 31] = (S_t){1,(char[]){'C'&31}}; 
    S['G' & 31] = (S_t){1,(char[]){'G'&31}}; 
    S['T' & 31] = (S_t){1,(char[]){'T'&31}}; 

    S['R'&31] = (S_t){2,(char[]){'A'&31,'G'&31}}; // dual ambigs
    S['Y'&31] = (S_t){2,(char[]){'C'&31,'T'&31}}; 
    S['S'&31] = (S_t){2,(char[]){'G'&31,'C'&31}}; 
    S['W'&31] = (S_t){2,(char[]){'A'&31,'T'&31}}; 
    S['K'&31] = (S_t){2,(char[]){'G'&31,'T'&31}}; 
    S['M'&31] = (S_t){2,(char[]){'A'&31,'C'&31}}; 

    S['B'&31] = (S_t){3,(char[]){'C'&31,'G'&31,'T'&31}}; // triple ambigs
    S['D'&31] = (S_t){3,(char[]){'A'&31,'G'&31,'T'&31}}; 
    S['H'&31] = (S_t){3,(char[]){'A'&31,'C'&31,'T'&31}}; 
    S['V'&31] = (S_t){3,(char[]){'A'&31,'C'&31,'G'&31}}; 
    
    S['N'&31] = (S_t){4,(char[]){'A'&31,'C'&31,'G'&31,'T'&31}}; // total ambigs
    //double gapP = -3;

    // Compute scoring matrix using IUPAC set theory. Score ~ 1 / max(set)
    double SM[32][32] = {0.};
    for (int i = 0; i < 32; ++i) {
        for (int j = i; j < 32; ++j) {
            double max_ix = S[i].L >= S[j].L ? S[i].L : S[j].L;
            if (!max_ix) continue; // no overlaps
            int ov = 0; // count the overlaps
            for (int L1 = 0; L1 < S[i].L; ++L1)
                for (int L2 = 0; L2 < S[j].L; ++L2) 
                    ov += S[i].M[L1] == S[j].M[L2];
            // Apply final score
            //if (i==31 || j==31) SM[i][j] = SM[j][i] = gapP;
            //else 
            if (lenient) SM[i][j] = SM[j][i] = (double)(ov)/max_ix;
            else SM[i][j] = SM[j][i] = (double)(ov>0)/max_ix;
        }
    }

    //SM[31][31] = 0; // don't penalize all-gap case

    // Debug: print the matrix
    /* for (int i = 0; i < 32; ++i) fprintf(outf,"\t%c",i+64);
    fprintf(outf,"\n");
    for (int i = 0; i < 32; ++i) {
        fprintf(outf,"%c",i+64);
        for (int j = 0; j < 32; ++j) {
            fprintf(outf,"\t%.4f",SM[i][j]);
        }
        fprintf(outf,"\n");
    } */

    // Tally the scores
    double score = 0, score_raw = 0, maxScore = 0, maxScore_raw = 0, global = 0, global_raw = 0; 
    int L = 0, skp = 0, gap = 0, max_ix = -1, max_L = 0, max_skp = 0, max_gap = 0,
        skp_global = 0, gap_global = 0;

    double mismatch = -match + 1, gapP = -1; // To adjust scores: 2*S - 1

    int stretch_st = 0, stretch_ed = 0, stretch_len = 0; double stretch_sc = 0;
    int longest_st = 0, longest_ed = 0, longest_len = 0; double longest_sc = 0;
    
    char *str = malloc(len+1), *str2 = malloc(len+1);
    fprintf(outf,"StartPos\tEndPos\tSeq1\tSeq2\tAlnLen\tScore\tpcntID\n");
    for (int i = 0; i < len; ++i) {
        double s = SM[s1[i]][s2[i]]; 
        int doskp = s1[i] == 31 && s2[i] == 31;
        int doGap = !doskp && (s1[i] == 31 || s2[i] == 31);
        score_raw += s; // the "raw" match score, used for ID calculation
        score += doskp? 0 : s*match + mismatch; 
        score += doGap ? gapP : 0; // the "adjusted" penalty used for scoring
        skp += doskp; gap += doGap; ++L;
        
        //printf("Pos %d: Score: %f, thisscore: %f\n",i,score,s); // Debug

        // Tally the longest positive stretch as an alternative report
        if (s > 0 || doskp || i == len - 1) {
            stretch_sc += s;
            stretch_ed = i;
            stretch_len += !doskp;
            if (stretch_sc > longest_sc) longest_sc = stretch_sc, longest_ed = i, 
                longest_len = stretch_len, longest_st = stretch_st; 
        } else  // streak interrupted
            stretch_st = i+1, stretch_len = 0, stretch_sc = 0;
        
        if (score > maxScore)  // update max score
            maxScore = score, maxScore_raw = score_raw,
            max_ix = i, max_L = L, max_skp = skp, max_gap = gap;
        
        if (score < cutoff || i == len-1) { // reset tally, write last max (if valid)
            if (maxScore >= cutoff) {
                int stix = 0; for (int k = max_ix-max_L+1; k <= max_ix; ++k) 
                    if (s1[k] != 31 || s2[k] != 31) {
                        int off = s1[k] == 31 || s2[k] == 31 || s1[k] == s2[k] ? 64 : 96;
                        str[stix] = s1[k]+off, str2[stix] = s2[k]+off, ++stix;
                    }
                str[stix] = str2[stix] = 0;
                fprintf(outf,"%d\t%d\t%s\t%s\t%d\t%f\t%f\n",max_ix-max_L+2,max_ix+1,str,str2,max_L-max_skp,
                    maxScore,maxScore_raw/(max_L-max_skp)*100);
                score = 0;
            }
            if (score <= 0) score = score_raw = maxScore = maxScore_raw = 0, L = skp = gap = 0;
        }
        global += doskp ? 0 : s*match + mismatch; 
        global_raw += s;
        skp_global += doskp;
        gap_global += doGap;
    }

    printf("Len=%d, Global score = %f [%.2f %%]; all-gap = %d\n",len,global, global_raw/(len-skp_global)*100,skp_global);
    printf("Longest contiguous chain: ID %.4f, st %d, ed %d, L=%d\n",longest_sc/longest_len*100, longest_st, longest_ed, longest_len);

    return 0;
}
/*
 * ml_fasta.c
 * 
 * ML: MinLookup
 * 1. Compute MinHash for each read
 * 2. Create dictionary of minHash -> kvec<readNum,pos>
 * 3. Output target reads with > threshold matches for each query
 *
 * Jeremy Wang
 * 20150624
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include "../klib/kseq.h" // FASTA/Q parser
#include "../klib/kvec.h" // C dynamic vector
#include "../klib/khash.h" // C hash table/dictionary

typedef struct {
  uint32_t hash;
  uint32_t pos;
} Min;

// MAXIMUM # READS = 2^31 (~2b)
typedef struct {
  uint32_t readNum; // last bit 0 if fw, 1 if rv (same as alternating fw/rv) [fw read n = n*2, rv read n = n*2+1]
  // in practice, we only hash the fw strand, and do lookup by fw and rv, so right now these will always end in 0 (readNum%2 == 0)
  uint32_t pos;
} readMin;

typedef struct {
  uint32_t qpos;
  uint32_t tpos;
} posPair;

typedef kvec_t(readMin) matchVec;

typedef kvec_t(posPair) offsetVec;

typedef kvec_t(Min*) minVec;

// init kseq struct
KSEQ_INIT(FILE*, fileread);

// creates uint32:kvec<readNum,pos> hash
KHASH_MAP_INIT_INT(minDict, matchVec);

// creates uint32:kvec<posPair> hash
KHASH_MAP_INIT_INT(overlapDict, offsetVec);

Min* minhash (char *s, int k, int h, uint32_t hash_seeds[], unsigned char reverse) {

  int i, j;
  int slen = strlen(s);

  // initialize minimums, one per hash seed
  Min* m = (Min*)malloc(sizeof(Min)*h);

  for(j = 0; j < h; j++) {
    m[j].hash = UINT32_MAX;
  }

  // set up kmer uint32
  uint32_t kmer = 0;
  int rev_shift = 2*(k-1);
  for(i = 0; i < k-1; i++) {
    if(reverse == 0) {
      kmer = (kmer << 2) + ((s[i] >> 1) & 3);
    } else {
      kmer = (kmer >> 2) + ((((s[i] >> 1) & 3) ^ 2) << rev_shift);
    }
  }
  for(i = 0; i <= slen-k; i++) {

    // just use uint32 kmer representation, with random xor
    if(reverse == 0) {
      kmer = (kmer << 2) + ((s[i+k-1] >> 1) & 3);
      if(k < 16) {
        kmer = kmer % (1 << (2*k));
      }
    } else {
      kmer = (kmer >> 2) + ((((s[i+k-1] >> 1) & 3) ^ 2) << rev_shift);
    }

    for(j = 0; j < h; j++) {
      kmer = kmer ^ hash_seeds[j];
      if(kmer < m[j].hash) {
        m[j].hash = kmer;
        m[j].pos = i;
      }
    }
  }

  return m;
}

// if reverse is true (1), forward and reverse signatures will be adjacent such that
// signature of read X forward is at index (2*X), and reverse is (2*X + 1)
void hash_fasta_signatures(char* fa, int k, int h, uint32_t hash_seeds[], khash_t(minDict) **mh_lookups, minVec *mh, int readLimit) {

  FILE* fp;
  kseq_t* seq;
  int l, i, ret_val;

  fp = fopen(fa, "r");
  seq = kseq_init(fp);
  //printf("Reading fasta file: %s\n", fa);

  khint_t bin; // hash bin (result of kh_put)

  uint32_t readNum = 0;
  while ((l = kseq_read(seq)) >= 0) {
    // name: seq->name.s, seq: seq->seq.s, length: l

    Min *m = minhash(seq->seq.s, k, h, hash_seeds, 0);
    if(mh_lookups != NULL) {
      for(i = 0; i < h; i++) {
        bin = kh_put(minDict, mh_lookups[i], m[i].hash, &ret_val);
        if(ret_val == 1) { // bin is empty (unset)
          kv_init(kh_value(mh_lookups[i], bin)); // kh_value should already be defined as type matchVec
        }
        readMin r;
        r.readNum = (readNum << 1);
        r.pos = m[i].pos;
        kv_push(readMin, kh_value(mh_lookups[i], bin), r);
      }
    }
    if(mh != NULL) {
      kv_push(Min*, *mh, m);
      m = minhash(seq->seq.s, k, h, hash_seeds, 1);
      kv_push(Min*, *mh, m);
    } else {
      free(m);
    }
    readNum++;

    if(readLimit > 0 && readNum >= readLimit) {
      break;
    }
  }

  kseq_destroy(seq);
  fclose(fp);
}

// have to reorder params to make this work with kseq
// so stupid
int fileread(FILE* f, char* buffer, int size) {
  return fread(buffer, 1, size, f);
}

int ml_fasta(char *query_fa, char *target_fa, int k, int h, int seed, int threshold, int max_kmers, int readLimit) {

  // construct array of hash seeds
  srand(seed); // seed random number generator
  // RAND_MAX must be at least 2**15, mine is 2**31
  // as long as we don't duplicate seeds, it should be fine
  uint32_t hash_seeds[h];
  int i;
  for(i = 0; i < h; i++) {
    hash_seeds[i] = rand(); // 0 to RAND_MAX
  }

  // ------------------------- Create target read hash -----------------------------

  time_t t0 = time(NULL);

  // read FASTA file, construct hash, including only forward direction
  // -- maybe assess doing this the opposite way at a later date, I'm not sure which will be faster
  khash_t(minDict) *mh_lookups[h];
  for(i = 0; i < h; i++) {
    mh_lookups[i] = kh_init(minDict); // allocate hash table
  }
  minVec mh;
  kv_init(mh);
  if(strcmp(query_fa, target_fa) == 0) { // only one input fasta
    hash_fasta_signatures(query_fa, k, h, hash_seeds, mh_lookups, &mh, readLimit); // load both lookup and hashes from query
  } else {
    hash_fasta_signatures(target_fa, k, h, hash_seeds, mh_lookups, NULL, readLimit); // load only lookup from target
    hash_fasta_signatures(query_fa, k, h, hash_seeds, NULL, &mh, readLimit); // load only hashes from query
  }

  time_t t1 = time(NULL);
  printf("# Loaded and processed %d target reads in %d seconds\n", mh.n/2, (t1-t0));
  t0 = t1;
  // -------------------------------------------------------------------------------

  // ---------------------------- Look up hashes in dictionary ------------------------------
  int q, qidx, qrev, m, ret_val;
  uint32_t target;
  khint_t bin; // hash bin (result of kh_put, kh_get)
  Min *qmin;
  matchVec matches;
  khash_t(overlapDict) *overlaps;
  for(q = 0; q < mh.n/2; q++) {
    for(qrev = 0; qrev <= 1; qrev++) {
      if(qrev == 1) {
        qidx = q*2 + 1;
      } else {
        qidx = q*2;
      }
      qmin = (Min*)(mh.a[qidx]);
      overlaps = kh_init(overlapDict); // allocate hash table
      for(i = 0; i < h; i++) {
        bin = kh_get(minDict, mh_lookups[i], qmin[i].hash);
        if(bin == kh_end(mh_lookups[i])) // key not found, IDK what happens if you don't test this
          continue;
        matches = kh_val(mh_lookups[i], bin);
        if(matches.n > max_kmers) { // repetitive, ignore it
          continue;
        }
        for(m = 0; m < matches.n; m++) {
          bin = kh_put(overlapDict, overlaps, matches.a[m].readNum>>1, &ret_val); // >>1 removes the fw/rv bit, which is always fw(0) right now
          if(ret_val == 1) { // bin is empty (unset)
            kv_init(kh_value(overlaps, bin)); // kh_value should already be defined as type matchVec
          }
          posPair ppair; // to store matching query/target positions
          ppair.qpos = qmin[i].pos;
          ppair.tpos = matches.a[m].pos;
          kv_push(posPair, kh_value(overlaps, bin), ppair);
        }
      }
      // count hits and compute offset for each target
      offsetVec offsets;
      khint_t iter;
      // from the definition of kh_foreach
      for (iter = kh_begin(overlaps); iter != kh_end(overlaps); ++iter) {
        if (!kh_exist(overlaps, iter)) continue;
        target = kh_key(overlaps, iter);
        offsets = kh_val(overlaps, iter);
        if(offsets.n >= threshold) {
          /*
           * To compute offset directly:
           *
          offset = 0;
          for(i = 0; i < offsets.n; i++) {
            offset = offset + (offset.a[i].tpos - offsets.a[i].qpos);
          }
          offset = offset / (int)offsets.n; // offsets.n is a size_t, does not play nicely
           *
           */
          printf("%d,%d,%d,%d", q, qrev, target, offsets.n);
          for(i = 0; i < offsets.n; i++)
            printf(",%d:%d", offsets.a[i].qpos, offsets.a[i].tpos);
          printf("\n");
        }
        kv_destroy(offsets);
      }
      kh_destroy(overlapDict, overlaps);
    }
  }

  t1 = time(NULL);
  printf("# Compared and output in %d seconds\n", (t1-t0));
  // ----------------------------------------------------------------------------------------

  for(i = 0; i < h; i++)
    kh_destroy(minDict, mh_lookups[i]);
  kv_destroy(mh);

  return 0;
}

int main(int argc, char *argv[]) {
  if(argc < 8) {
    printf("Usage: ml_fasta <query_fasta> <target_fasta> <k> <h> <seed> <threshold> <max_kmer_hits>\n");
    printf("  query_fasta: FASTA file containing reads\n");
    printf("  target_fasta: FASTA file containing reads\n");
    printf("  k: Size of k-mer to hash\n");
    printf("  h: Number of hash functions to apply\n");
    printf("  seed: Seed to random number generator\n");
    printf("  threshold: Minimum number of k-mers to declare a match\n");
    printf("  max_kmer_hits: Maximum occurrences of a k-mer before it is considered repetitive and ignored\n");
    printf("Not enough arguments.\n");
    return -1;
  }
  char *query_fasta = argv[1];
  char *target_fasta = argv[2];
  int k = atoi(argv[3]);
  if(k > 16) {
    printf("Current optimizations do not allow for k > 16");
    return -1;
  }
  int h = atoi(argv[4]);
  int seed = atoi(argv[5]);
  int threshold = atoi(argv[6]);
  int max_kmers = atoi(argv[7]);
  //printf("PROCESSING ONLY FIRST 1000 READS\n");
  return ml_fasta(query_fasta, target_fasta, k, h, seed, threshold, max_kmers, -1);
}

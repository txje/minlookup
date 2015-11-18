MinLookup
=======

Fast memory-efficient pre-alignment of long high-error sequences using locality-sensitive hashing


Installation
-------

    mkdir src
    git clone https://github.com/attractivechaos/klib.git klib
    git clone https://github.com/txje/minlookup.git src
    gcc src/ml_fasta.c -o minlookup


Usage
-----

Command line:

    minlookup <query_fasta> <target_fasta> <k> <h> <seed> <threshold> <max_kmer_count> <output>


ml_fasta.c
----------

Parameters: &lt;query_fasta> &lt;target_fasta> &lt;k> &lt;h> &lt;seed> &lt;threshold>

Query and target files may be the same, in which case it will do pairwise comparison

Results should be piped to an output file, and are in the form:
&lt;query_idx>,&lt;query_reverse?>,&lt;target_idx>,&lt;# matches>,&lt;qpos:tpos>,&lt;qpos:tpos>,&lt;qpos:tpos>,...

Where query and target are indices into their respective FASTA, query_reverse is 1(reverse complemented) or 0, # matches is the integer total of matched k-mers (out of a maximum of h), and the offset is the computed average offset between matched k-mers.

max_kmer_count is used to ignore minimum k-mer values with too many occurrences, ex. repetitive sequences, and should be on the order of your expected coverage

Values >16 for k are not allowed because I'm packing the k-mer into a uint32. This may change in the future.

Match positions reported on the reverse strand are the leftmost index of the matched k-mer relative to the forward strand

system "perl translate.pl";
system "perl hmm.pl > hmm_alignment.txt";
system "perl process_hmm.pl";
system "perl cluster_template.pl";


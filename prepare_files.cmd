For novel noncoding:
./check.pl -n -p alu > novel_noncoding_Alu.tab
cut -f1-5,9,10,12,14,16,18 novel_noncoding_Alu.tab > noncoding_novel_Alu.cut.tab
./check.pl -n -p line > novel_noncoding_LINE.tab
cut -f1-5,9,10,12,14,16,18 novel_noncoding_LINE.tab > noncoding_novel_LINE.cut.ta

For novel coding:
./check.pl -p alu > novel_coding_Alu.tab
cut -f1-5,9,11,13,15,17,19 novel_coding_Alu.tab > coding_novel_Alu.cut.tab


scp *.cut.tab salz:/ccb/salz4-1/gpertea/gtex_repeats/

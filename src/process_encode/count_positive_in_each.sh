# commands for counting positives

for i in $(seq 1 32);do cut calls_train.txt -f $i | grep -c 1 >> calls_train_npos.txt; done
for i in $(seq 1 32);do cut calls_valid.txt -f $i | grep -c 1 >> calls_valid_npos.txt; done
for i in $(seq 1 32);do cut calls_test.txt -f $i | grep -c 1 >> calls_test_npos.txt; done

echo -e 'cell\ttrain\tvalid\ttest' > positive_counts.txt
paste <(head -n 1 calls_train.txt | tr '\t' '\n') calls_train_npos.txt calls_valid_npos.txt calls_test_npos.txt >> positive_counts.txt
paste <(echo 'total_entries') <(tail -n +2 calls_train.txt | wc -l) <(tail -n +2 calls_valid.txt | wc -l) <(tail -n +2 calls_test.txt | wc -l) >> positive_counts.txt

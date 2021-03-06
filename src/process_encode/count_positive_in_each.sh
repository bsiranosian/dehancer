# commands for counting positives

#in the all calls set
cd /users/bsiranos/analysis/enhancer_conservation/encode_data/broad/window_calls
for i in $(seq 1 32);do cut all.calls -f $i | grep -c 1 >> calls_npos.txt; done
for i in $(seq 1 32);do cut all.calls -f $i | grep -c 0 >> calls_nneg.txt; done


# in training, validation and test data
cd /users/bsiranos/analysis/enhancer_conservation/encode_data/broad/ml_train_valid_test/
rm calls_train_npos.txt calls_valid_npos.txt calls_test_npos.txt calls_train_nneg.txt calls_valid_nneg.txt calls_test_nneg.txt 
for i in $(seq 1 32);do tail -n +2 calls_train.txt | cut -f $i | grep -c 1 >> calls_train_npos.txt; done
for i in $(seq 1 32);do tail -n +2 calls_valid.txt | cut -f $i | grep -c 1 >> calls_valid_npos.txt; done
for i in $(seq 1 32);do tail -n +2 calls_test.txt | cut -f $i | grep -c 1 >> calls_test_npos.txt; done
for i in $(seq 1 32);do tail -n +2 calls_train.txt | cut -f $i | grep -c 0 >> calls_train_nneg.txt; done
for i in $(seq 1 32);do tail -n +2 calls_valid.txt | cut -f $i | grep -c 0 >> calls_valid_nneg.txt; done
for i in $(seq 1 32);do tail -n +2 calls_test.txt | cut -f $i | grep -c 0 >> calls_test_nneg.txt; done

echo -e 'cell\ttrain\tvalid\ttest' > positive_counts.txt
paste <(head -n 1 calls_train.txt | tr '\t' '\n') calls_train_npos.txt calls_valid_npos.txt calls_test_npos.txt >> positive_counts.txt
echo -e 'cell\ttrain\tvalid\ttest' > negative_counts.txt
paste <(head -n 1 calls_train.txt | tr '\t' '\n') calls_train_nneg.txt calls_valid_nneg.txt calls_test_nneg.txt >> negative_counts.txt
# paste <(echo 'total_entries') <(tail -n +2 calls_train.txt | wc -l) <(tail -n +2 calls_valid.txt | wc -l) <(tail -n +2 calls_test.txt | wc -l) >> positive_counts.txt

# get weights from number of positive, negative compared to total
paste <(cut -f2 positive_counts.txt) <(cut -f2 negative_counts.txt) | tail -n +2 | awk '{print 1 /( $1 / ($1 + $2))}' > positive_call_weights.txt
paste <(cut -f2 negative_counts.txt) <(cut -f2 positive_counts.txt) | tail -n +2 | awk '{print 1 /( $1 / ($1 + $2))}' > negative_call_weights.txt

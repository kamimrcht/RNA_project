

# 9 100 41 36.59 96.13



echo "small_k	nb_small_k	nb_windows	recall	precision"> result_bench.txt
echo "small_k	nb_small_k	nb_windows	recall	precision"> result_genes_bench.txt
#~ small_k=13
#~ nb_small_k=1130
#~ nb_windows=19

#~ ~/RNA_project/long_read_connector.sh -f ~/RNA_project/tools/long_read_connector/scripts/read_simulator/simulatedReads.fa
						#~ ~/RNA_project/tools/long_read_connector/scripts/results_clustering.py long_read_connector_res.txt ~/RNA_project/tools/long_read_connector/scripts/read_simulator/simulatedReads.fa > result_tmp.txt
						#~ echo $small_k $nb_small_k $nb_windows $(cat result_tmp.txt) >> result_bench.txt
						#~ echo $small_k $nb_small_k $nb_windows $(cat /home/marchet/RNA_project/expe/global_results.txt) >> result_genes_bench.txt
#~ ~/RNA_project/long_read_connector.sh -f ~/RNA_project/expe/test.fa
						#~ ~/RNA_project/tools/long_read_connector/scripts/results_clustering.py long_read_connector_res.txt ~/RNA_project/tools/long_read_connector/scripts/read_simulator/simulatedReads.fa > result_tmp.txt -K ${small_k} -N ${nb_small_k} -W ${nb_windows}
						#~ echo $small_k $nb_small_k $nb_windows $(cat result_tmp.txt) >> result_bench.txt
						#~ echo $small_k $nb_small_k $nb_windows $(cat /home/marchet/RNA_project/expe/global_results.txt) >> result_genes_bench.txt


#~ for small_k in 5 7 9 11 13 15
	#~ do
		#~ for nb_windows in 5 10 15 20 25 30
			#~ do
				#~ case $small_k in
					#~ 5)
					
					#~ for nb_small_k in {800..1500..50}
					#~ do
						#~ ~/RNA_project/long_read_connector.sh -f ~/RNA_project/tools/long_read_connector/scripts/read_simulator/simulatedReads.fa -K ${small_k} -N ${nb_small_k} -W ${nb_windows}
						#~ ~/RNA_project/tools/long_read_connector/scripts/results_clustering.py long_read_connector_res.txt ~/RNA_project/tools/long_read_connector/scripts/read_simulator/simulatedReads.fa > result_tmp.txt
						#~ echo $small_k $nb_small_k $nb_windows $(cat result_tmp.txt) >> result_bench.txt
						#~ echo $small_k $nb_small_k $nb_windows $(cat /home/marchet/RNA_project/expe/global_results.txt) >> result_genes_bench.txt
					#~ done
					#~ ;;
					
					#~ 7)
			
					#~ for nb_small_k in {700..1500..50}
					#~ do
						#~ ~/RNA_project/long_read_connector.sh -f ~/RNA_project/tools/long_read_connector/scripts/read_simulator/simulatedReads.fa -K ${small_k} -N ${nb_small_k} -W ${nb_windows}
						#~ ~/RNA_project/tools/long_read_connector/scripts/results_clustering.py long_read_connector_res.txt ~/RNA_project/tools/long_read_connector/scripts/read_simulator/simulatedReads.fa > result_tmp.txt
						#~ echo $small_k $nb_small_k $nb_windows $(cat result_tmp.txt) >> result_bench.txt
						#~ echo $small_k $nb_small_k $nb_windows $(cat /home/marchet/RNA_project/expe/global_results.txt) >> result_genes_bench.txt
					#~ done
					#~ ;;
					
					#~ 9)
			
					#~ for nb_small_k in {700..1300..50}
					#~ do
						#~ ~/RNA_project/long_read_connector.sh -f ~/RNA_project/tools/long_read_connector/scripts/read_simulator/simulatedReads.fa -K ${small_k} -N ${nb_small_k} -W ${nb_windows}
						#~ ~/RNA_project/tools/long_read_connector/scripts/results_clustering.py long_read_connector_res.txt ~/RNA_project/tools/long_read_connector/scripts/read_simulator/simulatedReads.fa > result_tmp.txt
						#~ echo $small_k $nb_small_k $nb_windows $(cat result_tmp.txt) >> result_bench.txt
						#~ echo $small_k $nb_small_k $nb_windows $(cat /home/marchet/RNA_project/expe/global_results.txt) >> result_genes_bench.txt
					#~ done
					#~ ;;
					
					#~ 11)
			
					#~ for nb_small_k in {600..1300..50}
					#~ do
						#~ ~/RNA_project/long_read_connector.sh -f ~/RNA_project/tools/long_read_connector/scripts/read_simulator/simulatedReads.fa -K ${small_k} -N ${nb_small_k} -W ${nb_windows}
						#~ ~/RNA_project/tools/long_read_connector/scripts/results_clustering.py long_read_connector_res.txt ~/RNA_project/tools/long_read_connector/scripts/read_simulator/simulatedReads.fa > result_tmp.txt
						#~ echo $small_k $nb_small_k $nb_windows $(cat result_tmp.txt) >> result_bench.txt
						#~ echo $small_k $nb_small_k $nb_windows $(cat /home/marchet/RNA_project/expe/global_results.txt) >> result_genes_bench.txt
					#~ done
					#~ ;;
					
					#~ 13)
			
					#~ for nb_small_k in {500..1200..50}
					#~ do
						#~ ~/RNA_project/long_read_connector.sh -f ~/RNA_project/tools/long_read_connector/scripts/read_simulator/simulatedReads.fa -K ${small_k} -N ${nb_small_k} -W ${nb_windows}
						#~ ~/RNA_project/tools/long_read_connector/scripts/results_clustering.py long_read_connector_res.txt ~/RNA_project/tools/long_read_connector/scripts/read_simulator/simulatedReads.fa > result_tmp.txt
						#~ echo $small_k $nb_small_k $nb_windows $(cat result_tmp.txt) >> result_bench.txt
						#~ echo $small_k $nb_small_k $nb_windows $(cat /home/marchet/RNA_project/expe/global_results.txt) >> result_genes_bench.txt
					#~ done
					#~ ;;
					
					#~ 15)
			
					#~ for nb_small_k in {500..1000..50}
					#~ do
						#~ ~/RNA_project/long_read_connector.sh -f ~/RNA_project/tools/long_read_connector/scripts/read_simulator/simulatedReads.fa -K ${small_k} -N ${nb_small_k} -W ${nb_windows}
						#~ ~/RNA_project/tools/long_read_connector/scripts/results_clustering.py long_read_connector_res.txt ~/RNA_project/tools/long_read_connector/scripts/read_simulator/simulatedReads.fa > result_tmp.txt
						#~ echo $small_k $nb_small_k $nb_windows $(cat result_tmp.txt) >> result_bench.txt
						#~ echo $small_k $nb_small_k $nb_windows $(cat /home/marchet/RNA_project/expe/global_results.txt) >> result_genes_bench.txt
					#~ done
					#~ ;;
				#~ esac
			#~ done
	#~ done 
#~ rm result_tmp.txt





for small_k in 9 11 13
	do
		for nb_windows in 25 31 41 45 51
			do
				case $small_k in
					9)
					for nb_small_k in {80..200..20}
					do
						~/RNA_project/long_read_connector.sh -f ~/RNA_project/tools/long_read_connector/scripts/read_simulator/simulatedReads.fa -K ${small_k} -N ${nb_small_k} -W ${nb_windows}
						~/RNA_project/tools/long_read_connector/scripts/results_clustering.py long_read_connector_res.txt ~/RNA_project/tools/long_read_connector/scripts/read_simulator/simulatedReads.fa > result_tmp.txt
						echo $small_k $nb_small_k $nb_windows $(cat result_tmp.txt) >> result_bench.txt
						echo $small_k $nb_small_k $nb_windows $(cat /home/marchet/RNA_project/expe/global_results.txt) >> result_genes_bench.txt
					done
					;;

					13)
			
					for nb_small_k in {40..55..5}
					do
						~/RNA_project/long_read_connector.sh -f ~/RNA_project/tools/long_read_connector/scripts/read_simulator/simulatedReads.fa -K ${small_k} -N ${nb_small_k} -W ${nb_windows}
						~/RNA_project/tools/long_read_connector/scripts/results_clustering.py long_read_connector_res.txt ~/RNA_project/tools/long_read_connector/scripts/read_simulator/simulatedReads.fa > result_tmp.txt
						echo $small_k $nb_small_k $nb_windows $(cat result_tmp.txt) >> result_bench.txt
						echo $small_k $nb_small_k $nb_windows $(cat /home/marchet/RNA_project/expe/global_results.txt) >> result_genes_bench.txt
					done
					;;
				esac
			done
	done 
rm result_tmp.txt

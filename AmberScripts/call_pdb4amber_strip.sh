source /home/victormello/Programs/amber24/amber.sh

#pdbs_folder=../data/pdb/RuvB_Monomers/
pdbs_folder=../data/pdb/RuvB_Trimers/
out_folder=../data/amber/pdb4amber/

for str in $(ls $pdbs_folder)
do
	label=$(echo $str | sed 's/.pdb//')
	awk '$6 != 137 && $6 != 138 && $6 != 139 && $6 != 140 && $6 != 331' $pdbs_folder$str > $pdbs_folder$label.strip.pdb
	pdb4amber -i $pdbs_folder$label.strip.pdb -o $out_folder$label.strip.amber.pdb --dry --nohyd --prot
	wc $out_folder$label.strip.amber.pdb | cat
	rm $pdbs_folder$label.strip.pdb
done

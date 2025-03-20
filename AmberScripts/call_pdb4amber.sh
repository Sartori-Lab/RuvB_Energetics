source /home/victormello/Programs/amber24/amber.sh

pdbs_folder=../data/pdb/RuvB_Monomers/
#pdbs_folder=../data/pdb/RuvB_Trimers/
out_folder=../data/amber/pdb4amber/

for str in $(ls $pdbs_folder)
do
	label=$(echo $str | sed 's/.pdb//')
	pdb4amber -i $pdbs_folder$str -o $out_folder$label.amber.pdb --dry --nohyd --prot
done

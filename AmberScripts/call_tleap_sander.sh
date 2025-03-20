
source /home/victormello/Programs/amber24/amber.sh

#pdbs_folder=../data/pdb/RuvB_Monomers/
pdbs_folder=../data/pdb/RuvB_Trimers/
out_folder=../data/amber/pdb4amber/
energy_folder=../data/amber/sander/

for str in $(ls $pdbs_folder)
do
	label=$(echo $str | sed 's/.pdb//')
	echo $label

	content="\n
	source leaprc.protein.ff19SB\n
	#source leaprc.DNA.OL21\n
	#source leaprc.water.opc\n\n

	#mol = loadpdb $out_folder$label.amber.pdb\n
        mol = loadpdb $out_folder$label.strip.amber.pdb\n

	#addions mol Na+ 0\n
	#solvatebox mol TIP3PBOX 10.0\n

	saveamberparm mol $label.prmtop $label.inpcrd\n

	quit\n"

	echo -e $content > force_field_$label.txt

	tleap -f force_field_$label.txt > $energy_folder$label.tleap.out

	sander -O -i one_cycle.in \
			-o $energy_folder$label.full.out \
			-p $label.prmtop \
			-c $label.inpcrd \
			-inf $energy_folder$label.summary.out

	cat leap.log | grep 'total atoms in file:' | cat
	#tail $energy_folder$label.summary.out -n 5

	rm force_field_$label.txt
	rm $label.prmtop
	rm $label.inpcrd
	rm leap.log
	#rm restrt

done

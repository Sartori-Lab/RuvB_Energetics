
source /home/victormello/Programs/amber24/amber.sh

pdbs_folder=../data/pdb/RuvB_Monomers/
out_folder=../data/amber/pdb4amber/
energy_folder=../data/amber/sander/

for str in $(ls $pdbs_folder)
do
	label=$(echo $str | sed 's/.pdb//')
	echo $label

	content="\n
	source leaprc.protein.ff19SB\n
	#source leaprc.DNA.OL21\n
	#source leaprc.water.opc\n

	mol = loadpdb $out_folder$label.amber.pdb\n

	#solvatebox mol TIP3PBOX 10.0\n
	saveamberparm mol $label.prmtop $label.inpcrd\n

	quit\n"

	echo -e $content > force_field_$label.txt

	tleap -f force_field_$label.txt > $energy_folder$label.tleap.out

	sander -O -i minimize.in \
			-p $label.prmtop \
			-c $label.inpcrd \
			-r $label.rst

	sander -O -i one_cycle.in \
			-o $energy_folder$label.min.full.out \
			-p $label.prmtop \
			-c $label.rst \
			-inf $energy_folder$label.min.summary.out

	rm force_field_$label.txt
	rm $label.prmtop
	rm $label.inpcrd
	rm $label.rst
	#rm restrt

done

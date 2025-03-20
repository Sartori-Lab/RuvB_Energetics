
source /home/victormello/Programs/amber24/amber.sh

hex_folder=../data/pdb/RuvB_Hexamers/
pdbs_folder=../data/pdb/RuvB_Monomers/
out_folder=../data/amber/pdb4amber/
energy_folder=../data/amber/sander/

for str in $(ls $hex_folder)
do
	label=$(echo $str | sed 's/.pdb//')
	echo $label

	content="\n
	source leaprc.protein.ff19SB\n
	#source leaprc.DNA.OL21\n
	source leaprc.water.opc\n\n
	#source leaprc.gaff2\n\n

	chA = loadpdb $out_folder${label}_A.amber.pdb\n
        chB = loadpdb $out_folder${label}_B.amber.pdb\n
        chC = loadpdb $out_folder${label}_C.amber.pdb\n
        chD = loadpdb $out_folder${label}_D.amber.pdb\n
        chE = loadpdb $out_folder${label}_E.amber.pdb\n
        chF = loadpdb $out_folder${label}_F.amber.pdb\n
	hex = combine {chA chB chC chD chE chF}\n\n

	saveamberparm chA ${label}_A.prmtop ${label}_A.inpcrd\n
        saveamberparm chB ${label}_B.prmtop ${label}_B.inpcrd\n
        saveamberparm chC ${label}_C.prmtop ${label}_C.inpcrd\n
        saveamberparm chD ${label}_D.prmtop ${label}_D.inpcrd\n
        saveamberparm chE ${label}_E.prmtop ${label}_E.inpcrd\n
        saveamberparm chF ${label}_F.prmtop ${label}_F.inpcrd\n\n

        solvatebox hex TIP3PBOX 10.0\n\n

        saveamberparm hex ${label}_hex.prmtop ${label}_hex.inpcrd\n


	quit\n"

	echo -e $content > force_field_$label.txt

	#tleap -f force_field_$label.txt > $energy_folder$label.tleap.out

	#sander -O -i one_cycle.in \
	#pmemd -O -i one_cycle.in \
        #               	-o $energy_folder${label}_hex.full.out \
        #               	-p ${label}_hex.prmtop \
        #               	-c ${label}_hex.inpcrd \
       	#               	-inf $energy_folder${label}_hex.summary.out \
	#		-x traj.nc -L 99999999

	#MMPBSA.py -O -i decomp.in \
	#		-o $energy_folder${label}_A.hex.full.out \
	#		-sp ${label}_hex.prmtop \
	#		-cp ${label}_A.prmtop \
	#		-deo $energy_folder${label}_A.res.csv \
	#		-y traj.nc

	#MMPBSA.py -O -i one_cycle.in \
        #                -o $energy_folder${label}_B.hex.full.out \
        #                -sp ${label}_hex.prmtop \
        #                -cp ${label}_B.prmtop

	#MMPBSA.py -O -i one_cycle.in \
        #                -o $energy_folder${label}_C.hex.full.out \
        #                -sp ${label}_hex.prmtop \
        #                -cp ${label}_C.prmtop

	#MMPBSA.py -O -i one_cycle.in \
        #                -o $energy_folder${label}_D.hex.full.out \
        #                -sp ${label}_hex.prmtop \
        #                -cp ${label}_D.prmtop

        #MMPBSA.py -O -i one_cycle.in \
        #                -o $energy_folder${label}_E.hex.full.out \
        #                -sp ${label}_hex.prmtop \
        #                -cp ${label}_E.prmtop

        #MMPBSA.py -O -i one_cycle.in \
        #                -o $energy_folder${label}_F.hex.full.out \
        #                -sp ${label}_hex.prmtop \
        #                -cp ${label}_F.prmtop

	#sander -O -i one_cycle.in \
	#		-o $energy_folder$label.full.out \
	#		-p $label.prmtop \
	#		-c $label.inpcrd \
	#		-inf $energy_folder$label.summary.out

	rm force_field_$label.txt
	rm traj.nc
	rm *.prmtop
	rm *.inpcrd
	rm restrt

done


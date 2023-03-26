#! /bin/tcsh -f
#
#	convert the first F in an MTZ file into a symmetry-expanded text list
#
#		James Holton 6-27-17
#
#
set mtzfile = "$1"
if(! -e "$mtzfile") then
    if("$mtzfile" != "") echo "$mtzfile does not exist"
    cat << EOF
usage: $0 mtzfile.mtz
EOF
    exit 9
endif

if(! $?CCP4) goto try_phenix

mkdir -p ${CCP4_SCR}
set tempfile = ${CCP4_SCR}/tmp_dump$$

set firstF = `mtzdmp $mtzfile | grep -v "H H H " | grep "$2" | awk 'NF>5 && $(NF-1)=="F"{print $NF;exit}'`
set first2Gs = `mtzdmp $mtzfile | grep -v "H H H " | grep "$2" | awk 'NF>5 && $(NF-1)=="G"{print $NF}' | head -n 2`
if("$firstF" == "" && "$first2Gs" == "") then
    echo "ERROR: cannot find any Fs in $mtzfile"
    exit 9
endif

if("$firstF" != "") then

echo "selected $firstF"

cad hklin1 $mtzfile hklout ${tempfile}P1.mtz << EOF > /dev/null
labin file 1 E1=$firstF
outlim space 1
EOF

echo "nref -1\nFORMAT '(3i6,3g35.20)'" |\
mtzdump hklin ${tempfile}P1.mtz |\
awk 'substr($0,1,18)==sprintf("%6d%6d%6d",$1,$2,$3) && $4!="?"{\
      print $1,$2,$3,$4; print -$1,-$2,-$3,$4}' |\
awk '{printf("%4d %4d %4d %s\n",$1,$2,$3,$4)}' |\
cat >! P1.hkl

rm -f ${tempfile}P1.mtz

endif


if("$first2Gs" != "") then

echo "selected $first2Gs"

cad hklin1 $mtzfile hklout ${tempfile}P1.mtz << EOF > /dev/null
labin file 1 E1=$first2Gs[1] E2=$first2Gs[2]
outlim space 1
EOF

echo "nref -1\nFORMAT '(3i6,2g35.20)'" |\
mtzdump hklin ${tempfile}P1.mtz |\
awk 'substr($0,1,18)==sprintf("%6d%6d%6d",$1,$2,$3) && $4!="?" && $5!="?"{\
      print $1,$2,$3,$4; print -$1,-$2,-$3,$5}' |\
awk '{printf("%4d %4d %4d %s\n",$1,$2,$3,$4)}' |\
cat >! P1.hkl

rm -f ${tempfile}P1.mtz

endif

goto finish

tryphenix:

phenix.reflection_file_converter $mtzfile --shelx=shelx --expand_to_p1
awk '{print $1,$2,$3,$4}' shelx.shelx >! P1.hkl
rm -f shelx.shelx


finish:
set count = `cat P1.hkl | wc -l`
echo "$count h k l F lines output to P1.hkl"

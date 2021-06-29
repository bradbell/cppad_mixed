#! /bin/bash -e
cat << EOF > temp.sed
/a1_vector version of ran_likelihood/! b one
: loop_1
N
/}/! b loop_1
d
: one
/a1_vector template_ran_likelihood(/! b two
: loop_2
N
/)/! b loop_2
s|template_ran_likelihood|ran_likelihood|
s|)|) override|
:two
EOF
list=$(ls *.cpp)
for file in $list
do
	git checkout $file
	sed -i $file -f temp.sed
done

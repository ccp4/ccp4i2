#! /bin/bash

export template=$1
export target=$2


mkdir $target
mkdir $target/script

for extension in .def.xml .py _gui.py _report.py
do
ls $template/script/$template$extension
sed s/$template/$target/g $template/script/$template$extension > $target/script/$target$extension
done

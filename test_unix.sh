work=work

if [ ! -d $work ]; then
  mkdir $work
fi

labels=("OPT_uv" "OPT_fluor" "OPT_tadf")

for label in ${labels[@]}; do
  option=${label/_/ }
  cp Samples/ch2o.sdf ${work}/${label}.sdf
  python2 main.py "${work}/${label}.sdf" "${option}" > ${work}/${label}.out
done

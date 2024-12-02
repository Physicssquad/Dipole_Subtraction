output_dir=$(sed -n '7p' run.machine.dat | awk '{print $1}')
echo $output_dir

#------------------------------------------------------------------------------------#
#To view the ungoing process of a parallel code in one go

# fish script
for i in (seq 29 +1 37); cat output_2024_07_25_09_58_$i.Dipole; end

#shell script
for i in {37..29}; do cat "output_2024_07_25_09_58_${i}.Dipole"; done

#------------------------------------------------------------------------------------#

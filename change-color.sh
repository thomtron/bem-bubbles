files=/home/thomas/Documents/ma-thesis/videos/taib-cloud-data/interpolated/mesh-*.ply
output_dir=/home/thomas/Documents/ma-thesis/videos/taib-cloud-data/interpolated-col
mkdir $output_dir

for file in $files 
do 
    ./build/color $file /home/thomas/Documents/ma-thesis/code/first-tries/python-code/exterior-phi/binary.csv $output_dir -0.5 0.5
done
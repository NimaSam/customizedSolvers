# !bin/bash
pic_num=57
inputfile_prefix="V03"
text='U=0.3 m/s'

for num in $(seq 0 $pic_num); do
    # Your code here
    formatted_num=$(printf "%04d" "$num")
    inputfilename="$inputfile_prefix.$formatted_num.jpg"
    outputfilename="$inputfile_prefix.$formatted_num.edited.jpg"
    echo "$inputfilename"
    echo "$outputfilename"
    convert $inputfilename -pointsize 60 -fill black -annotate +900+300 $text $outputfilename
done
convert -loop 9 *.edited.jpg $inputfile_prefix.gif
rm -rf *.edited.jpg


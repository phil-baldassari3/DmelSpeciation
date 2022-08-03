#!/bin/bash

echo "Converting flybase r5 coordinates to r6 coordinates from csv files."
echo "Make sure coordinate_scraper_r5.py python script is ready for multiprocessing."

read -p "Press enter to continue"


cd /Users/philipbaldassari/Desktop/zim-cos_downsampled/null_condition

echo "scraping r5 coordinates from csvs"

python coordinate_scraper_r5.py

echo "done scraping r5 coordinates from csvs"
echo "starting conversion process"

for file in r5*.txt
do
	/Users/philipbaldassari/Desktop/dmel_r5_to_r6/dmel_r5_to_r6_converter.pl --input $file --output r6_$file
done

echo "bookkeeping..."

for file in r5*.txt
do
	rm $file
done

echo "The coordinate_scraper_r6.py python script needs to be prepared for multiprocessing"
echo "Here are the files that this script will be run on:"
echo "\n"

ls r6_r5*.txt

echo "\n"
echo "copy them into the file pool mapping list in the script and then continue this program"
echo "also take this opportunity to check .failed files to check for any major failures in the conversion process before proceeding"
echo "\n"

read -p "Press enter to continue"

echo "scraping r6 coordinates from r5 to r6 conversion files"

python coordinate_scraper_r6.py

echo "done scraping r6 coordinates from conversion files"
echo "bookkeeping..."

for file in *temp.csv
do
	rm $file
done

for file in r6_r5*
do
	rm $file
done


echo "\n"
echo "All processes complete. All flybase r5 coordinates have been converted to r6"











# A small script to rapidly rename WiFeS filename for easy use in PyWiFeS
# Test it before using it ... it will overwrite all files !
#
# To run it : in a terminal (in the location of the data) -> sh PyWiFeS_rename.sh

for i in *; do   
    # Print the command (for test purposes)
    #echo "$i" "->" "`echo $i | sed "s/T2m3w\(.*\)-\(.*\)-\(.*\).fits/\1\3.fits/"`"; 
    
    # Actually do it !
    mv "$i" "`echo $i | sed "s/T2m3w\(.*\)-\(.*\)-\(.*\).fits/\1\3.fits/"`";

done
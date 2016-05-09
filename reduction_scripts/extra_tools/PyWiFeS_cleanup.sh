
# A small script to rapidly all (but p11.fits) files created by PyWiFeS
# Test it before using it ... it will delete without warning !
# Also, make sure you are done with the data reduction ... or you will
# have to restart from scratch !
#
# To run it : in a terminal (in the location fo the data) -> sh PyWiFeS_cleanup.sh

for i in *; do   
    if [[ $i != *"p11.fits" ]]; then
	# Print the command (for test purposes)
	#echo "rm -f $i"

	# Actually do it !
	rm -f $i
    fi
done
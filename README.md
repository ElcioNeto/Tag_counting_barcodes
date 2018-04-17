# Tag counting barcodes
A Program in R for counting tags and identification of samples using barcodes 

## Tutorial

The program works as a Rscript. 

To see the possible entrances of the program, type:

      Rscript script_tag_count_and_sequence_barcode_removal.R --help

The following output will be generate

     Options:
	-f CHARACTER, --fast=CHARACTER
		A joined .fasta file [Required]

	-b CHARACTER, --bar=CHARACTER
		A barcode file [Required]

	-t TAG, --tag=TAG
		write YES if you want to see the header with tag pair [default= NO] 

	-n CHARACTER, --name=CHARACTER
		write YES if you want to see the names with tags in the sequence [default= NO] 

	-i CHARACTER, --inv=CHARACTER
		write YES if you want to invert the orientation of the read (recommended if the processment is taking to long) [default= NO] 

	-h, --help
		Show this help message and exit

The program needs obligatorily of a -f entrance, that is a .fasta file and the -b that is a .txt barcode file.

If only the obligated filleds be fullfilled, the program will produce 2 outputs files.

The options -t and -n are optional. If the user types YES in one (or both) of this options, one (or two) more file(s) will be created

The -i option will only change the methodology, but will not influence the amount of outputs. It may be and alternative when the processment is taking to long.   

### Autors
Daniel K. Morais (https://github.com/kdanielmorais)

Elcio A. O. Neto (https://github.com/ElcioNeto)

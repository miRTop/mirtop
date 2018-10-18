MINTcodes "license-plates"
------------------------------

1. General Information
----------------------
This code was created by Venetia Pliatsika, Isidore Rigoutsos, Jeffery Ma, Phillipe Loher
It can be used to create/encode molecular "license-plates" from sequences and to also decode the "license-plates"
back to sequences.  While initially created for tRFs (tRNA fragments), this tool can be used for 
any genomic sequences including but not limited to:  tRFs, isomiRs, reference miRNA, etc.
For more information on "license-plates", visit https://cm.jefferson.edu/MINTbase and 
refer to publications https://www.ncbi.nlm.nih.gov/pubmed/27153631/ and https://www.ncbi.nlm.nih.gov/pubmed/28220888/.
Contact us at: https://cm.jefferson.edu/contact-us/

2. Terms of Use
---------------
This code can be freely used for research, academic and other non-profit activities.
Only one instance of the code may be used at a time, and then for only one concurrent user. You may not
use the code to conduct any type of application service, service bureau or time-sharing operation or to
provide any remote processing, network processing, network telecommunications or similar services to
any person, entity or organization, whether on a fee basis or otherwise. The code can be copied and
compiled on any platform for the use authorized by these terms and conditions. All copies of the code
must be accompanied by this note. The code cannot be modified without the written permission of the
Computational Medicine Center of Thomas Jefferson University https://cm.jefferson.edu

Commercial use is strictly prohibited. If you wish to use these codes commercially please contact the
Computational Medicine Center of Thomas Jefferson University: https://cm.jefferson.edu/contact-us/

THE CODE IS PROVIDED “AS IS” WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESSED
OR IMPLIED. TO THE FULLEST EXTENT PERMISSIBLE PURSUANT TO APPLICABLE LAW. THOMAS JEFFERSON
UNIVERSITY, AND ITS AFFILIATES, DISCLAIM ALL WARRANTIES, EXPRESS OR IMPLIED, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF TITLE, MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NON-INFRINGEMENT.

NEITHER THOMAS JEFFERSON UNIVERSITY NOR ITS AFFILIATES MAKE ANY REPRESENTATION AS TO THE RESULTS
TO BE OBTAINED FROM USE OF THE CODE.

3. Usage Information
--------------------
Usage (Python 2 and 3 compatible):
	python MINTplates.py example_sequences_to_encode.txt en --p [prefix to add to license plate]
	python MINTplates.py example_license_plates_to_decode.txt de --p [optional prefix, not used in decoding]
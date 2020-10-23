Tool for de Bruijn graph construction
==============

### Version: 0.1

This is a tool for fast de Bruijn graph construction from PacBio HiFi reads.
We use combination of multiple known techniques such as bloom filters, sparse de Bruijn graphs and rolling hashs.
We construct de Bruijn graph of human genome from HiFi dataset with coverage 25 within 4 hours and using 32 threads on any walue of k from 250 to 7000.
For error-corrected reads only 30 minuts is sufficient.

For installation and running instructions please refer to [manual](manual.md)

License
-------

Flye is distributed under a BSD license. See the [LICENSE file](LICENSE) for details.


Credits
-------

This tool is developed by Anton Bankevich in [Pavel Pevzner's lab at UCSD](http://cseweb.ucsd.edu/~ppevzner/)

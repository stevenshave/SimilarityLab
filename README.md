# [SimilarityLab](https://similaritylab.bio.ed.ac.uk)
SimilarityLab; a website for running molecular similarity and target prediction.
---
[https://similaritylab.bio.ed.ac.uk](https://similaritylab.bio.ed.ac.uk)


SimilarityLab currently runs on the University of Edinburgh's Eleanor cloud service.
Styled by bootstrap, served by flask (Python) through Gnunicorn and Nginx, using redis and celery for message queuing and worker threads, SimilarityLab allows querying of the eMolecules database for 3D similars (USRCAT) of a given target protein.  Additionally, 3D similars to known actives in ChEMBL can be used to predict protein targets for small molecules.


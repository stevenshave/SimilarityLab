from pathlib import Path

class ChEMBLIDtoUniprotMapper:
    chemblid_to_uniprotid={}
    def __init__(self, mapfile_location:[Path, str]):
        if not isinstance(mapfile_location, Path):
            mapfile_location=Path(mapfile_location)
        mapfile_location=mapfile_location.expanduser()
        assert mapfile_location.exists(), "Couldnt open chembl to uniprot mapping file. Tried opening"+str(mapfile_location)
        for line in open(str(mapfile_location)).readlines()[1:]:
            if len(line)>2:
                tokens=line.strip().split("\t")
                self.chemblid_to_uniprotid[tokens[1]]=tokens[0]
    def query_chemblid(self,chemblid:str):
        if chemblid in self.chemblid_to_uniprotid.keys():
            return self.chemblid_to_uniprotid[chemblid]
        else:
            return None
import requests
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import fromstring as et_fromstring
import json
from pathlib import Path

class UniProtWebQuerier:
    _uniprotid_to_ecnumber_dict={}
    _tmp_file=None
    def __init__(self, custom_dat_file_name:[Path,str]=Path("dat_uniprot_ecnumbermapping.json")):
        self._tmp_file=custom_dat_file_name
        if isinstance(self._tmp_file, str):
            self._tmp_file=Path(custom_dat_file_name)
        assert isinstance(self._tmp_file, Path), "Cusstom_dat_file_name can be a path, or a string which will be converted to one"
        if self._tmp_file.exists():
            self._uniprotid_to_ecnumber_dict=json.load(open(str(self._tmp_file)))
            print(f"Loaded {len(self._uniprotid_to_ecnumber_dict.keys())} cached uniprot to ecNumber mappings")

    def uniprotid_to_ecnumber(self, uniprotid:str):
        if uniprotid in self._uniprotid_to_ecnumber_dict.keys():
            return self._uniprotid_to_ecnumber_dict[uniprotid]
        else:
            print(f"{uniprotid} not found, querying ...", end="")
            xml_string=requests.get(f"https://www.uniprot.org/uniprot/{uniprotid}.xml").text
            ecnumber_begin_pos=xml_string.find(">",xml_string.find("ecNumber"))+1
            ecnumber_end_pos=xml_string.find("</ecNumber", ecnumber_begin_pos)
            ecnumber=xml_string[ecnumber_begin_pos:ecnumber_end_pos]
            self._uniprotid_to_ecnumber_dict[uniprotid]=str(ecnumber)
            json.dump(self._uniprotid_to_ecnumber_dict, open(str(self._tmp_file),"w"))
            print(f"found : {ecnumber}")
            return ecnumber
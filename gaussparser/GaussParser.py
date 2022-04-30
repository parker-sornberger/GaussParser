"""
author: 
    Parker Sornberger
email:
    pdso223@uky.edu
"""


from types import MappingProxyType
import pandas as pd
import json

class TDDFTParser:
    def __init__(self, file):
        self.filename=file
        with open(file, "r") as file:
            self.contents = file.readlines()
        self.strongest_is_lowest_singlet = "Not yet found"
        self.strongest_is_lowest_triplet = "Not yet found"
        self._get_electrons()
        self._find_states()
        self._decorate_with_orbitals()
    def __repr__(self):
        return F"TDDFTParser at {hex(id(self))} for {self.filename}"
    def __str__(self):
        return F"TDDFTParser for {self.filename}"
    def __len__(self):
        return len(self.all_states)
    
    def __iter__(self):
        return self.all_states.__iter__()
    
    def items(self,kind=None):
        if kind == "s":
            return tuple(list(d.items() for d in self.singlets))
        elif kind == "t":
            return tuple(list(d.items() for d in self.triplets))
        else:
            return tuple(list(d.items() for d in self.all_states))
    def keys(self,kind=None):
        if kind == "s":
            return tuple(list(d.keys() for d in self.singlets))
        elif kind == "t":
            return tuple(list(d.keys() for d in self.triplets))
        else:
            return tuple(list(d.keys() for d in self.all_states))
    def values(self,kind=None):
        if kind == "s":
            return tuple(list(d.values() for d in self.singlets))
        elif kind == "t":
            return tuple(list(d.values() for d in self.triplets))
        else:
            return tuple(list(d.values() for d in self.all_states))
        
    def _get_electrons(self):
        alpha_beta = [line for line in self.contents if "alpha" in line and "beta" in line]
        
        alpha_index = alpha_beta[0].split().index("alpha")
        beta_index = alpha_beta[0].split().index("beta")
        
        alpha = int(alpha_beta[0].split()[alpha_index-1])
        beta = int(alpha_beta[0].split()[beta_index-1])
        if alpha == beta:
            self.homo_lumo_up = F"{alpha} ->{alpha+1}"
            self.homo_lumo_down = F"{alpha} <-{alpha+1}"
        else:
            self.homo_lumo_up = None
            self.homo_lumo_down = None
        
    def _find_states(self):
        self._all = []
        self._singlets = []
        self._triplets = []
        self._state_energy_indices = []
        def clean_lines(lines):
            clean = []
            bad = {"TD", "geom", "SavETr", "Version", "Version="}
            for line in lines:
                
                dont_use = False
                for b in bad:
                    if b in line:
                        dont_use = True
                if not dont_use:
                    clean.append(line)
                    
            return clean
        self.contents = clean_lines(self.contents)
        for i, line in enumerate(self.contents):
            
            if "State" in line:
                #print(line)
                self._state_energy_indices.append(i)
                
                values = line.split()
                name = " ".join(values[:3])
                subd = {}

                subd["name"] = name
                subd["state"] = values[3]
                subd["energy"] = float(values[4])
                subd["wavelength"] = float(values[6])
                subd["oscilator"] = float(values[8].split("=")[-1])
                self._all.append(subd)
        self._singlets = [d for d in self._all if "Singlet" in d["state"]]
        self._triplets = [d for d in self._all if "Triplet" in d["state"]]
        
        
        
    def _decorate_with_orbitals(self): 
        for i, line in enumerate(self.contents):
            if "->" in line or "<-" in line:
                between = self.contents[i-1]
                if "Excited" in between:
                    
                    transitions = []
                    temp_vals = between.split()
                    key = " ".join(temp_vals[:3])
                    
                    place, proxy_d = [(j, d) for (j,d) in enumerate(self._all) if d["name"] == key][0]
                
                td= {}
                
                which = "->" if "->" in line else "<-"
                
                brokenline = [sub.strip().split() for sub in line.split(which)]
                literal = F"{brokenline[0][0]} {which}{brokenline[1][0]}"
                td["name"] = literal
                td["direction"] = "up" if "->" in line else "down"
                td["lower"] = int(brokenline[0][0])
                td["upper"] = int(brokenline[1][0])
                td["strength"] = float(brokenline[1][1])
                
                if "->" in literal:
                    ishomolumo = literal == self.homo_lumo_up
                else:
                    ishomolumo = literal == self.homo_lumo_down
                td["is homo-lumo"] = ishomolumo
                
                transitions.append(MappingProxyType(td))
                proxy_d["transitions"] = tuple(transitions)
                self._all[place] = proxy_d
    
    @property
    def singlets(self):
        return tuple(list(map(MappingProxyType, self._singlets)))
    @property
    def triplets(self):
        return tuple(list(map(MappingProxyType, self._triplets)))
    
    @property
    def all_states(self):
        
        return tuple(list(map(MappingProxyType, self._all)))
    
    @property
    def lowest_singlet(self):
        if not self.singlets:
            return MappingProxyType({})
        energy = self.singlets[0]["energy"]
        oscilator = self.singlets[0]["oscilator"]
        transitions = self.singlets[0]["transitions"]
        return MappingProxyType(dict(name = self.singlets[0]["name"],
                                     energy=energy, oscilator=oscilator, 
                                     transitions=transitions))
    @property
    def lowest_triplet(self):
        if not self.triplets:
            return MappingProxyType({})
        transitions = self.triplets[0]["transitions"]
        return MappingProxyType(dict(name = self.triplets[0]["name"],
                                     energy=self.triplets[0]["energy"], 
                                     transitions = transitions))
    
    
    def homo_lumo_transitions(self):
        hlt = []
        for state in self.all_states:
            trans= state["transitions"]
            has_hlt = set(bool(t["is homo-lumo"]) for t in trans)
            if True in has_hlt:
                hlt.append(state)

        return tuple(hlt)
    
    def find_stongest_oscilator(self):
        which = sorted(self.all_states, key = lambda x : x["oscilator"], reverse=True)[0]
        low_t = self.lowest_triplet
        if low_t:
            if which["name"] == low_t["name"]:
                self.strongest_is_lowest_triplet = True
            else:
                self.strongest_is_lowest_triplet = False
        else:
            self.strongest_is_lowest_triplet = False
        low_s = self.lowest_singlet
        if low_s:
            if which["name"] == low_s["name"]:
                self.strongest_is_lowest_singlet = True
            else:
                self.strongest_is_lowest_singlet = False
        else:
            self.strongest_is_lowest_singlet = False
        return which
    
    
    def _json_dump(self, name, which):
        with open(F"{name}.json", "w") as out:
            json.dump(which, out)
    
    def _text_dump(self, name, which):
        pass
    
    def _pandas_csv_dump(self, name, which):
        #pass which as the item not the name
        index = list(map(lambda x : x["name"], which))
        df = pd.DataFrame(which, index = index).drop(["name"], axis = 1)
        df.to_csv(F"{name}.csv")
        
        
        
    
    def to_file(self, name, fmt, which="all"):
        self._json_dump("k", self.singlets)
        
        

        

class OmegaParser:
    def __init__(self, file):
        self.filename = file
        with open(file, "r") as fn:
            w_data = fn.readlines()[-2].split()[1]
        self.raw_omega = float("0.{}".format(w_data.split('.')[1]))
        self.iop_str = str(int(self.raw_omega* 1e4)).zfill(5) + "00000"
        self.route_param_str = F"iop(3/107={self.iop_str}, 3/108={self.iop_str})"
    
















# GuassParser
Easy TDDFT and Omega tuning parsers

**TDDFTParser**
* Works only with Gaussian TDDFT files. So far, only tested with outputs from Gaussian 16
* All lists are tuples, all dicts are MappingProxyType (it doesn't make sense for anything from this to be editable)
* Can obtain all triplet states, all singlet states, oscilator strenghts, and orbital transitions
* It can find the state which has the strongest transition oscilator strength, and informs if this transition is the lowest energy singlet state

*TDDFTParser Example*
```python
from gaussparser import TDDFTParser

parsed_results = TDDFTParser("MyAwesomeTDDFTFile.log")

print(F"My HOMO-LUMO Up transition is {parsed_results.homo_lumo_up}")
print(F"The state with the highest oscilator strength is {parsed_result.find_strongest_oscilator()}")

print(F" My list of singlet states is {parsed_results.singlets}")

print(F"These are the states that have a HOMO-LUMO or LUMO-HOMO transition {parsed_results.homo_lumo_transitions()}")

```
**TDDFTParser To-do**
* Make file writers (JSON, CSV, and Text)
* Sensible attribute getter : Index based or keyword based?


**OmegaParser**
* Works only with Gaussian log files for omega-tuning and tested only with outputs from Gaussian 16
* Pretty minimalist: You can get the raw omega value, and the key words for geometry optimization using this value

*OmegaParser Example*

```python

from gaussparser import OmegaParser

parsed_results = OmegaParser("MyCoolOmegaTuning.log")

print(F"The as-is omega value is {parsed_results.raw_omega}")

print(F"The keywords I need for a geometry optimization are {parsed_results.route_param_str}")
```

**OmegaParser To-do**
* More universal testing
* __repr__ and __str__ functions? Would anyone actually need this?








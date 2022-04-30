# GuassParser
Easy TDDFT and Omega tuning parsers

**TDDFTParser**
* Works only with Gaussian TDDFT files. So far, only tested with outputs from Gaussian 16



*Example*
```python
from gaussparser import TDDFTParser

parsed_results = TDDFTParser("MyAwesomeTDDFTFile.log")

print(F"My HOMO-LUMO Up transition is {parsed_results.homo_lumo_up}")
print(F"The state with the highest oscilator strength is {parsed_result.find_strongest_oscilator()}")

print(F" My list of singlet states is {parsed_results.singlets}")

print(F"These are the states that have a HOMO-LUMO or LUMO-HOMO transition {parsed_results.homo_lumo_transitions()}")

```




-Set up virtual environment:

```
$python3 -m venv venv
$source venv/bin/activate
$pip install -r requirements.txt
```

-Set up kernel channel in jupyter notebook (demo.ipynb):

```
$pip install ipykernel
$pip install notebook
$python3 -m ipykernel install --user --name=virtual_biopsy
```

-ASAP application: use "measurement" annotation
*Tabspace*+First-click point: starting point, where the biopsy core starts
*Tabspace*+Second-click point: direction point, used to dictate the biopsy direction
*Tabspace*+Number of annotations = number of biopsy cores

-Options: default by
*Tabspace*+Biopsy length: 1.0 mm (std = 0.35)
*Tabspace*+Biopy radius: 0.12 (std = 0.02)
*Tabspace*+They should be x10 but large to compute

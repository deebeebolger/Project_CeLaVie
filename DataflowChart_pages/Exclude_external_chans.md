## Exclude external channels

In the RELAX_cfg file, this parameter is called ElectrodesToDelete. This information should be based on electrodes noted in the cahier de manip (lab notebook). By default, it is blank. In the CLV implementation, the following electrodes are marked for exclusion:
- GSR1, GSR2, Erg1, Erg2, Resp, Plet, Temp, EXG1, EXG2, EXG3, EXG4, EXG5, EXG6, EXG7, EXG8

Those channels to reject are defined in the RELAX config structure as follows :

```
RELAX_cfg.ElectrodesToDelete = {'GSR1', 'GSR2', 'Erg1', 'Erg2', 'Resp', 'Plet', 'Temp', 'EXG1','EXG2', 'EXG3',...
        'EXG4', 'EXG5', 'EXG6', 'EXG7', 'EXG8'};
```


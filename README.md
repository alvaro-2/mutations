# Database of mutations in LLPS-related proteins

This repository contains necessary scripts and files to download, build and update the database.
Installation requirements:
- Python 3.5 or higher
- Pandas
- Snakemake

To run the pipeline to generate files, run

```
cd pipeline
snakemake --cores 8
```
or use the number of cores available in your system

-------------------------------

el pre-processing de cada database esta en c/u de las carpetas (clinvar, disgenet, uniprot).  
Luego, en el main, un script para generar la tabla de mutaciones y otro para la tabla de proteinas.

El orden para ver los scripts seria algo asi:  
- ~~/clinvar/clinvar_box1_mutaciones_en_proteinas_v2.ipynb  ~~
- /disgenet/disgenet_total.ipynb  
- /uniprot/uniprot_variants.ipynb  
- generar_tabla_mutaciones.ipynb  
- generar_tabla_proteinas.ipynb  
    
Por ultimo, se deberia ir guardando las tablas finales para la db en /db_tables

--------------------------

# Database Schema
![esquema](https://github.com/alvaro-2/mutations/blob/main/esquema.png)

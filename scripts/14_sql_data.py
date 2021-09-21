import pandas as pd
import os

DB_NAME = "mlo_mutations"
TABLES_ORDER_NUMTYPE = [
    #{"table": "protein", #for all proteins
    # "file": "protein",
    # "int": ["id_protein", "length"],
    # "float": ["disorder_content"]},
    {"table": "protein", #only the proteins with mutations
     "file": "protein_sel",
     "int": ["id_protein", "length"],
     "float": ["disorder_content"]},  
    {"table": "consequence",
     "file": "consequence",
     "int": ["id_consequence"]}, 
    {"table": "mutation",
     "file": "mutation",
     "int": ["id_mutation", "start_genomic", "end_genomic", "start_aa", "end_aa", 
             "id_protein", "id_consequence"]},     
    {"table": "disease",
     "file": "disease"},  
    #{"table": "disorder_region", #for all proteins
    # "file": "disorder_region",
    # "int": ["id_idr", "start", "end", "length", "id_protein"]},,  
    {"table": "disorder_region",
     "file": "disorder_region_sel", #only the proteins with mutations
     "int": ["id_idr", "start", "end", "length", "id_protein"]},  
    #{"table": "low_complexity",    #for all proteins 
    # "file": "low_complexity",
    # "int": ["id_lc", "start", "end", "length", "id_protein"]}, ,  
    {"table": "low_complexity",
     "file": "low_complexity_sel", #only the proteins with mutations
     "int": ["id_lc", "start", "end", "length", "id_protein"]}, 
    {"table": "pfam_domain",
     "file": "pfam_domain",
     "int": []},  
    {"table": "mutation_has_disease",
     "file": "mutation_has_disease",
     "int": ["id_mutation"]}, 
    {"table": "source",
     "file": "source",
     "int": ["id_source"]},
    {"table": "mlo",
     "file": "mlo",
     "int": ["id_mlo"]},  
    {"table": "rol",
     "file": "rol",
     "int": ["id_rol"]},  
    #{"table": "protein_has_pfam_domain", #for all proteins
    # "file": "protein_has_pfam_domain",
    # "int": ["id_protein", "start", "end", "length"]}, ,  
    {"table": "protein_has_pfam_domain",
     "file": "protein_has_pfam_domain_sel", #only the proteins with mutations
     "int": ["id_protein", "start", "end", "length"]}, 
    {"table": "mutation_has_low_complexity",
     "file": "mutation_has_low_complexity",
     "int": ["id_mutation", "id_lc"]},  
    {"table": "mutation_has_disorder_region",
     "file": "mutation_has_disorder_region",
     "int": ["id_mutation", "id_idr"]}, 
    {"table": "dataset",
     "file": "dataset",
     "int": ["id_dataset"]},     
    #{"table": "protein_has_mlo", #for all proteins
    # "file": "protein_has_mlo",
    # "int": ["id_proteinmlo", "id_protein", "id_dataset", "id_mlo", "id_rol"]},,     
    {"table": "protein_has_mlo",
     "file": "protein_has_mlo_sel", #only the proteins with mutations
     "int": ["id_proteinmlo", "id_protein", "id_dataset", "id_mlo", "id_rol"]},
    {"table": "mutation_has_pfam_domain",
     "file": "mutation_has_pfam_domain",
     "int": ["id_mutation", "id_protein", "start", "end"]}, 
    {"table": "mutation_has_source",
     "file": "mutation_has_source",
     "int": ["id_mutation", "id_source"]},
    {"table": "citation_source",
     "file": "citation_source",
     "int": ["id_citation_source"]},
    {"table": "mutation_has_citation",
     "file": "mutation_has_citation",
     "int": ["id_mutation", "id_citation_source"]},
    {"table": "cross_reference",
     "file": "cross_reference",
     "int": ["id_cross"]},
    {"table": "disease_has_synonyms",
     "file": "disease_has_synonyms",
     "int": ["id_cross"]},
    {"table": "type_ptm",
     "file": "type_ptm",
     "int": ["id_type"]},  
    {"table": "class_ptm",
     "file": "class_ptm",
     "int": ["id_class"]},
    #{"table": "ptm",
    # "file": "ptm", #only the proteins with mutations
    # "int": ["id_ptm", "id_type", "id_class", "pos_aa", "id_protein"]},   
    {"table": "ptm",
     "file": "ptm_sel", #only the proteins with mutations
     "int": ["id_ptm", "id_type", "id_class", "pos_aa", "id_protein"]}, 
    #{"table": "llps_region",
    #"file": "llps_region", #only the proteins with mutations
    #"int": ["id_llps", "start", "end", "length", "id_protein"]},   
    {"table": "llps_regions",
     "file": "llps_region_sel", #only the proteins with mutations
     "int": ["id_llps", "start", "end", "length", "id_protein"]}, 
    {"table": "mutation_has_llps",
     "file": "mutation_has_llps_region", #only the proteins with mutations
     "int": ["id_llps", "id_mutation"]},     
    ]

FOLDER_TABLES = "db_tables"
FOLDER_SQL = "sql"

def generate_all_sql():
    for i in TABLES_ORDER_NUMTYPE:
        generate_sql(i)
        
def generate_all_sources():    
    lines = ["source schema_mlo_mutations.sql;\n"] + ["source " + i['file'] + ".sql;\n" for i in TABLES_ORDER_NUMTYPE]
    with open(os.path.join("..", FOLDER_SQL, "00_schema_and_data_sources.sh"), 'w') as output_source:
        output_source.writelines(lines)
            

def generate_sql(table_info):
    print("Generating " + table_info['table'] + "\n")
    #load the table
    table_data = pd.read_csv(os.path.join("..", FOLDER_TABLES, table_info['file'] + ".tsv"), sep="\t")
    #read the colnames
    coln_names = list(table_data.columns)
        
    for i in coln_names:
        y = table_data[i].isna()
        #convert to string for the next step
        if i in table_info.get('int', []):
            table_data.loc[~y, i] = [str(int(x)) for x in table_data.loc[~y, i]]
        elif i in table_info.get('float', []):
            table_data.loc[~y, i] = [str(float(x)) for x in table_data.loc[~y, i]]
        else:            
            #convert string columns with double-quote
            table_data.loc[~y, i] = ["\"" + str(k) + "\"" for k in table_data[i][~y]] 
        #convert NULLs in string column       
        table_data.loc[y, i] = "NULL"
        
    #convert coln_names to SQL text
    coln_names = (", ").join(["`" + i + "`" for i in coln_names])
    coln_names = "(" + coln_names + ")"
    #generate alter table command
    table_name = "`" + DB_NAME + "`.`" + table_info['table'] + "`"
    with open(os.path.join("..", FOLDER_SQL, table_info['file'] + ".sql"), 'w') as output_sql:
        output_sql.write("LOCK TABLES " + table_name + " WRITE;\n" + "/*!40000 ALTER TABLE " + table_name + " DISABLE KEYS */;\n")
    
    table_data = list(', '.join(x for x in y) for y in table_data.values)
    #generate one insert per 9999 rows  
    for i in range(0, len(table_data), 9999): 
        lines = ["INSERT INTO " + table_name + " " + coln_names, "VALUES"]  
        max_line = min(len(table_data) - 1, i + 9998)
        lines = lines + ["(" + i + ")," for i in table_data[i:max_line]] + ["(" + table_data[max_line] + ");"]
        #add \n at the end of each line
        lines = [s + '\n' for s in lines]
        with open(os.path.join("..", FOLDER_SQL, table_info['file'] + ".sql"), 'a') as output_sql:
            output_sql.writelines(lines)
            
    #close alter table command
    with open(os.path.join("..", FOLDER_SQL, table_info['file'] + ".sql"), 'a') as output_sql:
        output_sql.write("/*!40000 ALTER TABLE " + table_name + " ENABLE KEYS */;\n" + "UNLOCK TABLES;\n")
    
generate_all_sql() 
generate_all_sources()